use bio::io::fastq::{self, FastqRead};
use clap::Parser;
use read_structure::ReadStructure;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::time::Instant;
use std::{
    fs::File, io::Error, os::unix::process::CommandExt, path::PathBuf, process::Command,
    str::FromStr,
};

/// CHANGE THIS
static SEQPROC_PATH: &str = "../seqproc/target/release/seqproc";

static SPLITCODE_SCI_RNA_SEQ_CONFIG: &str = "./sci-rna-seq3/sci-rna-seq3-config.txt";
static SPLITCODE_SPLITSEQ_SUB_CONFIG: &str = "./splitseq/splitseq-rt-bc.txt";
static SEQPROC_10X_EFGDL: &str = "./10x3v3/10x3v3.efgdl";
static SEQPROC_SCI_RNA_SEQ3_EFGDL: &str = "./sci-rna-seq3/sci-rna-seq3.efgdl";
static SEQPROC_SPLITSEQ_EFGDL: &str = "./splitseq/splitseq.efgdl";
static SEQPROC_SPLITSEQ_MAPPING: &str = "./splitseq/bc_mapping.tsv";

static MAX_NUM_READS: &str = "800000000";

// splitcode + sci-rna-seq3
// cargo run -- --program 2 --protocol 2 --r1 ./sci-rna-seq3/data/SRR7827206_1_head.fastq --r2 ./sci-rna-seq3/data/SRR7827206_2_head.fastq

// seqproc + sci-rna-seq3
// cargo run -- --program 1 --protocol 2 --r1 ./sci-rna-seq3/data/SRR7827206_1_head.fastq --r2 ./sci-rna-seq3/data/SRR7827206_2_head.fastq

/// Process single-cell sequencing protocols with specific programs.
#[derive(Parser, Debug)]
#[command(version, about)]
struct Cli {
    /// Program to process sequencing reads
    /// 1=seqproc
    /// 2=splitcode
    /// 3=fgbio
    /// 4=flexiplex
    #[arg(short = 'p', long = "program", default_value_t = 1)]
    program: usize,

    /// Protocol which sequencing reads derive from
    /// 1=10x3v3
    /// 2=sci-RNA-seq3
    /// 3=SPLiTseq
    #[arg(short = 'g', long = "protocol", default_value_t = 1)]
    protocol: usize,

    /// path to r1 file (.fastq/.fasta)
    #[arg(short = '1', long = "r1")]
    r1: std::path::PathBuf,
    /// path to r2 file (.fastq/.fasta)
    #[arg(short = '2', long = "r2")]
    r2: std::path::PathBuf,

    /// number of times to repeat processing
    #[arg(short = 'r', long = "repeat", default_value_t = 1)]
    repeat: usize,
}

enum Protocol {
    TenX,
    SciRNASeq3,
    SPLiTseq,
}

#[allow(dead_code)]
impl Protocol {
    fn from(num: usize) -> Self {
        match num {
            1 => Self::TenX,
            2 => Self::SciRNASeq3,
            3 => Self::SPLiTseq,
            _ => panic!("Unsupported protocol"),
        }
    }

    fn name(&self) -> &str {
        match self {
            Protocol::TenX => "10x3v3",
            Protocol::SciRNASeq3 => "sci-rna-seq3",
            Protocol::SPLiTseq => "splitseq",
        }
    }
}

#[derive(Clone, Debug)]
enum Program {
    Seqproc,
    Splitcode,
    Fgbio,
    Flexiplex,
}

impl Program {
    fn from(num: usize) -> Self {
        match num {
            1 => Self::Seqproc,
            2 => Self::Splitcode,
            3 => Self::Fgbio,
            4 => Self::Flexiplex,
            _ => panic!("Unsupport program"),
        }
    }

    fn protocol_supported(&self, protocol: &Protocol) -> bool {
        match self {
            Program::Seqproc => matches!(
                protocol,
                Protocol::SPLiTseq | Protocol::SciRNASeq3 | Protocol::TenX
            ),
            Program::Splitcode => matches!(
                protocol,
                Protocol::SPLiTseq | Protocol::SciRNASeq3 | Protocol::TenX
            ),
            Program::Fgbio => matches!(protocol, Protocol::TenX),
            Program::Flexiplex => matches!(protocol, Protocol::TenX),
        }
    }

    fn out_file(&self, protocol: &Protocol) -> &str {
        match &self {
            Program::Seqproc => match protocol {
                Protocol::TenX => "./10x3v3/seqproc_out/seqproc_10x3v3.fastq",
                Protocol::SciRNASeq3 => "sci-rna-seq3/seqproc_out/seqproc_sci_rna_seq3.fastq",
                Protocol::SPLiTseq => "splitseq/seqproc_out/seqproc_splitseq.fastq",
            },
            Program::Splitcode => match protocol {
                Protocol::TenX => "splitcode_10x3v3.fastq",
                Protocol::SciRNASeq3 => "splitcode_sci_rna_seq3.fastq",
                Protocol::SPLiTseq => "splitcode_splitseq.fastq",
            },
            Program::Fgbio => match protocol {
                Protocol::TenX => "./10x3v3/read_structures_out/read_structures_10x3v3.fastq",
                Protocol::SciRNASeq3 => todo!(),
                Protocol::SPLiTseq => todo!(),
            },
            Program::Flexiplex => match protocol {
                Protocol::TenX => todo!(),
                Protocol::SciRNASeq3 => todo!(),
                Protocol::SPLiTseq => todo!(),
            },
        }
    }

    fn time_file(&self, protocol: &Protocol) -> &str {
        match self {
            Program::Seqproc => match protocol {
                Protocol::TenX => "./10x3v3/seqproc_out/out.txt",
                Protocol::SciRNASeq3 => "./sci-rna-seq3/seqproc_out/out.txt",
                Protocol::SPLiTseq => "./splitseq/seqproc_out/out.txt",
            },
            Program::Splitcode => match protocol {
                Protocol::TenX => "./10x3v3/splitcode_out/out.txt",
                Protocol::SciRNASeq3 => "./sci-rna-seq3/splitcode_out/out.txt",
                Protocol::SPLiTseq => "./splitseq/splitcode_out/out.txt",
            },
            Program::Fgbio => match protocol {
                Protocol::TenX => "./10x3v3/read_structures_out/out.txt",
                Protocol::SciRNASeq3 => todo!(),
                Protocol::SPLiTseq => todo!(),
            },
            Program::Flexiplex => match protocol {
                Protocol::TenX => todo!(),
                Protocol::SciRNASeq3 => todo!(),
                Protocol::SPLiTseq => todo!(),
            },
        }
    }

    fn time(&self, protocol: &Protocol) {
        let file_path = self.time_file(protocol);

        if !PathBuf::from(file_path).exists() {
            return;
        }

        // x.xxuser x.xxsystem x:xx.xelsapsed x%CPU (xavgtext+xavgdata xmaxresident)k
        // xinputs+xoutputs (xmajor+xminor)pagefaults xswaps
        let file_contents =
            fs::read_to_string(file_path).expect("Failed to read executation statistics");

        let contents = file_contents
            .split('\n')
            .flat_map(|s| s.split(' '))
            .collect::<Vec<_>>();

        let wall_time = contents[2]
            .as_bytes()
            .iter()
            .filter(|c| !c.is_ascii_alphabetic())
            .copied()
            .collect::<Vec<_>>();

        if contents.is_empty() || wall_time.is_empty() {
            return;
        }

        let page_size = contents[5]
            .as_bytes()
            .iter()
            .filter(|c| !c.is_ascii_alphabetic() && **c != b')')
            .copied()
            .collect::<Vec<_>>();

        let mut alt_path = PathBuf::from(file_path);
        alt_path.set_file_name("final_timed.txt");

        let mut out = OpenOptions::new()
            .create(true)
            .append(true)
            .open(alt_path)
            .expect("Unable to open final time file.");

        out.write_all(&wall_time).expect("Failed to write");
        out.write_all(&[b'\n']).expect("Failed to write");
        out.write_all(&page_size).expect("Failed to write");
        out.write_all(&[b'\n']).expect("Failed to write");
    }

    #[allow(dead_code)]
    fn clean_up(&self, protocol: &Protocol) -> Result<(), Error> {
        std::fs::remove_file(self.out_file(protocol))
    }

    #[allow(dead_code)]
    fn parse_record_id<'a>(&self, id: &'a str) -> &'a str {
        match self {
            Program::Seqproc => id.split(' ').collect::<Vec<_>>().first().unwrap(),
            Program::Splitcode => id,
            Program::Fgbio => id,
            Program::Flexiplex => id,
        }
    }

    fn exec(&self, protocol: &Protocol, r1: &std::path::PathBuf, r2: &std::path::PathBuf) {
        let mut command = Command::new("gtime");
        command.args(["-o", self.time_file(protocol)]);

        match self {
            Program::Seqproc => match protocol {
                Protocol::TenX => {
                    let _ = command
                        .arg(SEQPROC_PATH)
                        .args(["-g", SEQPROC_10X_EFGDL])
                        .arg("-1")
                        .arg(r1)
                        .arg("-2")
                        .arg(r2)
                        .args(["-o", self.out_file(protocol)])
                        .args(["-t", "6"])
                        .spawn()
                        .expect("Failed seqproc processing on 10x data")
                        .wait();
                }
                Protocol::SciRNASeq3 => {
                    let _ = command
                        .arg(SEQPROC_PATH)
                        .args(["-g", SEQPROC_SCI_RNA_SEQ3_EFGDL])
                        .arg("-1")
                        .arg(r1)
                        .arg("-2")
                        .arg(r2)
                        .args(["-o", self.out_file(protocol)])
                        .args(["-t", "6"])
                        .spawn()
                        .expect("Failed seqproc processing on sci-rna-seq3 data")
                        .wait();
                }
                Protocol::SPLiTseq => {
                    let _ = command
                        .arg(SEQPROC_PATH)
                        .args(["-g", SEQPROC_SPLITSEQ_EFGDL])
                        .arg("-1")
                        .arg(r1)
                        .arg("-2")
                        .arg(r2)
                        .args(["-o", self.out_file(protocol)])
                        .args(["-t", "6"])
                        .args(["-a", SEQPROC_SPLITSEQ_MAPPING])
                        .spawn()
                        .expect("Failed seqproc processing on splitseq data")
                        .wait();
                }
            },
            Program::Splitcode => match protocol {
                Protocol::TenX => {
                    let _ = command
                        .arg("splitcode")
                        .args(["-x", "0:0<splitcode_10x3v3>0:16,0:16<splitcode_10x3v3>0:28"])
                        .arg("--x-only")
                        .arg("-nFastqs=2")
                        .args(["-n", MAX_NUM_READS])
                        .args(["-t", "6"])
                        .arg(r1)
                        .arg(r2)
                        .spawn()
                        .expect("Failed splitcode processing on 10x data")
                        .wait();
                }
                Protocol::SciRNASeq3 => {
                    let _ = command
                        .arg("splitcode")
                        .arg("-c")
                        .arg(SPLITCODE_SCI_RNA_SEQ_CONFIG)
                        .arg("--x-only")
                        .arg("-nFastqs=2")
                        .args(["-n", MAX_NUM_READS])
                        .args(["-t", "6"])
                        .arg(r1)
                        .arg(r2)
                        .spawn()
                        .expect("Failed splitcode processing on sci-RNA-seq3 data")
                        .wait();
                }
                Protocol::SPLiTseq => {
                    let _ = command
                        .arg("splitcode")
                        .arg("-c")
                        .arg(SPLITCODE_SPLITSEQ_SUB_CONFIG)
                        .arg("-nFastqs=2")
                        .args(["-n", MAX_NUM_READS])
                        .arg("--x-only")
                        .args(["-t", "6"])
                        .arg(r1)
                        .arg(r2)
                        .spawn()
                        .expect("Failed splitcode processing on SPLiTseq data bc map")
                        .wait();
                }
            },
            Program::Fgbio => match protocol {
                Protocol::TenX => {
                    let now = Instant::now();

                    {
                        let left_rs = ReadStructure::from_str("16B12M").unwrap();
                        let right_rs = ReadStructure::from_str("+T").unwrap();

                        read_validate_write_records(r1, self.out_file(protocol), left_rs);
                        read_validate_write_records(r2, "/dev/null", right_rs);
                    }

                    let mut time_file = PathBuf::from(self.time_file(protocol));
                    time_file.set_file_name("out.txt");

                    let mut out = OpenOptions::new()
                        .append(true)
                        .create(true)
                        .open(time_file)
                        .expect("Failed to open file");

                    let elapsed = now.elapsed();
                    out.write_all(format!("{:?}\n", elapsed.as_millis()).as_bytes())
                        .expect("Failed to write");
                }
                Protocol::SciRNASeq3 => {
                    panic!("fgbio fqtk cannot support the sci-rna-seq3 protocol")
                }
                Protocol::SPLiTseq => panic!("fgbio fqtk cannot support the SPLiTseq protocol"),
            },
            Program::Flexiplex => todo!(),
        }

        match self {
            Program::Splitcode | Program::Seqproc => {
                println!("{:?}", self);
                self.time(protocol);
            }
            _ => (),
        };
    }
}

fn main() {
    let args = Cli::parse();

    let program = Program::from(args.program);
    let protocol = Protocol::from(args.protocol);

    if !program.protocol_supported(&protocol) {
        panic!("Program doesn't support that protocol.")
    }

    for _ in 0..=args.repeat {
        program.exec(&protocol, &args.r1, &args.r2);
    }

    clean_up();
}

fn clean_up() {
    if std::path::Path::new("out").exists() {
        Command::new("rm").arg("-r").arg("out").exec();
    }
}

#[cfg(test)]
mod seqproc_tests {
    use super::*;
    use bio::io::fastq::Record;
    use std::io::Write;

    fn as_path_buf(_str: &str) -> PathBuf {
        std::path::PathBuf::from(_str)
    }

    #[test]
    fn read_structures() {
        let r1 = "./10x3v3/data/SRR10587809_1_head.fastq";
        let r2 = "./10x3v3/data/SRR10587809_2_head.fastq";

        compare_programs(Program::Seqproc, Program::Fgbio, Protocol::TenX, r1, r2)
    }

    #[test]
    fn sci_rna_seq3() {
        let r1 = "./sci-rna-seq3/data/SRR7827206_1_head.fastq";
        let r2 = "./sci-rna-seq3/data/SRR7827206_2_head.fastq";

        compare_programs(
            Program::Seqproc,
            Program::Splitcode,
            Protocol::SciRNASeq3,
            r1,
            r2,
        )
    }

    #[test]
    fn tenx() {
        let r1 = "./10x3v3/data/SRR10587809_1_head.fastq";
        let r2 = "./10x3v3/data/SRR10587809_2_head.fastq";

        compare_programs(Program::Seqproc, Program::Splitcode, Protocol::TenX, r1, r2)
    }

    #[test]
    fn splitseq() {
        let r1 = "./splitseq/data/SRR6750042_1_head.fastq";
        let r2 = "./splitseq/data/SRR6750042_2_head.fastq";

        compare_programs(
            Program::Seqproc,
            Program::Splitcode,
            Protocol::SPLiTseq,
            r1,
            r2,
        )
    }

    fn compare_programs(prog1: Program, prog2: Program, protocol: Protocol, r1: &str, r2: &str) {
        let r1 = as_path_buf(r1);
        let r2 = as_path_buf(r2);

        let _ = prog1.exec(&protocol, &r1, &r2);
        let outf1 = prog1.out_file(&protocol);

        let _ = prog2.exec(&protocol, &r1, &r2);
        let outf2 = prog2.out_file(&protocol);

        compare_fastq(outf1, outf2);
    }

    fn compare_fastq(fastq1: &str, fastq2: &str) {
        let fastq1_file = File::open(fastq1).expect("Could not open fastq");
        let fastq2_file = File::open(fastq2).expect("Could not open fastq");

        let fastq1_records = records(fastq1_file);
        let fastq2_records = records(fastq2_file);

        assert_eq!(fastq1_records.len(), fastq2_records.len())
    }

    fn records(fastq: File) -> Vec<Record> {
        // create FASTQ reader
        let mut reader = fastq::Reader::new(fastq);
        let mut record = fastq::Record::new();
        let mut records: Vec<Record> = vec![];

        reader.read(&mut record).expect("Failed to parse record");
        while !record.is_empty() {
            let check = record.check();
            if check.is_ok() {
                records.push(record.clone());
            }

            let mut res = reader.read(&mut record);
            while res.is_err() {
                res = reader.read(&mut record);
            }
            res.unwrap();
        }

        records
    }

    fn write_fastqs(outf: &mut File, records: Vec<Record>) -> Result<(), Error> {
        for record in records {
            let qual = record.qual();
            let id = record.id();
            let seq = record.seq();
            let res = write!(
                outf,
                "{id}\n{s}\n+\n{q}\n",
                s = String::from_utf8(seq.to_vec()).unwrap(),
                q = String::from_utf8(qual.to_vec()).unwrap()
            );

            if res.is_err() {
                return Err(res.err().unwrap());
            }
        }

        Ok(())
    }
}

fn read_validate_write_records(from: &PathBuf, to: &str, rs: ReadStructure) {
    let mut reader = fastq::Reader::new(File::open(from).expect("Could not open fastq"));
    let mut record = fastq::Record::new();
    let mut writer = fastq::Writer::new(File::create(to).expect("Could not open fastq"));

    reader.read(&mut record).expect("Failed to parse record");
    while !record.is_empty() {
        if let Ok(()) = record.check() {
            let segments = rs.segments();
            if segments.len() == 2 {
                writer
                    .write_record(&record)
                    .expect("Failed to write record")
            };
        };

        reader
            .read(&mut record)
            .expect("Failed to read next record");
    }
}
