use clap::Parser;
use std::{
    os::unix::process::CommandExt,
    process::{Child, Command},
};

static SPLITCODE_SCI_RNA_SEQ_CONFIG: &str = "./sci-rna-seq3/sci-rna-seq3-config.txt";
static FQTK_DEMUX_10X_METADATA: &str = "./10x3v3/10x-barcodes-metadata.tsv";
static SPLITCODE_SPLITSEQ_SUB_CONFIG: &str = "./splitseq/splitseq-rt-bc.txt";

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
}

enum Protocol {
    TenX,
    SciRNASeq3,
    SPLiTseq,
}

impl Protocol {
    fn from(num: usize) -> Self {
        match num {
            1 => Self::TenX,
            2 => Self::SciRNASeq3,
            3 => Self::SPLiTseq,
            _ => panic!("Unsupported protocol"),
        }
    }
}

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

    fn exec(
        &self,
        protocol: Protocol,
        r1: std::path::PathBuf,
        r2: std::path::PathBuf,
    ) -> Child {
        let mut command = Command::new("gtime");
        command.arg("-v");

        match self {
            Program::Seqproc => todo!(),
            Program::Splitcode => match protocol {
                Protocol::TenX => command
                    .arg("splitcode")
                    .arg("-x 0:0<10x-barcode-umi>0:16,0:16<10x-barcode-umi>0:28")
                    .arg("--x-only")
                    .arg("-nFastqs=2")
                    .arg(r1)
                    .arg(r2)
                    .spawn()
                    .expect("Failed splitcode processing on 10x data"),
                Protocol::SciRNASeq3 => command
                    .arg("splitcode")
                    .arg("-c")
                    .arg(SPLITCODE_SCI_RNA_SEQ_CONFIG)
                    .arg("-nFastqs=2")
                    .arg("-n 30000")
                    .arg("--x-only")
                    .arg(r1)
                    .arg(r2)
                    .spawn()
                    .expect("Failed splitcode processing on sci-RNA-seq3 data"),
                Protocol::SPLiTseq => command
                    .arg("splitcode")
                    .arg("-c")
                    .arg(SPLITCODE_SPLITSEQ_SUB_CONFIG)
                    .arg("-nFastqs=1")
                    .arg("-n 4")
                    .arg("--x-only")
                    .arg(r2)
                    .spawn()
                    .expect("Failed splitcode processing on SPLiTseq data bc map"),
            },
            Program::Fgbio => match protocol {
                Protocol::TenX => command
                    .args(["fqtk", "demux"])
                    .arg("--inputs")
                    .arg(r1)
                    .arg(r2)
                    .args(["--read-structures", "16B12M", "+T"])
                    .args(["--sample-metadata", FQTK_DEMUX_10X_METADATA])
                    .args(["--output", "out"])
                    .spawn()
                    .expect("Failed fqtk processing on 10x data"),
                Protocol::SciRNASeq3 => {
                    panic!("fgbio fqtk cannot support the sci-rna-seq3 protocol")
                }
                Protocol::SPLiTseq => panic!("fgbio fqtk cannot support the SPLiTseq protocol"),
            },
            Program::Flexiplex => todo!(),
        }
    }
}

fn main() {
    let args = Cli::parse();

    let program = Program::from(args.program);
    let protocol = Protocol::from(args.protocol);

    if !program.protocol_supported(&protocol) {
        panic!("Program doesn't support that protocol.")
    }

    let mut output = program.exec(protocol, args.r1, args.r2);
    let _ = output.wait();

    clean_up();
}

fn fgbio_10x() {
    let r1 = std::path::PathBuf::from("10x3v3/data/10x-r1.fastq");
    let r2 = std::path::PathBuf::from("10x3v3/data/10x-r2.fastq");

    let mut output = Program::Fgbio.exec(Protocol::TenX, r1, r2);
    let _ = output.wait();

    clean_up();
}

fn splitcode_10x() {
    let r1 = std::path::PathBuf::from("10x3v3/data/10x-r1.fastq");
    let r2 = std::path::PathBuf::from("10x3v3/data/10x-r2.fastq");

    let mut output = Program::Splitcode.exec(Protocol::TenX, r1, r2);
    let _ = output.wait();

    clean_up();
}

fn splitcode_sci_rna_seq3() {
    let r1 = std::path::PathBuf::from("sci-rna-seq3/data/SRR7827206_1_head.fastq");
    let r2 = std::path::PathBuf::from("sci-rna-seq3/data/SRR7827206_2_head.fastq");

    let mut output = Program::Splitcode.exec(Protocol::SciRNASeq3, r1, r2);
    let _ = output.wait();

    clean_up();
}

fn clean_up() {
    if std::path::Path::new("out").exists() {
        Command::new("rm").arg("-r").arg("out").exec();
    }
}
