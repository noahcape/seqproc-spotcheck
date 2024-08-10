# seqproc-spotcheck
Spotchecking single-cell sequencing methods

### 10x3v3 Barcodes
```
wget -P 10x3v3/data wget https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz
gunzip 10x3v3/data/3M-february-2018.txt.gz > 10x-barcodes.txt
```

Execute command similar to:
```
index=0
OUT_FILE=./10x3v3/data/10x-barcodes-metadata.tsv
echo "sample_id\tbarcode" >> $OUT_FILE

while read line; do
  echo "$index\t$line" >> $OUT_FILE;

  (( index = index + 1 ))
done <./10x3v3/data/10x-barcodes.txt
```

To convert `.tsv` file of barcodes to input type for `fgbio`'s `fqtk` processing method.

### 10x3V3 Data
```
wget -P 10x3v3/data -c \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/009/SRR10587809/SRR10587809_1.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/009/SRR10587809/SRR10587809_2.fastq.gz
```

### sci-RNA-seq3 Data
```
wget -P sci-rna-seq3/data \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR782/006/SRR7827206/SRR7827206_1.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR782/006/SRR7827206/SRR7827206_2.fastq.gz
```

### SPLiTseq Data
```
wget -P splitseq/data -c \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR675/002/SRR6750042/SRR6750042_1.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR675/002/SRR6750042/SRR6750042_2.fastq.gz
```

### Usage

By default `gtime -v` is set on and should print usage the stdout.

```
$ cargo run -- -h
Process single-cell sequencing protocols with specific programs

Usage: seqproc_paper_benchmarking [OPTIONS] --r1 <R1> --r2 <R2>

Options:
  -p, --program <PROGRAM>    Program to process sequencing reads 1=seqproc 2=splitcode 3=fgbio 4=flexiplex [default: 1]
  -g, --protocol <PROTOCOL>  Protocol which sequencing reads derive from 1=10x3v3 2=sci-RNA-seq3 3=SPLiTseq [default: 1]
  -1, --r1 <R1>              path to r1 file (.fastq/.fasta)
  -2, --r2 <R2>              path to r2 file (.fastq/.fasta)
  -h, --help                 Print help
  -V, --version              Print version
```