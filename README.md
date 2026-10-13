# RSBio-Seq

[![Cargo tests](https://github.com/anuradhawick/rsbio-seq/actions/workflows/rust_test.yml/badge.svg)](https://github.com/anuradhawick/rsbio-seq/actions/workflows/rust_test.yml)
[![Downloads](https://static.pepy.tech/badge/rsbio-seq)](https://pepy.tech/project/rsbio-seq)
[![PyPI - Version](https://img.shields.io/pypi/v/rsbio-seq)](https://pypi.org/project/rsbio-seq/)
[![Upload to PyPI](https://github.com/anuradhawick/rsbio-seq/actions/workflows/pypi.yml/badge.svg)](https://github.com/anuradhawick/rsbio-seq/actions/workflows/pypi.yml)
[![License: Apache 2.0 OR GPL v3](https://img.shields.io/badge/License-Apache%202.0%20OR%20GPLv3-blue.svg)](LICENSE)

<div align="center">
<pre>
██████  ███████ ██████  ██  ██████        ███████ ███████  ██████  
██   ██ ██      ██   ██ ██ ██    ██       ██      ██      ██    ██ 
██████  ███████ ██████  ██ ██    ██ █████ ███████ █████   ██    ██ 
██   ██      ██ ██   ██ ██ ██    ██            ██ ██      ██ ▄▄ ██ 
██   ██ ███████ ██████  ██  ██████        ███████ ███████  ██████  
                                                              ▀▀   
</pre>
</div>

RSBio-Seq reads and writes FASTA and FASTQ files in plain text (`fasta`, `fa`, `fna`, `fastq`, `fq`) and compressed formats (`.gz`). Indexed reading and writing are supported for both FASTA and FASTQ, including BGZF-compressed files.

The sequential reader uses a buffered Rust parser in `src/io.rs`, inspired by
[klib's kseq.h](https://github.com/attractivechaos/klib/blob/master/kseq.h),
with reusable byte buffers and support for wrapped FASTA/FASTQ records, LF/CRLF
line endings, and concatenated gzip members (including BGZF). FASTQ quality
length must match sequence length; malformed input raises `OSError` during
iteration. The library does not depend on the Rust `bio` crate.

## Installation

### 1. From PyPI (Recommended)

Use the following command to install from PyPI.

```bash
pip install rsbio-seq
```

### 2. Build and install from source

To build from source, make sure you have the following programs installed.

- Rust - https://www.rust-lang.org/tools/install
- Maturin - https://www.maturin.rs/installation
- Python environment with Python >=3.9 - https://www.python.org/downloads/

To build and install the development version of the wheel.

```bash
maturin develop # this installs the development version in the env
maturin develop --release # this installs a release version in the env
```

To build a release mode wheel for installation, use this command.

```bash
maturin build --release
```

You will find the `whl` file inside the `target/wheels` directory. Your `whl` file will have a name depicting your python environment and CPU architecture. The built wheel can be installed using this command.

```bash
pip install target/wheels/*.whl
```

## Usage

Once installed you can import the library and use as follows.

### Reading

```python
from rsbio_seq import SeqReader, Sequence, ascii_to_phred

# each seq entry is of type Sequence
seq: Sequence

for seq in SeqReader("path/to/seq.fasta.gz"):
    print(seq.id)
    print(seq.seq)
    # for fastq quality line
    print(seq.qual) # prints IIII
    print(ascii_to_phred(seq.qual)) # prints [40, 40, 40, 40]
    # optional description attribute
    print(seq.desc)
```

### Reading from FASTA or FASTQ (`fai`/`fai+gzi`) indexes

The indexed reader retrieves FASTA and FASTQ records by ID from plain text or BGZF-compressed files. Both formats require a `.fai` index; FASTQ indexes also contain quality offsets.

```python
from rsbio_seq import SeqReaderIndexed, Sequence

seqs = SeqReaderIndexed(
    "path/to/seq.fa",
    "path/to/seq.fa.fai"
)
seq: Sequence = seqs["Record_1"]
print(seq.id)
print(seq.seq)
print(seq.desc)

"Record_2" in seqs # returns a boolean
```

For FASTQ, use the same API with a `.fq` or `.fastq` file and its index. The returned record includes its quality string.

```python
seqs = SeqReaderIndexed("path/to/seq.fq", "path/to/seq.fq.fai")
seq = seqs["Read_1"]
print(seq.seq)
print(seq.qual)
```

BGZF-compressed FASTA and FASTQ files also require a `.gzi` index. Indexed access requires BGZF compression rather than ordinary gzip.

```python
seqs = SeqReaderIndexed(
    "path/to/seq.fa.gz",
    "path/to/seq.fa.gz.fai",
    "path/to/seq.fa.gz.gzi",
)
seq: Sequence = seqs["Record_1"]
print(seq.id)
print(seq.seq)
print(seq.desc)

"Record_2" in seqs # returns a boolean

# BGZF-compressed FASTQ
reads = SeqReaderIndexed(
    "path/to/seq.fq.gz",
    "path/to/seq.fq.gz.fai",
    "path/to/seq.fq.gz.gzi",
)
read = reads["Read_1"]
print(read.seq)
print(read.qual)
```

Using an invalid key will result in `KeyError`. 

### Writing

```python
from rsbio_seq import SeqWriter, Sequence, phred_to_ascii

# writing fasta
seq = Sequence("id", "desc", "ACGT") # id, description, sequence
writer = SeqWriter("out.fasta")
writer.write(seq)
writer.close()

# writing fastq
seq = Sequence("id", "desc", "ACGT", "IIII") # id, description, sequence, quality
writer = SeqWriter("out.fastq")
writer.write(seq)
writer.close()

# writing gzipped
seq = Sequence("id", "desc", "ACGT", "IIII") # id, description, sequence, quality
writer = SeqWriter("out.fq.gz")
writer.write(seq)
writer.close()

# writing gzipped with phred score translation
qual = phred_to_ascii([40, 40, 40, 40])
seq = Sequence("id", "desc", "ACGT", qual) # id, description, sequence, quality
writer = SeqWriter("out.fq.gz")
writer.write(seq)
writer.close()
```

Note: `close()` is only required if you want to read the file again in the same function/code scope. Closing opened files is a good practice either way.

We provide two utility functions for your convenience.

* `phred_to_ascii` - convert phred scores list of numbers to a string
* `ascii_to_phred` - convert the quality string to a list of numbers

RSBio-Seq reads and write quality string in ascii format only. Please use these helper functions to translate if you intend to read them.

### Writing FASTA or FASTQ with indexes (`fai`/`fai+gzi`)

Set `index=True` to write a FASTA or FASTQ file with a `.fai` index. A `.gz` suffix selects BGZF compression and also creates a `.gzi` index. Index paths are inferred from the output filename. Call `close()` to finish writing the file and indexes before reading them.

```python
# Plain text
seq = Sequence("id", "desc", "ACGT") # id, description, sequence
writer = SeqWriter("out.fasta", index=True)
writer.write(seq)
writer.close()

# Compressed
seq = Sequence("id", "desc", "ACGT") # id, description, sequence
writer = SeqWriter("out.fa.gz", index=True)
writer.write(seq)
writer.close()
```

FASTQ indexes include offsets for both sequences and quality strings.

```python
from rsbio_seq import SeqWriter, SeqReaderIndexed, Sequence

seq = Sequence("Read_1", "desc", "ACGT", "IIII")

# Creates out.fq and out.fq.fai
writer = SeqWriter("out.fq", index=True)
writer.write(seq)
writer.close()

# Creates out.fq.gz, out.fq.gz.fai, and out.fq.gz.gzi
writer = SeqWriter("out.fq.gz", index=True)
writer.write(seq)
writer.close()

reads = SeqReaderIndexed("out.fq.gz", "out.fq.gz.fai", "out.fq.gz.gzi")
read = reads["Read_1"]
print(read.seq)   # ACGT
print(read.qual)  # IIII
```

For wrapped output, pass a line width of at least 10 to `writer.write(seq, wrap=60)`. FASTQ sequences and quality strings are wrapped at the same width.

## Authors

- Anuradha Wickramarachchi [https://anuradhawick.com](https://anuradhawick.com)
- Vijini Mallawaarachchi [https://vijinimallawaarachchi.com](https://vijinimallawaarachchi.com)

## Support and contributions

Please get in touch via author websites or GitHub issues. Thanks!

## License

RSBio-Seq is dual licensed under the [Apache License, Version 2.0](LICENSE-APACHE)
or the [GNU General Public License, Version 3 only](LICENSE-GPL), at your option.

SPDX license expression: `Apache-2.0 OR GPL-3.0-only`.
