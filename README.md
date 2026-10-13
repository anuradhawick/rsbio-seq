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

RSBio-Seq reads and writes FASTA (`.fasta`, `.fa`, `.fna`) and FASTQ
(`.fastq`, `.fq`) files. Append `.gz` for compressed files. File format is
inferred from these lowercase filename extensions; `.bgz` and `.bgzf` suffixes
are not supported.

| Operation | Plain text | Ordinary gzip (`.gz`) | BGZF (`.gz`) |
| --- | --- | --- | --- |
| Sequential reading | Supported | Supported | Supported |
| Indexed reading | Requires `.fai` | Unsupported | Requires `.fai` and `.gzi` |
| Building indexes for existing files | Creates `.fai` | Unsupported | Creates `.fai` and `.gzi` |
| Writing | Supported | Use `index=False` (default) | Use `index=True`; also creates `.fai` and `.gzi` |

The sequential reader uses a buffered Rust parser in `src/io.rs`, inspired by
[klib's kseq.h](https://github.com/attractivechaos/klib/blob/master/kseq.h),
with direct reads into record buffers that become owned strings, and support
for wrapped FASTA/FASTQ records, LF/CRLF
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

From a checkout of this repository, create and activate a virtual environment,
then install Maturin and build the extension into that environment:

```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install "maturin>=1.9.3,<2.0"
maturin develop
# Optional optimized build:
maturin develop --release
```

To build a release mode wheel for installation, use this command.

```bash
maturin build --release
```

The wheel is written to `target/wheels`. Install the wheel matching your platform
and Python version (replace the placeholder with the actual filename):

```bash
python -m pip install "target/wheels/<wheel-filename>.whl"
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
    # FASTA records have an empty quality string.
    print(seq.qual)  # prints an empty line
    # optional description attribute
    print(seq.desc)
```

For FASTQ, `qual` contains the quality string. For example, `"IIII"` decodes
to four Phred scores of 40:

```python
from rsbio_seq import SeqReader, ascii_to_phred

for seq in SeqReader("path/to/seq.fastq.gz"):
    print(seq.qual)
    print(ascii_to_phred(seq.qual))

assert ascii_to_phred("IIII") == [40, 40, 40, 40]
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

A missing record ID raises `KeyError`. Errors opening the file or loading its
indexes raise `OSError`.

### Indexing existing files

```python
from rsbio_seq import build_index, SeqReaderIndexed

build_index("reads.fa")  # creates reads.fa.fai
build_index("reads.fq.gz")  # BGZF input; creates .fai and .gzi
reads = SeqReaderIndexed("reads.fq.gz", "reads.fq.gz.fai", "reads.fq.gz.gzi")
read = reads["Read_1"]
```

`build_index(path)` indexes an existing file without rewriting it. Compressed
input must be BGZF; ordinary gzip raises `OSError`. Existing indexes are replaced
after input validation. Record IDs are the first whitespace-delimited token
following `>` or `@` in the header. Indexing requires unique IDs, nonempty
sequences, no blank lines, and sequence lines without whitespace. Each record
must use consistent line widths and LF or CRLF endings, except that the final
line may be shorter or have no newline. FASTQ quality length and wrapping must
match the sequence. Invalid or truncated input raises `OSError`. Rebuild the
indexes if the input file changes.

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

Call `writer.close()` after the last write to finalize the output and any
indexes, especially the `.gzi` index for BGZF output. Close a writer only once
and do not write to it afterward. Creating a writer replaces any existing file
at the output path.

For FASTQ, supply a quality string with exactly one character per sequence base;
the writer does not validate these lengths. Use ASCII sequence and quality
strings without embedded line breaks, and IDs without whitespace. For indexed
output, use unique IDs and nonempty sequences. FASTA output ignores `qual`.

We provide two utility functions for your convenience.

* `phred_to_ascii` - convert phred scores list of numbers to a string
* `ascii_to_phred` - convert the quality string to a list of numbers

Quality strings use ASCII characters. These helpers convert between quality
strings and numerical scores using Phred+33 encoding. For standard printable
FASTQ qualities, use scores from 0 to 93.

### Writing FASTA or FASTQ with indexes (`fai`/`fai+gzi`)

Set `index=True` to write a FASTA or FASTQ file with a `.fai` index. With
`index=True`, a `.gz` suffix selects BGZF compression and also creates a `.gzi`
index. With `index=False`, `.gz` output uses ordinary gzip. Index paths are
inferred from the output filename. Call `close()` to finish writing the file and
indexes before reading them.

```python
from rsbio_seq import SeqWriter, Sequence

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
