# RSBio-Seq

[![Cargo tests](https://github.com/anuradhawick/rsbio-seq/actions/workflows/rust_test.yml/badge.svg)](https://github.com/anuradhawick/rsbio-seq/actions/workflows/rust_test.yml)
[![Downloads](https://static.pepy.tech/badge/rsbio-seq)](https://pepy.tech/project/rsbio-seq)
[![PyPI - Version](https://img.shields.io/pypi/v/rsbio-seq)](https://pypi.org/project/rsbio-seq/)
[![Upload to PyPI](https://github.com/anuradhawick/rsbio-seq/actions/workflows/pypi.yml/badge.svg)](https://github.com/anuradhawick/rsbio-seq/actions/workflows/pypi.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

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

RSBio-Seq intends to provide reading/writing facility on common sequence formats (FASTA/FASTQ) in both raw (`fasta`, `fa`, `fna`, `fastq`, `fq`) and compressed formats (`.gz`).

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
maturin develop --rust # this installs a release version in the env
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

## Authors

- Anuradha Wickramarachchi [https://anuradhawick.com](https://anuradhawick.com)
- Vijini Mallawaarachchi [https://vijinimallawaarachchi.com](https://vijinimallawaarachchi.com)

## Support and contributions

Please get in touch via author websites or GitHub issues. Thanks!
