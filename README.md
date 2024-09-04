# RSBio-Seq

[![Cargo tests](https://github.com/anuradhawick/rsbio-seq/actions/workflows/rust_test.yml/badge.svg)](https://github.com/anuradhawick/rsbio-seq/actions/workflows/rust_test.yml)
[![Downloads](https://static.pepy.tech/badge/rsbio-seq)](https://pepy.tech/project/rsbio-seq)
[![PyPI - Version](https://img.shields.io/pypi/v/rsbio-seq)](https://pypi.org/project/rsbio-seq/)
[![Upload to PyPI](https://github.com/anuradhawick/rsbio-seq/actions/workflows/pypi.yml/badge.svg)](https://github.com/anuradhawick/rsbio-seq/actions/workflows/pypi.yml)

RSBio-Seq intends to provide just reading facility on common sequence formats (FASTA/FASTQ) in both raw and compressed formats.

## Installation 

### 1. From PyPI (Recommended)

Simple use the following command

```bash
pip install rsbio-seq
```

### 2. Build and install from source

To build you need to have the following installed.

- Rust - https://www.rust-lang.org/tools/install
- Maturin - https://www.maturin.rs/installation
- Python environment with Python >=3.9

```bash
maturin develop # this installs the development version in the env
maturin develop --rust # this installs a release version in the env
```

To build a wheel for installation

```bash
maturin build --release
```

You will find the `whl` file inside the `target/wheels` directory. Your `whl` file will have a name depicting your python environment and CPU architecture.


## Usage

Once installed you can import the library and use as follows.

### Reading

```python
from rsbio_seq import SeqReader, SeqWriter, Sequence

# reading
for seq in SeqReader("path/to/seq.fasta.gz"):
    print(seq.id)
    print(seq.seq)
    # for fastq quality line
    print(seq.qual)
    # optional description attribute
    print(seq.desc)
```

### Writing

```python
# writing fasta
seq = Sequence("id", "desc", "ACGT") # id, description, sequence
writer = SeqWriter("out.fasta")
writer.write(seq)
writer.close()

# writing fastq
seq = Sequence("id", "desc", "ACGT") # id, description, sequence
writer = SeqWriter("out.fastq")
writer.write(seq)
writer.close()

# writing gzipped
seq = Sequence("id", "desc", "ACGT", "IIII") # id, description, sequence, quality
writer = SeqWriter("out.fq.gz")
writer.write(seq)
writer.close()
```

Note: `close()` is needed if you want to read within the same program scope. Otherwise, rust will automatically do this for you.

## Authors

- Anuradha Wickramarachchi [https://anuradhawick.com](https://anuradhawick.com)
- Vijini Mallawaarachchi [https://vijinimallawaarachchi.com](https://vijinimallawaarachchi.com)

## Support and contributions

Please get in touch via author websites or GitHub issues. Thanks!
