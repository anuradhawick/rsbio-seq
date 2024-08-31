# RSBio-Seq

[![Cargo tests](https://github.com/anuradhawick/rsbio-seq/actions/workflows/rust_test.yml/badge.svg)](https://github.com/anuradhawick/rsbio-seq/actions/workflows/rust_test.yml)

RSBio intends to provide just reading facility on common sequence formats (FASTA/FASTQ) in both raw and compressed formats.

## Build and install from source

To build you need to have the following installed.

- Rust - https://www.rust-lang.org/tools/install
- Maturin - https://www.maturin.rs/installation
- Python environment with Python >=3.9

```bash
maturin develop # this installs the development version in the env
maturin develop --rust # this installs a release version in the env
```

To build a wheel

```bash
maturin build --release
```

You will find the `whl` file inside the `target/wheels` directory.

## Install from PyPI

Simple use the following command

```bash
pip install rsbio-seq
```

## Usage

Once installed you can import the library and use as follows.

```python
from rsbio_seq import SeqIO

for seq in SeqIO("path/to/seq.fasta.gz"):
    print(seq.id)
    print(seq.seq)
    # for fastq quality line
    print(seq.qual)
    # optional description attribute
    print(seq.desc)
```

## Timeline

- [x] Reading fasta/fastq formats raw and gzipped.
- [ ] Writing fasta/fastq formats raw and gzipped.

## Authors

- Anuradha Wickramarachchi [https://anuradhawick.com](https://anuradhawick.com)
- Vijini Mallawaarachchi [https://vijinimallawaarachchi.com](https://vijinimallawaarachchi.com)

## Support and contributions

Please get in touch via author websites or GitHub issues. Thanks!
