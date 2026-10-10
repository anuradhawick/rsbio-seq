# Indexed FASTQ fixtures

`reads_indexable.fq` contains the same seven record IDs, descriptions, and
sequences as `reads_indexable.fa`. Each record uses the four-line FASTQ format.
Quality strings repeat `!5I@+` (Phred+33 scores 0, 20, 40, 31, and 10), truncated
to the sequence length.

`reads_indexable.fq.gz` is BGZF-compressed, with multiple blocks to exercise
random access across block boundaries. The `.fai` files use the six-column
FASTQ index format, including the quality offset; `.gzi` indexes BGZF blocks.

To regenerate the indexes and compressed copy after changing the plain fixture:

```sh
samtools fqidx test_data/reads_indexable.fq
bgzip --force --keep --index test_data/reads_indexable.fq
samtools fqidx test_data/reads_indexable.fq.gz
```

These fixtures prepare for FASTQ indexed-reader support, which is not yet
implemented. The current sequential `SeqReader` uses a single-member gzip
decoder and cannot read this multi-block BGZF fixture completely.
