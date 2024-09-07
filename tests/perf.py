import gzip

from Bio import Seq
from Bio import SeqIO as BioSeqIO
from Bio import SeqRecord
from rsbio_seq import SeqReader, Sequence, SeqWriter, phred_to_ascii
from tqdm import tqdm


# Generates a specified number of Sequence objects with a repeated "ACGT" sequence.
def test_seqs_rs(count: int):
    seq = "ACGT" * 1_000_000
    for c in range(count):
        yield Sequence(f"rec_id_{c+1}", f"description_{c+1}", seq)


# Generates a specified number of SeqRecord objects with a repeated "ACGT" sequence.
def test_seqs_bio(count: int):
    seq = "ACGT" * 1_000_000
    for c in range(count):
        yield SeqRecord.SeqRecord(
            id=f"rec_id_{c+1}", seq=Seq.Seq(seq), description=f"description_{c+1}"
        )


# Generates a specified number of Sequence objects with a repeated "ACGT" sequence and quality scores.
def test_fq_seqs_rs(count: int):
    seq = "ACGT" * 1_000_000
    qual_phred = [40, 40, 40, 40] * 1_000_000
    for c in range(count):
        qual = phred_to_ascii(qual_phred)
        yield Sequence(f"rec_id_{c+1}", f"description_{c+1}", seq, qual)


# Generates a specified number of SeqRecord objects with a repeated "ACGT" sequence and quality scores.
def test_fq_seqs_bio(count: int):
    seq = "ACGT" * 1_000_000
    qual = [40, 40, 40, 40] * 1_000_000
    for c in range(count):
        yield SeqRecord.SeqRecord(
            id=f"rec_id_{c+1}",
            seq=Seq.Seq(seq),
            letter_annotations={"phred_quality": qual},
            description=f"description_{c+1}",
        )


if __name__ == "__main__":
    # writing
    print("Writing normal fasta files (wrapped)")
    writer = SeqWriter("out_rs.fasta")
    for record in tqdm(test_seqs_rs(100), desc="Rs-writer"):
        writer.write(record, 60)

    BioSeqIO.write(
        tqdm(test_seqs_bio(100), desc="Bio-writer"), "out_bio.fasta", "fasta"
    )

    print("Writing 2-line fasta files")
    writer = SeqWriter("out_rs.fasta")
    for record in tqdm(test_seqs_rs(100), desc="Rs-writer"):
        writer.write(record)

    BioSeqIO.write(
        tqdm(test_seqs_bio(100), desc="Bio-writer"), "out_bio.fasta", "fasta-2line"
    )

    print("Writing fastq files")
    writer = SeqWriter("out_rs.fastq")
    for record in tqdm(test_fq_seqs_rs(100), desc="Rs-fq-writer"):
        writer.write(record)

    BioSeqIO.write(
        tqdm(test_fq_seqs_bio(100), desc="Bio-fq-writer"), "out_bio.fastq", "fastq"
    )

    print("Writing fastq gzipped files")
    writer = SeqWriter("out_rs.fq.gz")
    for record in tqdm(test_fq_seqs_rs(100), desc="Rs-fq-writer"):
        writer.write(record)

    with gzip.open("out_bio.fq.gz", "wt") as handle:
        BioSeqIO.write(
            tqdm(test_fq_seqs_bio(100), desc="Bio-fq-writer"),
            handle,
            "fastq",
        )

    # reading
    print("Reading fasta files")
    for seq in tqdm(
        SeqReader("out_bio.fasta"),
        desc="Rs-iteration",
    ):
        pass

    for seq in tqdm(
        BioSeqIO.parse("out_bio.fasta", "fasta"),
        desc="Bio-iteration",
    ):
        pass

    print("Reading fastq gzipped files")
    for seq in tqdm(
        SeqReader("out_bio.fq.gz"),
        desc="Rs-gz-iteration",
    ):
        pass

    with gzip.open("out_bio.fq.gz", "rt") as handle:
        for seq in tqdm(
            BioSeqIO.parse(handle, "fastq"),
            desc="Bio-gz-iteration",
        ):
            pass


# Writing normal fasta files (wrapped)
# Rs-writer:  100it [00:00, 332.09it/s]
# Bio-writer: 100it [00:01, 97.90it/s]

# Writing 2-line fasta files
# Rs-writer:  100it [00:00, 522.58it/s]
# Bio-writer: 100it [00:00, 547.41it/s]

# Writing fastq files
# Rs-fq-writer:  100it [00:08, 11.14it/s]
# Bio-fq-writer: 100it [00:12,  7.88it/s]

# Writing fastq gzipped files
# Rs-fq-writer:  100it [00:10,  9.99it/s]
# Bio-fq-writer: 100it [00:14,  6.76it/s]

# Reading fasta files
# Rs-iteration:  100it [00:00, 353.14it/s]
# Bio-iteration: 100it [00:00, 206.77it/s]

# Reading fastq gzipped files
# Rs-gz-iteration:  100it [00:01, 75.07it/s]
# Bio-gz-iteration: 100it [00:06, 14.39it/s]
