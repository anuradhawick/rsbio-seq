import pathlib

from rsbio_seq import SeqReader, Sequence

dir = pathlib.Path(__file__).parent


def test_read_fa():
    seqs = SeqReader(dir.joinpath("../test_data/reads.fa").as_posix())
    seq: Sequence = next(seqs)
    assert seq.id == "Record_1"
    assert (
        seq.seq
        == "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA"
    )
    seq = next(seqs)
    assert seq.id == "Record_2"
    assert (
        seq.seq
        == "GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT"
    )


def test_read_fa_gz():
    seqs = SeqReader(dir.joinpath("../test_data/reads.fa.gz").as_posix())
    seq: Sequence = next(seqs)
    assert seq.id == "Record_1"
    assert (
        seq.seq
        == "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA"
    )
    seq: Sequence = next(seqs)
    assert seq.id == "Record_2"
    assert (
        seq.seq
        == "GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT"
    )


def test_read_fq():
    seqs = SeqReader(dir.joinpath("../test_data/reads.fq").as_posix())
    seq: Sequence = next(seqs)
    assert seq.id == "Read_1"
    assert (
        seq.seq
        == "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA"
    )
    assert (
        seq.qual
        == "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIIIIIII>IIIIII/"
    )
    seq: Sequence = next(seqs)
    assert seq.id == "Read_2"
    assert (
        seq.seq
        == "GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT"
    )
    assert (
        seq.qual
        == "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGII>IIIII-I)8I"
    )


def test_read_fq_gz():
    seqs = SeqReader(dir.joinpath("../test_data/reads.fq.gz").as_posix())
    seq: Sequence = next(seqs)
    assert seq.id == "Read_1"
    assert (
        seq.seq
        == "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA"
    )
    assert (
        seq.qual
        == "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9ICIIIIIIIIIIIIIIIIIIIIDIIIIIII>IIIIII/"
    )
    seq: Sequence = next(seqs)
    assert seq.id == "Read_2"
    assert (
        seq.seq
        == "GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT"
    )
    assert (
        seq.qual
        == "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGII>IIIII-I)8I"
    )
