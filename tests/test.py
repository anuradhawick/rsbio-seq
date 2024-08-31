from rsbio_seq import SeqIO, Sequence
import pathlib

dir = pathlib.Path(__file__).parent


def test_fa():
    seqs = SeqIO(dir.joinpath("../test_data/reads.fa").as_posix())
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


def test_fa_gz():
    seqs = SeqIO(dir.joinpath("../test_data/reads.fa.gz").as_posix())
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


def test_fq():
    seqs = SeqIO(dir.joinpath("../test_data/reads.fq").as_posix())
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


def test_fq_gz():
    seqs = SeqIO(dir.joinpath("../test_data/reads.fq.gz").as_posix())
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
