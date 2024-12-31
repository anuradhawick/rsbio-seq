import gzip
import pathlib

from rsbio_seq import Sequence, SeqWriter

dir = pathlib.Path(__file__).parent


def test_write_fa():
    writer = SeqWriter(dir.joinpath("../test_data/out.fa").as_posix())
    seq = Sequence(
        "Record_1",
        "Desc",
        "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA",
    )
    writer.write(seq)
    writer.close()

    assert (
        open(dir.joinpath("../test_data/out.fa").as_posix()).read()
        == """>Record_1 Desc
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA
"""
    )


def test_write_fa_index():
    writer = SeqWriter(dir.joinpath("../test_data/out_ix.fa").as_posix(), index=True)
    seq = Sequence(
        "Record_1",
        "Desc 1",
        "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA",
    )
    writer.write(seq, 20)
    seq = Sequence(
        "Record_2",
        "Desc 2",
        "TTAACAACTTAAGGGTTTTCAAATAGAGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCC",
    )
    writer.write(seq, 20)
    writer.close()

    assert (
        open(dir.joinpath("../test_data/out_ix.fa.fai").as_posix()).read()
        == """Record_1	72	17	20	21
Record_2	72	110	20	21
"""
    )


def test_write_fa_gz_index():
    writer = SeqWriter(dir.joinpath("../test_data/out_ix.fa.gz").as_posix(), index=True)
    seq = Sequence(
        "Record_1",
        "Desc 1",
        "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA",
    )
    writer.write(seq, 20)
    seq = Sequence(
        "Record_2",
        "Desc 2",
        "TTAACAACTTAAGGGTTTTCAAATAGAGGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCC",
    )
    writer.write(seq, 20)
    writer.close()

    assert (
        open(dir.joinpath("../test_data/out_ix.fa.gz.fai").as_posix()).read()
        == """Record_1	72	17	20	21
Record_2	72	110	20	21
"""
    )
    assert dir.joinpath("../test_data/out_ix.fa.gz.gzi").exists()


def test_write_fa_gz():
    writer = SeqWriter(dir.joinpath("../test_data/out.fa.gz").as_posix())
    seq = Sequence(
        "Record_1",
        "Desc",
        "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA",
    )
    writer.write(seq)
    writer.close()

    assert (
        gzip.open(dir.joinpath("../test_data/out.fa.gz").as_posix(), "rt").read()
        == """>Record_1 Desc
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA
"""
    )


def test_write_fq():
    writer = SeqWriter(dir.joinpath("../test_data/out.fq").as_posix())
    seq = Sequence(
        "Record_1",
        "Desc",
        "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA",
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGII>IIIII-I)8I",
    )
    writer.write(seq)
    writer.close()

    assert (
        open(dir.joinpath("../test_data/out.fq").as_posix()).read()
        == """@Record_1 Desc
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGII>IIIII-I)8I
"""
    )


def test_write_fq_gz():
    writer = SeqWriter(dir.joinpath("../test_data/out.fq.gz").as_posix())
    seq = Sequence(
        "Record_1",
        "Desc",
        "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA",
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGII>IIIII-I)8I",
    )
    writer.write(seq)
    writer.close()

    assert (
        gzip.open(dir.joinpath("../test_data/out.fq.gz").as_posix(), "rt").read()
        == """@Record_1 Desc
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBIIIIIIIIIIIIIIIIIIIIIIIGII>IIIII-I)8I
"""
    )
