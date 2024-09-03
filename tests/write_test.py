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
