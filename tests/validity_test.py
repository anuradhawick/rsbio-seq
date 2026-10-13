"""Reader validity against Biopython, plus malformed-input and round-trip tests.

Requires pytest and biopython. Run against a built/installed RSBio-Seq module:
    python -m pytest tests/validity_test.py -q

Biopython includes the identifier in record.description and stores FASTQ
qualities as Phred integers; normalize those representations before comparison.
Comparisons use valid ASCII sequence data, not permissive malformed-input rules.
"""
import gzip
import random
import zlib
from pathlib import Path

import pytest
from Bio import SeqIO

from rsbio_seq import SeqReader, Sequence, SeqWriter


def read_rsbio(path):
    return [(r.id, r.desc, r.seq, r.qual) for r in SeqReader(str(path))]


def read_biopython(path, fmt):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="ascii") as handle:
        result = []
        for record in SeqIO.parse(handle, fmt):
            header = record.description.split(maxsplit=1)
            desc = header[1] if len(header) == 2 else ""
            qual = "".join(chr(q + 33) for q in record.letter_annotations.get("phred_quality", []))
            result.append((record.id, desc, str(record.seq), qual))
        return result


def save_input(tmp_path, fmt, payload, compressed=False):
    path = tmp_path / ("reads." + fmt + (".gz" if compressed else ""))
    path.write_bytes(gzip.compress(payload) if compressed else payload)
    return path


def wrap(value, width):
    return [value[i:i + width] for i in range(0, len(value), width)] or [""]


@pytest.mark.parametrize("fmt", ["fasta", "fastq"])
@pytest.mark.parametrize("compressed", [False, True], ids=["plain", "gzip"])
@pytest.mark.parametrize("newline", ["\n", "\r\n"], ids=["LF", "CRLF"])
@pytest.mark.parametrize("final_newline", [False, True], ids=["no-final-LF", "final-LF"])
def test_generated_records_match_biopython(tmp_path, fmt, compressed, newline, final_newline):
    # Deterministic, independently known records spanning the 64 KiB read buffer.
    rng = random.Random(20261013)
    expected = []
    lines = []
    lengths = [0, 1, 2, 79, 255, 256, 257, 65535, 65536, 65537, 150000, 3]
    for index, length in enumerate(lengths):
        ident = "read_" + str(index)
        desc = "description with spaces " + str(index) if index % 2 else ""
        seq = "".join(rng.choices("ACGTNacgtnRYKMSWBDHV", k=length))
        qual = "".join(rng.choices("!5I@+>~", k=length)) if fmt == "fastq" else ""
        expected.append((ident, desc, seq, qual))
        title = ident + (" " + desc if desc else "")
        lines.append((">" if fmt == "fasta" else "@") + title)
        # Different wrap widths exercise qualities whose line lengths differ
        # from sequence lines, as permitted by FASTQ.
        lines.extend(wrap(seq, [1, 71, 65537][index % 3]))
        if fmt == "fastq":
            lines.append("+" + (title if index % 2 else ""))
            lines.extend(wrap(qual, [13, 97, 70001][index % 3]))
    payload = (newline.join(lines) + (newline if final_newline else "")).encode("ascii")
    path = save_input(tmp_path, fmt, payload, compressed)
    assert read_biopython(path, fmt) == expected
    assert read_rsbio(path) == expected


@pytest.mark.parametrize("fmt", ["fasta", "fastq"])
@pytest.mark.parametrize("length", [65535, 65536, 65537, 200000])
def test_unwrapped_long_header_sequence_and_quality(tmp_path, fmt, length):
    desc = "x" * length
    seq = ("ACGT" * ((length + 3) // 4))[:length]
    qual = ("@+!~" * ((length + 3) // 4))[:length] if fmt == "fastq" else ""
    payload = (">" if fmt == "fasta" else "@") + "long " + desc + "\n" + seq
    if fmt == "fastq":
        payload += "\n+\n" + qual
    path = save_input(tmp_path, fmt, payload.encode("ascii"))
    expected = [("long", desc, seq, qual)]
    assert read_biopython(path, fmt) == expected
    assert read_rsbio(path) == expected


@pytest.mark.parametrize("fmt", ["fasta", "fastq"])
@pytest.mark.parametrize("compressed", [False, True])
def test_empty_file_and_repeated_eof(tmp_path, fmt, compressed):
    path = save_input(tmp_path, fmt, b"", compressed)
    assert read_biopython(path, fmt) == []
    reader = SeqReader(str(path))
    assert iter(reader) is reader
    for _ in range(3):
        with pytest.raises(StopIteration):
            next(reader)


@pytest.mark.parametrize("fmt,payload,expected", [
    ("fasta", b"\n>one\tdescription with spaces\r\nAC\r\n\r\nGT\r\n>empty\r\n>last\nN", [
        ("one", "description with spaces", "ACGT", ""), ("empty", "", "", ""), ("last", "", "N", "")]),
    ("fastq", b"@one\tdescription\r\nAC\r\nGT\r\n+\r\n@+\r\n!I\r\n@two\nN\n+\n#", [
        ("one", "description", "ACGT", "@+!I"), ("two", "", "N", "#")]),
    ("fasta", b">duplicate\nA\n>duplicate\nT\n", [
        ("duplicate", "", "A", ""), ("duplicate", "", "T", "")]),
    ("fastq", b"@duplicate\nA\n+\n!\n@duplicate\nT\n+\n~\n", [
        ("duplicate", "", "A", "!"), ("duplicate", "", "T", "~")]),
])
def test_special_records_match_biopython(tmp_path, fmt, payload, expected):
    path = save_input(tmp_path, fmt, payload)
    assert read_biopython(path, fmt) == expected
    assert read_rsbio(path) == expected


@pytest.mark.parametrize("payload", [
    b"bad\n", b"@one\nAC\n", b"@one\nAC\n@next\n",
    b"@one\nAC\n+", b"@one\nAC\n+\n", b"@one\nAC\n+\n!",
    b"@one\nA\n+\n!!\n", b"@one\nACGT\n+\n!!\n!\n",
], ids=["bad-header", "missing-plus", "next-header-before-plus", "unfinished-plus",
        "missing-quality", "short-quality", "long-quality", "short-wrapped-quality"])
@pytest.mark.parametrize("compressed", [False, True])
def test_malformed_fastq_rejected_by_both(tmp_path, payload, compressed):
    path = save_input(tmp_path, "fastq", payload, compressed)
    with pytest.raises(ValueError):
        read_biopython(path, "fastq")
    reader = SeqReader(str(path))
    with pytest.raises(OSError):
        next(reader)
    for _ in range(2):
        with pytest.raises(StopIteration):
            next(reader)


@pytest.mark.parametrize("fmt,payload", [("fasta", b">\nAC\n"), ("fastq", b"@\nA\n+\n!\n")])
def test_empty_identifier_is_rejected(tmp_path, fmt, payload):
    # RSBio-Seq intentionally validates identifiers more strictly than Biopython.
    path = save_input(tmp_path, fmt, payload)
    with pytest.raises(OSError, match="identifier is empty"):
        next(SeqReader(str(path)))


def test_valid_record_before_malformed_record(tmp_path):
    path = save_input(tmp_path, "fastq", b"@good\nAC\n+\n!I\n@bad\nACGT\n+\n!!\n")
    reader = SeqReader(str(path))
    first = next(reader)
    assert (first.id, first.seq, first.qual) == ("good", "AC", "!I")
    with pytest.raises(OSError, match="truncated"):
        next(reader)
    assert (first.id, first.seq, first.qual) == ("good", "AC", "!I")
    with pytest.raises(StopIteration):
        next(reader)


@pytest.mark.parametrize("fmt", ["fasta", "fastq"])
def test_gzip_members_split_inside_records(tmp_path, fmt):
    payload = b">one desc\nACGT\n>two\nN\n" if fmt == "fasta" else b"@one desc\nACGT\n+\n@+!I\n@two\nN\n+\n~\n"
    path = tmp_path / ("members." + fmt + ".gz")
    # A gzip member boundary need not coincide with a line or record boundary.
    path.write_bytes(b"".join(gzip.compress(payload[i:i + 3]) for i in range(0, len(payload), 3)))
    expected = [("one", "desc", "ACGT", "@+!I" if fmt == "fastq" else ""),
                ("two", "", "N", "~" if fmt == "fastq" else "")]
    assert read_biopython(path, fmt) == expected
    assert read_rsbio(path) == expected


@pytest.mark.parametrize("damage", ["truncated", "checksum", "bad-header"])
@pytest.mark.parametrize("fmt", ["fasta", "fastq"])
def test_damaged_gzip_raises_io_error(tmp_path, damage, fmt):
    payload = b">one\nACGT\n" if fmt == "fasta" else b"@one\nACGT\n+\nIIII\n"
    data = bytearray(gzip.compress(payload))
    if damage == "truncated":
        data = data[:-5]
    elif damage == "checksum":
        data[-8] ^= 0xFF
    else:
        data[0] = 0
    path = tmp_path / ("broken." + fmt + ".gz")
    path.write_bytes(data)
    with pytest.raises((OSError, EOFError, zlib.error)):
        read_biopython(path, fmt)
    with pytest.raises(OSError):
        read_rsbio(path)


@pytest.mark.parametrize("fmt", ["fasta", "fastq"])
@pytest.mark.parametrize("storage", ["plain", "gzip", "bgzf"])
@pytest.mark.parametrize("width", [None, 10, 79])
def test_writer_round_trip_matches_biopython(tmp_path, fmt, storage, width):
    expected = [("long", "description", "ACGTN" * 20000, "@+!I~" * 20000 if fmt == "fastq" else ""),
                ("short", "", "N", "!" if fmt == "fastq" else "")]
    path = tmp_path / ("roundtrip." + fmt + (".gz" if storage != "plain" else ""))
    writer = SeqWriter(str(path), index=storage == "bgzf")
    try:
        for record in expected:
            writer.write(Sequence(*record), wrap=width)
    finally:
        writer.close()
    assert read_biopython(path, fmt) == expected
    assert read_rsbio(path) == expected


@pytest.mark.parametrize("fmt,suffix", [("fasta", "fa"), ("fastq", "fq")])
@pytest.mark.parametrize("compressed", [False, True])
def test_existing_indexable_fixtures_match_biopython(fmt, suffix, compressed):
    path = Path(__file__).resolve().parents[1] / "test_data" / ("reads_indexable." + suffix + (".gz" if compressed else ""))
    records = read_rsbio(path)
    assert len(records) == 7
    assert records == read_biopython(path, fmt)
