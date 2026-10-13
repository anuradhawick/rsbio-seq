from pathlib import Path
import shutil

import pytest

from rsbio_seq import SeqReader, SeqReaderIndexed, build_index


@pytest.mark.parametrize("suffix", ["fa", "fq", "fa.gz", "fq.gz"])
def test_build_existing_index(tmp_path, suffix):
    source = Path("test_data") / f"reads_indexable.{suffix}"
    path = tmp_path / source.name
    shutil.copyfile(source, path)
    assert build_index(str(path)) is None
    assert Path(f"{path}.fai").read_bytes() == Path(f"{source}.fai").read_bytes()
    gzi = f"{path}.gzi" if suffix.endswith(".gz") else None
    if gzi:
        assert Path(gzi).read_bytes() == Path(f"{source}.gzi").read_bytes()
    reader = SeqReaderIndexed(str(path), f"{path}.fai", gzi)
    for expected in reversed(list(SeqReader(str(path)))):
        actual = reader[expected.id]
        assert (actual.id, actual.desc, actual.seq, actual.qual) == (
            expected.id, expected.desc, expected.seq, expected.qual
        )


def test_invalid_input_preserves_index(tmp_path):
    path = tmp_path / "bad.fq"
    path.write_bytes(b"@x\nAA\n+\n!\n")
    index = Path(f"{path}.fai")
    index.write_bytes(b"existing index")
    with pytest.raises(OSError):
        build_index(str(path))
    assert index.read_bytes() == b"existing index"
