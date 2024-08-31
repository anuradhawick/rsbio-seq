from typing import Iterator

class Sequence:
    """
    Sequence entry.

    Attributes
        id (str): record id of the sequence.
        desc (str): description of the record if available.
        seq (str): content of the sequence.
        qual (str): quality string of the sequence if available.
    """

    id: str
    desc: str
    seq: str
    qual: str

class SeqIO:
    """
    Sequence reader.
    """

    def __init__(self, path: str) -> None:
        """
        Initialise the reader with path of the file.

        Args:
            path (str): The path to file fasta, fastq, fa, fq and compressed formats with gz are supported.
        """
        ...

    def __iter__(self) -> Iterator[Sequence]:
        """
        Return an iterator of Sequence objects.

        Returns:
            Iterator[Sequence]: An iterator over sequences in the file.
        """
        ...

__all__ = ["Sequence", "SeqIO"]
