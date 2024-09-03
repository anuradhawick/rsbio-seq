from typing import Iterator, Union

class Sequence:
    """
    Sequence entry.

    Attributes:
        id (str): record id of the sequence.
        desc (str): description of the record if available.
        seq (str): content of the sequence.
        qual (str): quality string of the sequence if available.
    """

    id: str
    desc: str
    seq: str
    qual: str

    def __init__(self, id: str, desc: str, seq: str, qual: str = "") -> None:
        """
        Initialise the sequence instance.

        Args:
            id (str): record id of the sequence.
            desc (str): description of the record if available.
            seq (str): content of the sequence.
            qual (str): optional quality string of the sequence if available.
                        this is discarded if the writing format is fasta.
        """
        ...

class SeqReader:
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

class SeqWriter:
    """
    Sequence writer.
    """

    def __init__(self, path: str) -> None:
        """
        Initialise the reader with path of the file.

        Args:
            path (str): The path to file fasta, fastq, fa, fq and compressed formats with gz are supported.
                        The file path is used to infer the intentded file format.
        """
        ...

    def write(self, seq: Sequence, wrap: Union[int, None] = None) -> None:
        """
        Writes a sequence to file.

        Args:
            seq (Sequence): the Sequence instance you want to write
            wrap (Union[int, None]): wrap length of the written sequence (default: None)
        """
        ...

    def close(self) -> None:
        """
        End writing the sequence to file.
        """
        ...

__all__ = ["Sequence", "SeqReader", "SeqWriter"]
