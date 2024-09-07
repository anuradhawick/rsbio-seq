from typing import Iterator, Union, List

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

def phred_to_ascii(scores: List[int]) -> str:
    """
    Convert a list of Phred scores to a string of ASCII quality values.

    Each Phred score is converted by adding 33 to it (Phred+33 encoding),
    and then transforming it into the corresponding ASCII character.

    Args:
        scores (List[int]): A list of integers where each integer represents a Phred score.

    Returns:
        str: A string where each character represents an ASCII quality value corresponding to the input Phred scores.
    """
    ...

def ascii_to_phred(qual: str) -> List[int]:
    """
    Convert a string of ASCII quality values back to a list of Phred scores.

    Each ASCII character is converted back to a Phred score by subtracting 33
    from its ASCII value. This function assumes standard Phred+33 encoding is used.

    Args:
        qual (str): A string of ASCII characters representing quality scores.

    Returns:
        List[int]: A list of integers where each integer is a Phred score derived from the input ASCII characters.
    """
    ...

__all__ = ["Sequence", "SeqReader", "SeqWriter", "phred_to_ascii", "ascii_to_phred"]
