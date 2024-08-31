from tqdm import tqdm
from Bio import SeqIO as BioSeqIO
import sys
from rsbio_seq import SeqIO as RsSeqIO


if __name__ == "__main__":
    assert len(sys.argv) == 3, "Must provide a valid sequence file and format"
    for seq in tqdm(
        RsSeqIO(sys.argv[1]),
        desc="Rs-iteration",
    ):
        pass

    for seq in tqdm(
        BioSeqIO.parse(sys.argv[1], sys.argv[2]),
        desc="Bio-iteration",
    ):
        pass

# Rs-iteration:  28963522it [01:16, 379747.32it/s]
# Bio-iteration: 28963522it [02:43, 177329.69it/s]
