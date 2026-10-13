"""Compare SeqIO.parse with old/new release-built RSBio-Seq modules.

Usage: python benchmarks/readers.py --bio-path /path/to/old/module \
    --latest-path /path/to/new/module
Each module path must contain an importable rsbio_seq extension/package.
"""
import argparse
import gc
import json
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def worker(args):
    if args.worker == "biopython":
        import Bio
        from Bio import SeqIO
        version = Bio.__version__

        def scan():
            with open(args.input) as handle:
                return sum(1 for _ in SeqIO.parse(handle, args.format))
    else:
        sys.path.insert(0, args.module_path)
        import rsbio_seq
        version = "1.0.0"

        def scan():
            return sum(1 for _ in rsbio_seq.SeqReader(args.input))

    assert scan() == args.records  # Warm file cache and parser before timing.
    elapsed = []
    for _ in range(args.runs):
        gc.collect()
        start = time.perf_counter()
        count = scan()
        elapsed.append(time.perf_counter() - start)
        assert count == args.records
    median = statistics.median(elapsed)
    print(json.dumps({"version": version, "seconds": median,
                      "mb_per_second": Path(args.input).stat().st_size / median / 1e6,
                      "records_per_second": args.records / median,
                      "trials_seconds": elapsed}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bio-path")
    parser.add_argument("--latest-path")
    parser.add_argument("--records", type=int, default=200_000)
    parser.add_argument("--runs", type=int, default=7)
    parser.add_argument("--worker", choices=["biopython", "rsbio"])
    parser.add_argument("--module-path")
    parser.add_argument("--input")
    parser.add_argument("--format", choices=["fasta", "fastq"])
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    if not args.bio_path or not args.latest_path:
        parser.error("--bio-path and --latest-path are required")
    results = {}
    with tempfile.TemporaryDirectory(prefix="rsbio-benchmark-") as directory:
        for fmt, suffix, record in [
            ("fasta", "fa", ">read description\n" + "ACGT" * 37 + "AC\n"),
            ("fastq", "fq", "@read description\n" + "ACGT" * 37 + "AC\n+\n" + "I" * 150 + "\n"),
        ]:
            path = Path(directory) / ("reads." + suffix)
            path.write_text(record * args.records)
            results[fmt] = {"bytes": path.stat().st_size, "implementations": {}}
            for name, mode, module in [
                ("Biopython", "biopython", None),
                ("Bio crate", "rsbio", args.bio_path),
                ("Latest implementation", "rsbio", args.latest_path),
            ]:
                command = [sys.executable, str(Path(__file__).resolve()),
                           "--worker", mode, "--input", str(path), "--format", fmt,
                           "--records", str(args.records), "--runs", str(args.runs)]
                if module:
                    command += ["--module-path", str(Path(module).resolve())]
                output = subprocess.check_output(command, text=True)
                results[fmt]["implementations"][name] = json.loads(output)
    print(json.dumps({"records": args.records, "bases_per_record": 150,
                      "runs": args.runs, "results": results}, indent=2))


if __name__ == "__main__":
    main()
