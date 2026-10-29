"""Shared checksum and output-table helpers; existing outputs are protected."""
import csv
import hashlib
from pathlib import Path


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_csv(path, fields, rows):
    with Path(path).open('x', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
