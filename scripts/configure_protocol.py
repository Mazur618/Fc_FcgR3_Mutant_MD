#!/usr/bin/env python3
"""Write a new protocol with paths relative to its location; no trajectory reads."""
import argparse
import json
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--template', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'analysis_v2/protocol.json')
    args = parser.parse_args()
    cfg = json.loads(args.template.read_text())
    cfg['prep_root'] = os.path.relpath(args.data_root.resolve(), args.out.resolve().parent)
    cfg['reference_pdb'] = os.path.relpath(args.reference.resolve(), args.out.resolve().parent)
    with args.out.open('x') as fh:
        json.dump(cfg, fh, indent=2)
        fh.write('\n')
    print(f'Protocol written: {args.out}')


if __name__ == '__main__':
    main()
