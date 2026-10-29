#!/usr/bin/env python3
"""Write a new protocol with paths relative to its location; no trajectory reads."""
import argparse
import json
import os
from pathlib import Path


def configure_protocol(template, data_root, reference, out):
    """Write a portable protocol into a new file."""
    cfg = json.loads(template.read_text())
    cfg['prep_root'] = os.path.relpath(data_root.resolve(), out.resolve().parent)
    cfg['reference_pdb'] = os.path.relpath(reference.resolve(), out.resolve().parent)
    with out.open('x') as fh:
        json.dump(cfg, fh, indent=2)
        fh.write('\n')
    print(f'Protocol written: {out}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--template', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'config/analysis_protocol.json')
    args = parser.parse_args()
    configure_protocol(args.template, args.data_root, args.reference, args.out)


if __name__ == '__main__':
    main()
