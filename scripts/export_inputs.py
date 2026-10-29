#!/usr/bin/env python3
"""Materialize documented inputs in a NEW dataset directory; run large copies in an allocation."""
import argparse
import hashlib
import json
from pathlib import Path


def selected_inputs(manifest, source_root, groups):
    source_root = source_root.resolve()
    selected = []
    seen = set()
    for entry in manifest['entries']:
        if entry['group'] not in groups:
            continue
        relative = Path(entry['path'])
        if relative.is_absolute() or '..' in relative.parts or relative.as_posix() in seen:
            raise ValueError('unsafe or duplicate destination: ' + str(relative))
        seen.add(relative.as_posix())
        source_relative = Path(entry['source'])
        if source_relative.is_absolute() or '..' in source_relative.parts:
            raise ValueError('unsafe source: ' + str(source_relative))
        source = (source_root / source_relative).resolve()
        source.relative_to(source_root)
        if not source.is_file() or source.stat().st_size != entry['size_bytes']:
            raise ValueError('missing or size-mismatched input: ' + entry['source'])
        selected.append((entry, source))
    if not selected:
        raise ValueError('no inputs selected')
    return selected


def export(manifest, source_root, out, groups):
    if out.exists():
        raise FileExistsError('choose a new dataset directory: ' + str(out))
    inputs = selected_inputs(manifest, source_root, groups)
    out.mkdir()  # Parent must already exist; preserve partial outputs after failures.
    copied = []
    for entry, source in inputs:
        target = out / entry['path']
        target.parent.mkdir(parents=True, exist_ok=True)
        digest = hashlib.sha256()
        size = 0
        with source.open('rb') as src, target.open('xb') as dst:
            for block in iter(lambda: src.read(1024 * 1024), b''):
                dst.write(block)
                digest.update(block)
                size += len(block)
        actual = digest.hexdigest()
        if size != entry['size_bytes'] or (entry.get('sha256') and actual != entry['sha256']):
            raise ValueError('input changed or checksum mismatch: ' + entry['source'])
        copied.append(dict(path=entry['path'], source=entry['source'], group=entry['group'],
                           size_bytes=size, sha256=actual))
        print('Copied ' + entry['path'], flush=True)
    with (out / 'dataset-manifest.json').open('x') as fh:
        json.dump(dict(schema_version=1, entries=copied), fh, indent=2)
        fh.write('\n')
    return copied


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'docs/input-manifest.json')
    parser.add_argument('--source-root', required=True, type=Path,
                        help='original Fc_FcR3 directory containing prep and original WT files')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--groups', nargs='+', choices=['analysis', 'preparation', 'simulation'],
                        default=['analysis', 'preparation'])
    args = parser.parse_args()
    export(json.loads(args.manifest.read_text()), args.source_root, args.out, set(args.groups))


if __name__ == '__main__':
    main()
