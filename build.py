#!/usr/bin/env python3
# Run in the PyRosetta env (Rosetta is only imported if you pass --mut).
import argparse, math, sys
import gemmi

ap = argparse.ArgumentParser()
ap.add_argument('--name', required=True, help='system name, used for all output files')
ap.add_argument('--protein', default='WT.pdb',
                help='starting protein PDB (chains A,B,C; already carries the allotype mutation)')
ap.add_argument('--cif', default='7URU.cif',
                help='source of the glycans and CD16a NAG stubs only')
ap.add_argument('--mut', action='append', default=[],
                help='CHAIN:NUM:AA3, e.g. C:158:PHE or Fc:332:GLU (Fc = both A and B). Repeatable.')
ap.add_argument('--his', action='append', default=[],
                help='CHAIN:NUM:HID|HIE|HIP, e.g. A:310:HIP. Repeatable.')
ap.add_argument('--radius', type=float, default=8.0,
                help='repack only residues within this distance (A) of mutated sites; 0 = repack all')
ap.add_argument('--buffer', type=float, default=8.0, help='solvent buffer (A)')
ap.add_argument('--vac', action='store_true', help='also save a vacuum prmtop/rst7 (off by default)')
ap.add_argument('--ions2', action='store_true', help='use faster addIons2 instead of addIons')
args = ap.parse_args()

# ---------------- system definition (edit only if the glycoform changes) ----------------
NAG_RENAME   = {'C7': 'C2N', 'C8': 'CME', 'O7': 'O2N'}
GLYCAN_CODES = {1: 'UYB', 2: '4YB', 3: 'VMB', 4: '2MA', 5: '0YB', 6: '2MA', 7: '0YB', 8: '0fA'}
GLYCOSIDIC   = {2: (1, 'O4'), 3: (2, 'O4'), 4: (3, 'O3'), 5: (4, 'O2'),
                6: (3, 'O6'), 7: (6, 'O2'), 8: (1, 'O6')}      # child: (parent, parent O)
ASN_GLYCAN   = {('A', 297): 'D', ('B', 297): 'E'}
STUBS        = {('C', 162): 1001, ('C', 45): 1002}
DISULFIDES   = [('A', 261, 321), ('A', 367, 425), ('B', 261, 321),
                ('B', 367, 425), ('C', 29, 71), ('C', 110, 154)]
PROTEIN_CHAINS = ['A', 'B', 'C']
# -----------------------------------------------------------------------------------------

cys = {(c, n) for c, a, b in DISULFIDES for n in (a, b)}
asn = set(ASN_GLYCAN) | set(STUBS)

def parse_sites(specs):
    out = []
    for s in specs:
        ch, num, val = s.split(':')
        for c in (['A', 'B'] if ch.lower() == 'fc' else [ch]):
            out.append((c, int(num), val.upper()))
    return out

muts = parse_sites(args.mut)
his = {(c, n): v for c, n, v in parse_sites(args.his)}
for v in his.values():
    assert v in ('HID', 'HIE', 'HIP'), f'bad histidine state {v}'

for c, n, aa in muts:
    if (c, n) in cys or (c, n) in asn:
        sys.exit(f'ERROR: {c}{n} is a disulfide Cys or glycosylated Asn. '
                 f'Edit DISULFIDES/ASN_GLYCAN/STUBS if you really mean it.')

# ---------------- step 1: mutate and repack (PyRosetta, only if --mut) ----------------
protein_file = args.protein
if muts:
    import pyrosetta
    from pyrosetta import rosetta
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import (
        RestrictToRepacking, OperateOnResidueSubset, PreventRepackingRLT)

    pyrosetta.init('-mute all -ignore_unrecognized_res')
    pose = pyrosetta.pose_from_pdb(args.protein)
    idxs = []
    for c, n, aa in muts:
        i = pose.pdb_info().pdb2pose(c, n)
        if i == 0:
            sys.exit(f'ERROR: {c}{n} not found in {args.protein}')
        print(f'mutating {c}{n}: {pose.residue(i).name3()} -> {aa} (pose index {i})')
        rosetta.protocols.simple_moves.MutateResidue(i, aa).apply(pose)
        idxs.append(i)

    sfxn = rosetta.core.scoring.get_score_function()
    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())
    if args.radius > 0:
        sel = NeighborhoodResidueSelector(
            ResidueIndexSelector(','.join(map(str, idxs))), args.radius, True)
        tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), sel, True))
    task = tf.create_task_and_apply_taskoperations(pose)
    rosetta.protocols.minimization_packing.PackRotamersMover(sfxn, task).apply(pose)
    protein_file = f'{args.name}_protein.pdb'
    pose.dump_pdb(protein_file)
    with open(f'{args.name}_mutations.txt', 'w') as f:
        f.write('\n'.join(f'{c}{n}{aa}' for c, n, aa in muts) + '\n')

print(f'protein source: {protein_file}   glycan source: {args.cif}')

# ---------------- step 2: assemble protein + glycans ----------------
cif = gemmi.read_structure(args.cif)
cif.remove_hydrogens(); cif.remove_alternative_conformations()
cmodel = cif[0]
pst = gemmi.read_structure(protein_file)
pst.remove_hydrogens()                      # tleap rebuilds all hydrogens
pmodel = pst[0]

residues, index = [], {}
def add(key, resname, res, rename=None, ter=False):
    atoms = {(rename or {}).get(a.name, a.name): (a.pos.x, a.pos.y, a.pos.z, a.element.name)
             for a in res}
    residues.append({'name': resname, 'atoms': atoms, 'ter': ter})
    index[key] = len(residues)

prot_by_chain = {}
for cname in PROTEIN_CHAINS:
    prot = [r for r in pmodel[cname] if r.het_flag == 'A' and r.name != 'HOH']
    prot_by_chain[cname] = prot
    print(f'chain {cname}: {len(prot)} residues')
    for i, r in enumerate(prot):
        key = (cname, r.seqid.num)
        if key in asn:
            assert r.name == 'ASN', f'{key} is {r.name}, not ASN'
            resname = 'NLN'
        elif key in cys:
            assert r.name == 'CYS', f'{key} is {r.name}, not CYS'
            resname = 'CYX'
        elif key in his:
            assert r.name in ('HIS', 'HID', 'HIE', 'HIP'), f'{key} is {r.name}, not histidine'
            resname = his[key]
        else:
            resname = r.name
        add(key, resname, r, ter=(i == len(prot) - 1))

for gch in sorted(set(ASN_GLYCAN.values())):
    for r in cmodel[gch]:
        add((gch, r.seqid.num), GLYCAN_CODES[r.seqid.num], r,
            rename=NAG_RENAME if r.name == 'NAG' else None, ter=True)
for (cname, n), sn in STUBS.items():
    r = [x for x in cmodel[cname] if x.seqid.num == sn and x.name == 'NAG'][0]
    add(('stub', n), '0YB', r, rename=NAG_RENAME, ter=True)

for (pch, n) in ASN_GLYCAN:
    r3 = [r for r in prot_by_chain[pch] if r.seqid.num == n + 2][0]
    if r3.name not in ('THR', 'SER'):
        print(f'WARNING: {pch}{n + 2} is {r3.name}; the N{n} sequon is broken but the glycan is still attached')

# ---------------- step 3: bonds ----------------
bonds = []
for c, a, b in DISULFIDES:
    bonds.append((index[(c, a)], 'SG', index[(c, b)], 'SG'))
for (pch, n), gch in ASN_GLYCAN.items():
    bonds.append((index[(pch, n)], 'ND2', index[(gch, 1)], 'C1'))
for gch in sorted(set(ASN_GLYCAN.values())):
    for child, (parent, oatom) in GLYCOSIDIC.items():
        bonds.append((index[(gch, parent)], oatom, index[(gch, child)], 'C1'))
for (c, n) in STUBS:
    bonds.append((index[(c, n)], 'ND2', index[('stub', n)], 'C1'))

print('== bond distances (A) ==')
bad = 0
for i, a, j, b in bonds:
    d = math.dist(residues[i-1]['atoms'][a][:3], residues[j-1]['atoms'][b][:3])
    limit = 2.8 if (a == 'SG' and b == 'SG') else 2.0
    flag = d > limit
    bad += flag
    print(f"{residues[i-1]['name']}{i}.{a} - {residues[j-1]['name']}{j}.{b}  {d:.2f}"
          + ('   <-- CHECK' if flag else ''))
if bad:
    print(f'\nWARNING: {bad} bond(s) beyond tolerance. Asn-NAG or glycosidic bonds of many A '
          f'mean the protein is in a different frame from {args.cif}; long S-S bonds are '
          f'usually just strained and relax in minimization.')

# ---------------- step 4: write PDB and tleap input ----------------
n = args.name
with open(f'{n}_glycam.pdb', 'w') as f:
    serial = 0
    for k, r in enumerate(residues, 1):
        for an, (x, y, z, el) in r['atoms'].items():
            serial += 1
            nm = (' ' + an.ljust(3)) if len(an) < 4 else an
            f.write(f"ATOM  {serial:5d} {nm} {r['name']:<3s} A{k:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n")
        if r['ter']:
            f.write('TER\n')
    f.write('END\n')

ion = 'addIons2' if args.ions2 else 'addIons'
with open(f'build_{n}.in', 'w') as f:
    # GLYCAM first: it loads TIP3P ion params, which the OPC3 ion frcmod below must override
    f.write('source leaprc.protein.ff19SB\n')
    f.write('source leaprc.GLYCAM_06j-1\n')
    f.write('source leaprc.water.opc3\n')
    f.write('loadAmberParams frcmod.ionslm_126_opc3\n\n')
    f.write(f'mol = loadPdb {n}_glycam.pdb\n\n')
    for i, a, j, b in bonds:
        f.write(f'bond mol.{i}.{a} mol.{j}.{b}\n')
    f.write('\ncheck mol\n')
    if args.vac:
        f.write(f'saveAmberParm mol {n}_vac.prmtop {n}_vac.rst7\n')
    f.write(f'\nsolvateOct mol OPC3BOX {args.buffer}\n')
    f.write(f'{ion} mol Cl- 0\n')
    f.write(f'{ion} mol Na+ 0\n')
    f.write(f'saveAmberParm mol {n}.prmtop {n}.rst7\n')
    f.write(f'savePdb mol {n}_solv.pdb\n')
    f.write('quit\n')
print(f'\n{len(residues)} residues, {len(bonds)} bonds -> build_{n}.in')
