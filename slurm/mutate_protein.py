#!/usr/bin/env python3
"""Mutate and repack a protein using an explicitly invoked PyRosetta workflow."""
import argparse


def argument_parser():
    """Describe the existing mutation command-line options."""
    ap = argparse.ArgumentParser()
    ap.add_argument("inpdb")
    ap.add_argument("--resi", type=int, help="pose residue index (1-based)")
    ap.add_argument("--chain", help="chain ID (alternative to --resi)")
    ap.add_argument("--num", type=int, help="PDB residue number (with --chain)")
    ap.add_argument("--aa", required=True, help="new AA three-letter code")
    ap.add_argument("--radius", type=float, default=0.0,
                    help="repack only residues within this distance (A) of the "
                         "mutation; 0 = repack everything")
    ap.add_argument("-o", "--out", default="mut_repacked.pdb")
    return ap


def parse_args(argv=None):
    return argument_parser().parse_args(argv)


def resolve_residue(pose, args):
    """Resolve the target pose index or report an invalid chain/residue request."""
    parser = argument_parser()
    if args.resi is None:
        if args.chain is None or args.num is None:
            parser.error("give --resi, or both --chain and --num")
        args.resi = pose.pdb_info().pdb2pose(args.chain, args.num)
        if args.resi == 0:
            parser.error(f"{args.chain}{args.num} not found in {args.inpdb}")
    return args.resi


def mutate_protein(args):
    """Load, mutate, repack, and write the protein with the original settings."""
    import pyrosetta
    from pyrosetta import rosetta
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import (
        RestrictToRepacking, OperateOnResidueSubset, PreventRepackingRLT)

    pyrosetta.init("-mute all")
    pose = pyrosetta.pose_from_pdb(args.inpdb)
    args.resi = resolve_residue(pose, args)
    print(f"Mutating pose residue {args.resi} "
          f"({pose.residue(args.resi).name3()} -> {args.aa})")

    rosetta.protocols.simple_moves.MutateResidue(args.resi, args.aa).apply(pose)

    scorefxn = rosetta.core.scoring.get_score_function()
    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())
    if args.radius > 0:
        sel = NeighborhoodResidueSelector(
            ResidueIndexSelector(str(args.resi)), args.radius, True)
        # flip_subset=True: freeze everything NOT in the neighborhood
        tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), sel, True))
    task = tf.create_task_and_apply_taskoperations(pose)

    rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, task).apply(pose)
    pose.dump_pdb(args.out)


def main():
    mutate_protein(parse_args())


if __name__ == '__main__':
    main()
