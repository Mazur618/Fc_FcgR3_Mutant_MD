"""Extract the protein reference from the source CIF."""


def extract_protein(source='7URU.cif', out='protein_only.pdb'):
    """Remove hydrogens, alternate conformations, ligands, and waters."""
    import gemmi
    structure = gemmi.read_structure(source)
    structure.remove_hydrogens()
    structure.remove_alternative_conformations()
    structure.remove_ligands_and_waters()
    structure.write_pdb(out)


def main():
    extract_protein()


if __name__ == '__main__':
    main()
