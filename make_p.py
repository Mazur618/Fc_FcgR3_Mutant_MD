# make_protein.py
import gemmi
st = gemmi.read_structure('7URU.cif')
st.remove_hydrogens()
st.remove_alternative_conformations()
st.remove_ligands_and_waters()
st.write_pdb('protein_only.pdb')
