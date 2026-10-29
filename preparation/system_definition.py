"""Shared glycan connectivity, chain identities, and disulfide definitions."""

NAG_RENAME   = {'C7': 'C2N', 'C8': 'CME', 'O7': 'O2N'}
GLYCAN_CODES = {1: 'UYB', 2: '4YB', 3: 'VMB', 4: '2MA', 5: '0YB', 6: '2MA', 7: '0YB', 8: '0fA'}
GLYCOSIDIC   = {2: (1, 'O4'), 3: (2, 'O4'), 4: (3, 'O3'), 5: (4, 'O2'),
                6: (3, 'O6'), 7: (6, 'O2'), 8: (1, 'O6')}      # child: (parent, parent O)
ASN_GLYCAN   = {('A', 297): 'D', ('B', 297): 'E'}
STUBS        = {('C', 162): 1001, ('C', 45): 1002}
DISULFIDES   = [('A', 261, 321), ('A', 367, 425), ('B', 261, 321),
                ('B', 367, 425), ('C', 29, 71), ('C', 110, 154)]
PROTEIN_CHAINS = ['A', 'B', 'C']
