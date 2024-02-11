import os
import sys
from hsr.pre_processing import *
sys.path.append(os.path.abspath('../'))
from csr import *

cwd = os.getcwd()

molecules = load_molecules_from_sdf(f'{cwd}/lambda_delta_isomers.sdf', removeHs=False, sanitize=False)

n_molecules = len(molecules)
print(f'\n(in-house) CSR Similarity of inorganic enantiomers:')

similarity = compute_similarity(molecules[0], molecules[1])
print(f"CSR: {similarity}")


