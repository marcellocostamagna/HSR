import os
from hsr.pre_processing import *
from experiments.usr_optiso import *

cwd = os.getcwd()

molecules = load_molecules_from_sdf(f'{cwd}/experiments/Chirality/lambda_delta_isomers.sdf', removeHs=False, sanitize=False)

print(f'\n(in-house) USR:OptIso Similarity of inorganic enantiomers:')
similarity = compute_similarity(molecules[0], molecules[1])
       
print(f"USR:OptIso : {similarity}")


