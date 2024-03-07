import argparse
from .pre_processing import *
from .fingerprint import *
from .pca_transform import *
from .utils import *
from .similarity import *
from .version import __version__

def main():
    
    def parse_dict_or_none(argument):
        if argument == 'None':
            return None
        elif argument == 'PROTON_FEATURES':
            return PROTON_FEATURES
        elif argument == 'DEFAULT_FEATURES':  
        # Ensure DEFAULT_FEATURES is a defined dictionary of features
            return DEFAULT_FEATURES
        else:
            print(f"Warning: Unrecognized feature set '{argument}'. Using DEFAULT_FEATURES.")
            return DEFAULT_FEATURES
    
    parser = argparse.ArgumentParser(description="HSR (Hyper Shape Recognition) CLI for molecule comparison.",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s', '--similarity', nargs=2, metavar=('mol1', 'mol2'), 
                        help='Calculate similarity between two molecule files. Requires exactly two molecule files.')
    parser.add_argument('-d', '--distance', nargs=2, metavar=('mol1', 'mol2'), 
                        help='Calculate distance between two molecule files. Requires exactly two molecule files.')
    parser.add_argument('-f', '--fingerprint', nargs='+', metavar='molecule', 
                        help='Generate fingerprint for one or more molecule files.')
    
    parser.add_argument('-removeH', action='store_true', help='Remove hydrogen atoms')
    parser.add_argument('-chirality', action='store_true', help='Consider chirality')
    parser.add_argument('-features', default='DEFAULT_FEATURES', type=parse_dict_or_none,
                        help="Defines the use of additional features for the comparison. "
                             "Available options: DEFAULT_FEATURES for 6D representation (default), 'None' for only spatial coordinates, and"
                             "PROTON_FEATURES for 4D representation including proton feature besides spatial coordinates. "
                             "For different features, it is necessary to define new ones in the package (see documentation).")
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}',
                        help='Show the version of the HSR package and exit.')

    args = parser.parse_args()

    # Handling of each operation
    if args.similarity:
        mol1, mol2 = args.similarity
        # Implement the similarity calculation using the provided molecule files
        print(f"Calculating similarity for {mol1} and {mol2}")
        if args.removeH:
            mol1 = read_mol_from_file(mol1, removeHs=True, sanitize=True)
            mol2 = read_mol_from_file(mol2, removeHs=True, sanitize=True)
        else:
            mol1 = read_mol_from_file(mol1)
            mol2 = read_mol_from_file(mol2)
        if mol1 is None or mol2 is None:
            print("Error reading molecule files")
            return
        # Compute the similarity
        similarity = compute_similarity(mol1, mol2, features=args.features, chirality=args.chirality)
        print(f"Similarity: {similarity:.4f}")
        
    elif args.distance:
        mol1, mol2 = args.distance
        # Implement the distance calculation
        print(f"Calculating distance for {mol1} and {mol2}")
        if args.removeH:
            mol1 = read_mol_from_file(mol1, removeHs=True, sanitize=True)
            mol2 = read_mol_from_file(mol2, removeHs=True, sanitize=True)
        else:
            mol1 = read_mol_from_file(mol1)
            mol2 = read_mol_from_file(mol2)
        if mol1 is None or mol2 is None:
            print("Error reading molecule files")
            return
        # Compute the distance
        distance, _ = compute_distance(mol1, mol2, features=args.features, removeHs=args.removeH, chirality=args.chirality)
        print(f"Manhattan Distance: {distance:.4f}")

    elif args.fingerprint:
        # Implement the fingerprint generation for one or more molecule files
        for molecule in args.fingerprint:
            print(f"Generating fingerprint for {molecule}")
            if args.removeH:
                mol = read_mol_from_file(molecule, removeHs=args.removeH, sanitize=True)
            else:
                mol = read_mol_from_file(molecule)
            if mol is None:
                print(f"Error reading molecule file {molecule}")
                continue
            # Generate the fingerprint
            fingerprint = generate_fingerprint_from_molecule(mol, features=args.features, removeHs=args.removeH, chirality=args.chirality)
            print(f"Fingerprint for {molecule}: {fingerprint}")
        

if __name__ == '__main__':
    main()