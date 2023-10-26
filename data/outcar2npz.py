import argparse
import numpy as np
from dpdata import LabeledSystem


ATOM_DICT = {
    'H': 1,
    'B': 5,
    'C': 6,
    'N': 7,
    'O': 8,
    'F': 9,
    'Ti': 22,
}


def name2numb(names):
    """ convert a list of atom symbols to a list of atom numbers
    """
    numbs = []
    for name in names:
        numbs.append(ATOM_DICT[name])
    return numbs


def outcar2npz(fpath, fname):
    data = LabeledSystem(fpath, fmt="vasp/outcar").data
    cells = data['cells']
    coords = data['coords']
    forces = data['forces']
    energies = data['energies']

    atom_list = np.array(name2numb(data['atom_names']))
    # atom numbers for one frame, because the structure of each frame is the same.
    atom_numbs = atom_list[data['atom_types']]

    np.savez(f'{fname}.npz', cells=cells, coords=coords, forces=forces,
             energies=energies, atom_numbs=atom_numbs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fpath', type=str, default='OUTCAR', help='file path of OUTCAR')
    parser.add_argument('--fname', type=str, default='data', help='file name of output npz')
    args = parser.parse_args()
    outcar2npz(args.fpath, args.fname)
