import os
from ase.io import read, write


def find_files(fname):
    """ Recursively traverse folders and subfolders to find files you want
    """
    # Set the directory path to the current directory
    directory_path = os.getcwd()

    fpath_list = []
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file == fname:
                fpath_list.append(os.path.join(root, file))
    return fpath_list


if __name__ == "__main__":
    outcar_list = find_files('OUTCAR')
    combined_trajectory = []
    for outcar in outcar_list:
        vasp_trajectory = read(outcar, index=':')
        combined_trajectory.extend(vasp_trajectory)

    write('all_data.extxyz', combined_trajectory, format='extxyz')
