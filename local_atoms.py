import numpy as np
from ase import atoms
from ase import atom
import math
import csv
from ase.io import read, write
import time as t


def find_local_atoms(atoms, selected_atom, local_atom_type, threshold_distance=8):
    distances_for_indexes = {}
    for atom in atoms:
        if atom.symbol == local_atom_type:
            one_distances = atoms.get_distance(
                selected_atom.index, atom.index)
            distances_for_indexes[atom.index] = one_distances
            # append((one_distances, atom.index))

    # NOTE: If you also want to get the nearest_atom_object with this function
    # nearest_atom_object = []
    # for index in distances_for_indexes:
    #     if distances_for_indexes[index] < 3 and distances_for_indexes[index] != 0:
    #         nearest_atom_object.append(atoms[index])

    local_atoms = []
    for index in distances_for_indexes:
        if distances_for_indexes[index] < threshold_distance and distances_for_indexes[index] != 0:
            local_atoms.append(atoms[index])

    # NOTE: If you also want to get the local li_concentration with this function
    # li_count = 0
    # x_count = 0
    # for atom in local_atoms:
    #     if atom.symbol == 'Li':
    #         li_count += 1
    #     if atom.symbol == 'X':
    #         x_count += 1
    # li_concentration = li_count/(li_count+x_count)

    return local_atoms


# atoms = read(filename="small_new.xsf")
# selected_atom = atoms[1410]
# print(selected_atom)
# print("")
# loc = find_local_atoms(atoms, selected_atom,
#                        local_atom_type='Li', threshold_distance=8)
# for item in loc:
#     print(item)
# atoms[selected_atom.index].symbol = 'Cs'
# for item in loc:
#     atoms[item.index].symbol = 'Mg'
# write(filename=f"before_structure1.xsf", images=atoms)

