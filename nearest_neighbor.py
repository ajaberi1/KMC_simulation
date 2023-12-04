# import numpy as np
# from ase import atoms
# from ase import atom
# import math
# from ase.io import read, write
# import time as t
# # import concurrent.futures
# from ase.neighborlist import NeighborList, natural_cutoffs
# import pandas as pd
# import json


def find_neighbor_from_DataFrame(all_neighbors_DataFrame, selected_atom):
    selected_atom_tag = selected_atom.tag
    neigbors_plus_nones = all_neighbors_DataFrame[selected_atom_tag].values
    nearest_atom_object = [
        item for item in neigbors_plus_nones if item != None]
    return nearest_atom_object


def find_nearest_neighbor_with_PBC(atoms, selected_atom):
    # NOTE: To get the NeighborList object you should have these four following lines outside of your function:
    # from ase.neighborlist import NeighborList, natural_cutoffs
    # my_cutoffs = natural_cutoffs(atoms)
    # nl = NeighborList(cutoffs=my_cutoffs, bothways=True, self_interaction=False)
    # nl.update(atoms)
    nearest_atom_object = []
    # indices, offsets = nl.get_neighbors(selected_atom.index)
    # for index in indices:
    #     if atoms[index].symbol == 'Li':
    #         nearest_atom_object.append(atoms[index])

    return nearest_atom_object


def find_nearest_neighbor(atoms, selected_atom):
    # print('Selected atom:  ', selected_atom)
    # atoms = atom_selectedatom[0]
    # selected_atom = atom_selectedatom[1]
    distances_for_indexes = {}
    for atom in atoms:
        # first if to check if they are in the same plane
        if abs(atom.position[2] - selected_atom.position[2]) < 0.001:
            if abs(atom.position[1] - selected_atom.position[1]) < 5 and abs(atom.position[0] - selected_atom.position[0]) < 5:
                if atom.symbol == 'Li' or atom.symbol == 'X':
                    one_distances = atoms.get_distance(
                        selected_atom.index, atom.index)
                    distances_for_indexes[atom.index] = one_distances
                    # append((one_distances, atom.index))
    nearest_atom_object = []
    for index in distances_for_indexes:
        if distances_for_indexes[index] < 3 and distances_for_indexes[index] != 0:
            nearest_atom_object.append(atoms[index])

    return nearest_atom_object


def find_nearest_neighbor_old(atoms, selected_atom):

    # print('Selected atom:  ', selected_atom)
    distances_for_indexes = {}
    for atom in atoms:
        if atom.symbol == 'Li' or atom.symbol == 'X':
            one_distances = atoms.get_distance(
                selected_atom.index, atom.index)
            distances_for_indexes[atom.index] = one_distances
            # append((one_distances, atom.index))

    nearest_atom_object = []
    for index in distances_for_indexes:
        if distances_for_indexes[index] < 3 and distances_for_indexes[index] != 0:
            nearest_atom_object.append(atoms[index])

    return nearest_atom_object


# atoms = read(filename="small_new.xsf")
# selected_atom = atoms[917]
# print("selected_atom: ", selected_atom)
# print("")
# my_cutoffs = natural_cutoffs(atoms)
# nl = NeighborList(cutoffs=my_cutoffs, bothways=True, self_interaction=False)
# nl.update(atoms)
# my_neighbor = find_nearest_neighbor_new(atoms, selected_atom)
# for i in my_neighbor:
#     print(i)
# print("")
# nearest_atom_object = find_nearest_neighbor(atoms, selected_atom)
# for i in nearest_atom_object:
#     print(i)
# time1 = t.time()
# for atom in atoms:
#     if atom.symbol == 'Li':
#         my_neighbor = find_nearest_neighbor_old(atoms, atom)
# print("time for the first one: ", t.time() - time1)
# print("")
# # print(indices)
# time2 = t.time()
# for atom in atoms:
#     if atom.symbol == 'Li':
#         nearest_atom_object = find_nearest_neighbor(atoms, selected_atom)
# print("time for the second one: ", t.time() - time2)


# df = pd.DataFrame()
# nearest_atom_object = find_nearest_neighbor(atoms, selected_atom)
# dic = {selected_atom: nearest_atom_object}
# print(dic)
# print(df)
# print("")
# print(df["1"])

# for item in new_nearest_atom:
#     item.symbol = 'Cs'
# atoms[selected_atom.index].symbol = 'Mg'
# write(filename='neighbors1.xsf', images=atoms)

# print(nearest_atom_object)
# all_li = [item for item in atoms if item.symbol == 'Li']
# start_time1 = t.time()
# for Li in all_li:
#     start_time0 = t.time()
#     nearest_atom_object = find_nearest_neighbor(atoms, Li)
#     # print("time for one: ", t.time() - start_time0)
#     print(len(nearest_atom_object))
# print("time to find all the neighbors: ", t.time() - start_time1)
# start_time2 = t.time()
# with concurrent.futures.ThreadPoolExecutor() as executor:
#     results = [executor.submit(find_nearest_neighbor, [atoms, Li])
#                for Li in all_li]
#     for f in concurrent.futures.as_completed(results):
#         # print(f.result())
#         pass
# # for item in results:
# #     print(item.result)
# print("time to find all the neighbors with thread: ", t.time() - start_time2)

# for item in nearest_atom_object:
#     print(item.index)
#     atoms[item.index].symbol = 'Cs'
# atoms[selected_atom.index].symbol = 'Mg'
# write(filename='neighbors.xsf', images=atoms)

# start_time = t.time()
# for atom in atoms:
#     if atom.symbol == "Li":
#         near = find_nearest_neighbor(atoms, atom)
# print("time: ", t.time() - start_time)
# print("")

# start_time = t.time()
# for atom in atoms:
#     if atom.symbol == "Li":
#         near = find_nearest_neighbor_old(atoms, atom)
# print("time: ", t.time() - start_time)
# print("Done")
