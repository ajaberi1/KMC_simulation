import numpy as np
import math
import random
# from ase import atoms
# from ase import atom
# import csv
# from ase.io import read, write
# import time as t

# NOTE: each element of list_all_rates has the following items: ['rate', 'atom',  'atom.tag', 'target_end_point','target_end_point.tag', 'hopping_mechanism', 'barrier', 'kra']


def find_jump_rate_atom_end_timestp_type_barrier(atoms, temp, list_all_rates):
    # list_all_rates = np.zeros((0, 6))
    # with open(f"list_all_rates.csv") as file:
    #     reader = csv.reader(file)
    #     for row in reader:
    #         list_all_rates = np.append(list_all_rates,
    #                                    [[row[0], row[1], row[2], row[3], row[4], row[5]]], axis=0)

    all_rates = np.array(list_all_rates[:, 0], dtype=float)
    total_rate = np.sum(all_rates)

    time_step = -(1/total_rate)*math.log(random.random())
    # print(time_step)
    sum_rates = 0.0
    rndm = random.random()
    for k in range(len(list_all_rates)):
        if k != 0:
            sum_rates += float(list_all_rates[k-1, 0])
        lower = (1/total_rate)*sum_rates
        upper = (1/total_rate)*(sum_rates+float(list_all_rates[k, 0]))
        if lower < rndm <= upper:
            selected_event = list(list_all_rates[k])
            # print("Lower: ", lower)
            # print("upper: ", upper)
            # print("Random: ", rndm)
            # print(selected_event[3])
            # print("############")
            break

    # print("selected_event:")
    # print(selected_event)
    # print("##############")

    # NOTE: each element of list_all_rates has the following items: ['rate', 'atom',  'atom.tag', 'target_end_point','target_end_point.tag', 'hopping_mechanism', 'barrier', 'kra']

    selected_rate = selected_event[0]

    selected_atom_tag = int(selected_event[2])
    for atom in atoms:
        if atom.tag == selected_atom_tag:
            selected_atom = atom

    selected_end_tag = int(selected_event[4])
    for atom in atoms:
        if atom.tag == selected_end_tag:
            selected_end = atom
    mechanism = selected_event[5]
    barrier = selected_event[6]
    kra = selected_event[7]
    jump_rate_atom_end_timestp_type_barrier = np.array(
        [selected_rate, selected_atom, selected_end, time_step, mechanism, barrier, kra], dtype=object)

    return jump_rate_atom_end_timestp_type_barrier


# atoms = read(filename='stable_config_initial.xsf')


# def taging_initial_config(atoms):
#     for atom in atoms:
#         atom.tag = atom.index


# taging_initial_config(atoms)
# jump_rate_atom_end_timestp_type = find_jump_rate_atom_end_timestp_type(
#     atoms, 300)

# print(jump_rate_atom_end_timestp_type)
