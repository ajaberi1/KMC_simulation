from local_atoms import find_local_atoms
import numpy as np
from clease import NewStructures
from clease.calculator import attach_calculator
from clease.tools import update_db, wrap_and_sort_by_position
from clease.settings import  Concentration, ClusterExpansionSettings
from clease.calculator.clease import Clease
from ase.calculators.espresso import Espresso
import pandas as pd
import json
from ase.io import read, write
from nearest_neighbor import find_nearest_neighbor, find_neighbor_from_DataFrame
from jump_rates import find_rates
import random
import csv
import os
import sys
from displacement import displacement
from select_jump import find_jump_rate_atom_end_timestp_type_barrier
from monte_carlo_implemented import monte_carlo
import time as t
from clease.montecarlo import Montecarlo
from clease.montecarlo.constraints import FixedElement
import subprocess
from joblib import Parallel, delayed
import concurrent.futures
import multiprocessing as mp
import math
# from jump_rates import calculate_energy
# from nearest_neighbor import find_nearest_neighbor_old
# from compare_lists_all_rates import compare
# from perform_jumps import perform_jump
# from ase.db import connect
# from clease.tools import update_db
# from ase import Atoms
# from ase import atoms
# from ase import atom
# import math
# from clease.calculator import attach_calculator, get_ce_energy

######### NOTE: User input #########
target_li_conc = 0.2
number_of_ensembles = 10
steps_MC = 5000
temp = 300
project_name = "KMC_V3.6.3"
max_num_KMC = 500
max_num_jumps = 100000
parallelization = False
use_C = False
####################################

def count_atoms(atoms):
    atom_count = {key: 0 for key in atoms.symbols}
    for atom in atoms:
        atom_count[atom.symbol] += 1
    return atom_count


def taging_initial_config(atoms):
    for atom in atoms:
        atom.tag = atom.index


# Specify CE settings
conc = Concentration(basis_elements=[['Co'], ['Co'], ['Co'], ['Li', 'X'], ['Li', 'X'], ['Li', 'X'], ['Li', 'X'], ['Li', 'X'], ['Li', 'X'], ['Li', 'X'], ['Li', 'X'], ['Li', 'X'], [
                     'Mn'], ['Mn'], ['Mn'], ['Ni'], ['Ni'], ['Ni'], ['O'], ['O'], ['O'], ['O'], ['O'], ['O'], ['O'], ['O'], ['O'], ['O'], ['O'], ['O'], ['O'], ['O'], ['O'], ['O'], ['O'], ['O']], grouped_basis=[[0, 1, 2], [3, 4, 5, 6, 7, 8, 9, 10, 11], [12, 13, 14], [15, 16, 17], [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]])
conc.set_conc_formula_unit(
    formulas=["Co<1>", "Li<x>X<1-x>", "Mn<1>", "Ni<1>", "O<1>"], variable_range={"x": (0, 1)})

atoms = wrap_and_sort_by_position(read("NMC.xsf"))

### NOTE: when we have p1 spacegroup (not worry about the symmetry and don't want to use the samllest primitive cell), we use this option to define setting for CLEASE ###
settings = ClusterExpansionSettings(atoms, conc, max_cluster_dia=[7, 7, 7], db_name="NMC.db",
                                    size=[2, 2, 1])

# Get the ecis
with open('eci_l2.json') as f:
    eci_file = json.load(f)
eci = eci_file

atoms = settings.atoms.copy()*(3, 3, 3)
atoms = attach_calculator(settings, atoms=atoms, eci=eci)

calc = Clease(settings=settings, eci=eci)
atoms.set_calculator(calc)

######### NOTE: Initializing the structure #########
print("...Making the fisrt stable structure with Monte Carlo simulation...")
monte_carlo(atoms, temp, target_li_conc)
write(filename='stable_config_initial_ensemble_0.xsf', images=atoms)

taging_initial_config(atoms)
print("Initialization...taging done...")

if use_C == True:
    print("Compiling the C++ code..")
    os.system("g++ main.cpp -o neighbors_by_C")

count = count_atoms(atoms)
number_of_Li = count['Li']
number_of_X = count['X']
li_conc = (number_of_Li)/(number_of_Li+number_of_X)
a_file = open("structure_info.txt", "w")
print("Number of Li: ", number_of_Li, file=a_file)
print("Number of vacancies: ", number_of_X, file=a_file)
print("Li concentration: ", li_conc, file=a_file)
a_file.close()

# This part store all the position of all the transition metals
transition_metals_list = []
for atom in atoms:
    if atom.symbol == "Mn" or atom.symbol == "Ni" or atom.symbol == "Co":
        transition_metals_list.append([atom.symbol, atom.position[0], atom.position[1], atom.position[2]])

transition_metals_np = np.empty((0, 4))
for i in range(len(transition_metals_list)):
    # print(i, transition_metals_list[i])
    transition_metals_np = np.vstack([transition_metals_np, transition_metals_list[i]])

print("")
print(f"[-- start KMC for Li conc: {round(li_conc,2)} --]")
print("")

# display(all_neighbors_DataFrane)

# def make_copies():
#     #for parrallelization
#     global atoms_copy
#     global atoms_copy_calcs
#     new_atoms = atoms.copy()
#     copy_calc = Clease(settings=settings, eci=eci)
#     new_atoms.set_calculator(copy_calc)
#     atoms_copy_calcs.append(copy_calc)
#     atoms_copy.append(new_atoms)

ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK', default=1))
if parallelization:
    print(f"______Run in parallel: {ncpus} CPUs______")
    print(f"making Atoms copies...")
    atoms_copy = []
    atoms_copy_calcs = []
    for i in range(ncpus):
        new_atoms = atoms.copy()
        copy_calc = Clease(settings=settings, eci=eci)
        new_atoms.set_calculator(copy_calc)
        atoms_copy_calcs.append(copy_calc)
        atoms_copy.append(new_atoms)
else:
    print("______Run with Unparallel option______")

######################################################################################################################################################################################################
######################################################################################################################################################################################################
# This section contains all the fuctions required for KMC loops
# NOTE: temp and num_KMCS come from user input
######################################################################################################################################################################################################
######################################################################################################################################################################################################

# def change_atoms(atoms, system_changes):
#     atoms[system_changes[0][0]].symbol = system_changes[0][2]
#     atoms[system_changes[1][0]].symbol = system_changes[1][2]

#     atoms[system_changes[0][0]].tag, atoms[system_changes[1][0]
#                                            ].tag = atoms[system_changes[1][0]].tag, atoms[system_changes[0][0]].tag


def update_neighbors_DataFrane(atoms, swept_atoms):
    # NOTE: Since Atoms object is a globlen object, you do not need to update neighbors DataFrame for all atoms around the jumped atom. Only the columns of the swept atoms shuold be updated.
    for atom in swept_atoms:
        nearest_atom_object = find_nearest_neighbor(atoms, atom)
        nearest_atom_object.extend([None]*(6-len(nearest_atom_object)))
        all_neighbors_DataFrane[atom.tag] = nearest_atom_object


def exchange_tags(atoms, system_changes):
    atoms[system_changes[0][0]].tag, atoms[system_changes[1][0]
                                           ].tag = atoms[system_changes[1][0]].tag, atoms[system_changes[0][0]].tag
    return atoms


def perform_jump(atoms, temp, list_all_rates, calc, parallelization):
    selected_jump = find_jump_rate_atom_end_timestp_type_barrier(
        atoms, temp, list_all_rates)

    # rate = selected_jump[0]
    selected_atom = selected_jump[1]
    end_point = selected_jump[2]
    jump_time = selected_jump[3]
    jump_type = selected_jump[4]
    barrier = selected_jump[5]
    kra = selected_jump[6]
    t_metal_type = selected_jump[7]

    system_changes = [(selected_atom.index, selected_atom.symbol, end_point.symbol),
                      (end_point.index, end_point.symbol, selected_atom.symbol)]

    jumped_atom_tag = selected_atom.tag

    # jumped_atom_tag = int(atoms[system_changes[0][0]].tag)
    # print("jumped_atom_tag", jumped_atom_tag)
    # print('jumped atom: ', selected_atom)
    # print('end: ', end_point)
    # print('jump type: ', jump_type)

    exchange_tags(atoms, system_changes)
    new_energy = calc.get_energy_given_change(
        system_changes, keep_changes=True)
    # write(filename=f"structure_{system_changes[0][0]}.xsf", images=atoms)

    # NOTE: note that after the jump, the formerly jumped atom is now in the endpoint (jumped_atom=end_point). The order of items in swept_jump_list does not matter.
    swept_atoms = [end_point, selected_atom]
    update_neighbors_DataFrane(atoms, swept_atoms)

    if parallelization:
        global atoms_copy
        global atoms_copy_calcs
        for i in range(len(atoms_copy)):
            exchange_tags(atoms_copy[i], system_changes)
            atoms_copy_calcs[i].get_energy_given_change(
                system_changes, keep_changes=True)
            # atoms_copy[i].get_calculator().clear_history()
            # change_atoms(atoms_item, system_changes)

        # for i in range(len(atoms)):
        #     if atoms[i].tag != atoms_copy[0][i].tag:
        #         print("@@@@@@@@@@")
        #         print("@@@@@@@@@@")
        #         print("@@@@@@@@@@")
        #         print("Something is wrong in atoms global variable")
        #         print("@@@@@@@@@@")
        #         print("@@@@@@@@@@")
        #         print("@@@@@@@@@@")

    # NOTE: from now on the atoms object is completely changed everywhere including the neigbors DataFrame because it is a class level atrribiute of Atoms class
    # print(new_energy, jump_time, jumped_atom_tag, jump_type)
    one_jump = [new_energy, jump_time,
                jumped_atom_tag, jump_type, barrier, kra, t_metal_type]
    return one_jump


def get_rates(tag_iteration):
    # NOTE: This function will run in parrallel
    # NOTE: We will get the rates with the copy of atoms objects not the original atoms object

    # print("I'm here 4")

    atom_tag = tag_iteration[0]
    iteration = tag_iteration[1]
    global atoms_copy
    new_atoms = atoms_copy[iteration-1]

    for atom in new_atoms:
        if atom.tag == atom_tag:
            my_atom = atom
            break

    with open("energy.csv") as file:
        reader = csv.reader(file)
        readers_list = list(reader)
        initial_position_energy = float(readers_list[-1][0])
        # print("initial_position_energy for getting all_rats: ",initial_position_energy)

    global temp
    global all_neighbors_DataFrane
    possible_jump = []
    nearest_neigbors = find_neighbor_from_DataFrame(
        all_neighbors_DataFrane, my_atom)

    end_point = np.empty((0, 1))
    for atom1 in nearest_neigbors:
        if atom1.symbol == "X":
            end_point = np.append(end_point, atom1)
    if len(end_point) == 0:
        pass
    else:
        for target_end_point in end_point:
            hopping_mechanism, rate, barrier, kra, t_metal = find_rates(
                new_atoms, nearest_neigbors, my_atom, target_end_point, initial_position_energy, temp, all_neighbors_DataFrane, transition_metals_np)
            possible_jump.append(
                [rate, my_atom, target_end_point, hopping_mechanism, barrier, kra, t_metal])

    possible_jump_file_suffix = None
    # For cases where all neighbors are filled with Li,  len(possible_jump) equals 0: thus we need this following if statement:
    if len(possible_jump) != 0:
        with open(f"possible_jump_{atom_tag}.csv", 'w') as file:
            possible_jump_file_suffix = atom_tag
            for item in possible_jump:
                writer = csv.writer(file)
                writer.writerow([item[0], item[1], item[1].tag,
                                 item[2], item[2].tag, item[3], item[4], item[5]])
        # else:
        #     print(
        #         f"got a possible_jump with 0 length for iteration: {iteration}")

    # print("I'm here 4-2")
    # print("possible_jump_file_suffix: ",possible_jump_file_suffix)

    return possible_jump_file_suffix


# NOTE: This function will get the rates of all the atoms that were passed in the first parmeter...it can be eigther the atoms-object or a list of atom
def all_rates(atoms, temp, initial_position_energy, parallelization, atoms_list):
    # print("I am here 2")
    list_all_rates = np.zeros((0, 9))

    # NOTE: No parallelization - method 0
    if not parallelization:
        for atom in atoms_list:
            if atom.symbol == "Li":
                selected_atom = atom
                nearest_neigbors = find_neighbor_from_DataFrame(
                    all_neighbors_DataFrane, selected_atom)
                end_point = np.empty((0, 1))
                for atom1 in nearest_neigbors:
                    if atom1.symbol == "X":
                        end_point = np.append(end_point, atom1)
                # print("len(end_point)= ", len(end_point))
                if len(end_point) == 0:
                    pass
                else:
                    for target_end_point in end_point:
                        hopping_mechanism, rate, barrier, kra, t_metal = find_rates(atoms, nearest_neigbors, selected_atom, target_end_point, initial_position_energy, temp, all_neighbors_DataFrane,transition_metals_np)
                        list_all_rates = np.append(
                            list_all_rates, [[rate, atom,  atom.tag, target_end_point, target_end_point.tag, hopping_mechanism, barrier, kra, t_metal]], axis=0)

    if parallelization:

        with open("energy.csv", 'a') as file:
            writer = csv.writer(file)
            writer.writerow(
                [initial_position_energy])

        # NOTE: This 'paralelization_groups' paramater indicates the Li tag and its parrallel number (iteration). The parrallel number is used to indicate which one of the atoms_copy object is being used for each cpu.
        paralelization_groups = []
        iteration = 0
        each_group = []
        for atom in atoms_list:
            if atom.symbol == "Li":
                iteration += 1
                each_group.append([int(atom.tag), iteration])
                if iteration == ncpus:
                    paralelization_groups.append(each_group)
                    iteration = 0
                    each_group = []
        if iteration != 0:
            paralelization_groups.append(each_group)

        all_suffixes = []
        possible_jump_file_all_suffixes = []

        # NOTE: parallelization - method 1
        # my_time = t.time()
        for each_group in paralelization_groups:
            with concurrent.futures.ProcessPoolExecutor(max_workers=ncpus) as executor:
                possible_jump_file_all_suffixes = []
                results = [executor.submit(get_rates, tag_iteration)
                           for tag_iteration in each_group]
                for f in concurrent.futures.as_completed(results, timeout=None):
                    a = f.result()
                    possible_jump_file_all_suffixes.append(a)
            for sufix in possible_jump_file_all_suffixes:
                all_suffixes.append(sufix)
        # print("parallel time: ", t.time() - my_time)
        # results = executor.map(get_rates, data) #Use this one if the sequence of outputs matters
        # for result in results:
        #     possible_jump_file_all_suffixes.append(result)

        # NOTE: parallelization - method 2
        # pool = mp.Pool(processes=ncpus)
        # for each_group in paralelization_groups:
        #     results = [pool.apply_async(get_rates, args=(
        #         tag_iteration,)) for tag_iteration in each_group]
        #     possible_jump_file_all_suffixes = [p.get() for p in results]
        #     for sufix in possible_jump_file_all_suffixes:
        #         all_suffixes.append(sufix)

        # NOTE: parallelization - method 3
        # for each_group in paralelization_groups:
        #     results = Parallel(n_jobs=ncpus, backend='multiprocessing')(delayed(get_rates)(tag_iteration)
        #                                                                 for tag_iteration in each_group)
        #     for result in results:
        #         possible_jump_file_all_suffixes.append(result)

        #     for sufix in possible_jump_file_all_suffixes:
        #         all_suffixes.append(sufix)

        # After parallelization we have to read the list_all_rates from all the possible_jump_* files
        for suffix in all_suffixes:
            if suffix != None:
                with open(f"possible_jump_{suffix}.csv") as file:
                    reader = csv.reader(file)
                    for row in reader:
                        list_all_rates = np.append(list_all_rates,
                                                   [[row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7]]], axis=0)

    return list_all_rates


def update_all_rates(atoms, temp, initial_position_energy, jumped_atom_tag, list_all_rates, parallelization):
    new_all_rates = np.zeros((0, 10))
    # print("jumped_atom_tag: ", jumped_atom_tag)
    for atom in atoms:
        if atom.tag == jumped_atom_tag:
            jumped_atom = atom
            break

    # NOTE: Option1: update the rates of large number of lithium atoms in the same plane of the jumped atom and also the rates of lithiums at the adjacant planes CLOSE to the jumped atom
    # upading_atoms = []
    # for atom in atoms:
    #     if atom.symbol == "Li":
    #         if abs(atom.position[2] - jumped_atom.position[2]) < 0.001:
    #             upading_atoms.append(atom)

    # # Find tha atoms close to jumped-atom at the planes just above and below the current plane
    # local_atoms_of_jumped_atom = find_local_atoms(
    #     atoms, jumped_atom, local_atom_type='Li', threshold_distance=10)
    # for atom in local_atoms_of_jumped_atom:
    #     if abs(atom.position[2] - jumped_atom.position[2]) > 0.001:
    #         upading_atoms.append(atom)

    # # NOTE: Option2: update the rates of all local atoms with the same distance from the jumped atom
    upading_atoms = [jumped_atom]
    local_atoms_of_jumped_atom = find_local_atoms(
        atoms, jumped_atom, local_atom_type='Li', threshold_distance=10)
    for atom1 in local_atoms_of_jumped_atom:
        upading_atoms.append(atom1)

    # Fore checking purposs:
    # all_indexes = []
    # for item in upading_atoms:
    #     all_indexes.append(item.tag)
    # print("upading_atoms tag:", all_indexes)
    # print("number of upading_atoms :", len(upading_atoms))

    new_all_rates = all_rates(
        atoms, temp, initial_position_energy, parallelization, upading_atoms)

    row_to_delete = []
    for atom in upading_atoms:
        row = np.where(np.array(list_all_rates[:, 2]) == atom.tag)
        row_to_delete.extend(list(row[0]))
    # print("row_to_delete:", row_to_delete)
    reduced_list_all_rates = np.delete(
        list_all_rates, row_to_delete, axis=0)

    new_list_all_rates = np.concatenate(
        (reduced_list_all_rates, new_all_rates), axis=0)

    return new_list_all_rates


##################################################
##################################################
################ NOTE KMC Loops ##################
##################################################
##################################################
# can not be displayed: Copy Right
