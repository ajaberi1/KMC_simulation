from local_atoms import find_local_atoms
import numpy as np
from clease import NewStructures
from clease.settings import CECrystal, Concentration
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
import time as t
from clease.montecarlo import Montecarlo
from clease.montecarlo.constraints import FixedElement
import subprocess
from joblib import Parallel, delayed
import concurrent.futures
import multiprocessing as mp
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
number_of_ensembles = 10
steps_MC = 5000
temp = 300
project_name = "KMC_V3.6"
max_num_KMC = 500
max_num_jumps = 100000
parallelization = False
####################################

print("...Making the fisrt stable structure with Monte Carlo simulation...")
complete = subprocess.run(["python", "monte_carlo.py"],
                          capture_output=True, text=True)

print("stdout", complete.stdout)
atoms = read(filename='stable_config_initial_ensemble_0.xsf')
# atoms = read(filename='stable_config_initial.xsf')

######### NOTE: Initializing the structure #########


def count_atoms(atoms):
    atom_count = {key: 0 for key in atoms.symbols}
    for atom in atoms:
        atom_count[atom.symbol] += 1
    return atom_count


def taging_initial_config(atoms):
    for atom in atoms:
        atom.tag = atom.index


taging_initial_config(atoms)
print("Initialization...taging done...")


count = count_atoms(atoms)
number_of_Li = count['Li']
number_of_X = count['X']
li_conc = (number_of_Li)/(number_of_Li+number_of_X)
a_file = open("structure_info.txt", "w")
print("Number of Li: ", number_of_Li, file=a_file)
print("Number of vacancies: ", number_of_X, file=a_file)
print("Li concentration: ", li_conc, file=a_file)
a_file.close()

print("")
print(f"[-- start KMC for Li conc: {round(li_conc,2)} --]")
print("")

conc = Concentration(basis_elements=[['Li', 'X'], ['Co'], ['O']])


# Specify CE settings
settings = CECrystal(
    concentration=conc,
    spacegroup=166,
    basis=[
        (0.00000000,        0.00000000,        0.00000000),
        (-0.00000000,        -0.00000000,        0.50000000),
        (0.00000000,        0.00000000,        0.23958700)
    ],
    cellpar=[2.84289827, 2.84289827, 14.14561550,
             90.00000000, 90.00000000, 120.00000000],
    size=[1, 1, 1],
    db_name="LiCoO2.db",
    max_cluster_size=4,
    max_cluster_dia=[7.0, 7.0, 6.0])


# Get the ecis
with open('eci_l2.json') as f:
    eci_file = json.load(f)
eci = eci_file

calc = Clease(settings=settings, eci=eci)
atoms.set_calculator(calc)


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
                jumped_atom_tag, jump_type, barrier, kra]
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
            hopping_mechanism, rate, barrier, kra = find_rates(
                new_atoms, nearest_neigbors, my_atom, target_end_point, initial_position_energy, temp, all_neighbors_DataFrane)
            possible_jump.append(
                [rate, my_atom, target_end_point, hopping_mechanism, barrier, kra])

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
    list_all_rates = np.zeros((0, 8))

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
                if len(end_point) == 0:
                    pass
                else:
                    for target_end_point in end_point:
                        hopping_mechanism, rate, barrier, kra = find_rates(
                            atoms, nearest_neigbors, selected_atom, target_end_point, initial_position_energy, temp, all_neighbors_DataFrane)
                        list_all_rates = np.append(
                            list_all_rates, [[rate, atom,  atom.tag, target_end_point, target_end_point.tag, hopping_mechanism, barrier, kra]], axis=0)

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
        row = np.where(np.array(list_all_rates[:,2])==atom.tag)
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

for i in range(number_of_ensembles):
    os.system(f"mkdir ensemble_{i}")

for ensemble in range(number_of_ensembles):
    start_time = t.time()
    print(f"Performing KMC on ensemble number {ensemble} .....")
    if ensemble == 0:
        pass
    else:
        if parallelization:
            del atoms_copy
            del atoms_copy_calcs

        mc = Montecarlo(atoms, temp)
        cnst1 = FixedElement('Co')
        mc.generator.add_constraint(cnst1)
        cnst2 = FixedElement('O')
        mc.generator.add_constraint(cnst2)
        mc.run(steps=steps_MC)
        write(
            filename=f"stable_config_initial_ensemble_{ensemble}.xsf", images=atoms)

        taging_initial_config(atoms)
        print("Initialization...taging done...")

        atoms_copy = []
        atoms_copy_calcs = []
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

    arr = np.empty((0, 6))
    for atom in atoms:
        arr = np.append(
            arr, [[atom.symbol, atom.position[0], atom.position[1], atom.position[2], atom.tag, atom.index]], axis=0)

    with open(f"initial_ensemble_{ensemble}.csv", "w") as file:
        writer = csv.writer(file)
        for row in arr:
            writer.writerow(row)

    with open(f"results_ensemble_{ensemble}.csv", 'a') as file:
        writer = csv.writer(file)
        writer.writerow(['KMCS', 'jump_count', 'time',
                         'energy', 'log_D', 'mechanism', 'cor_func1', 'cor_func2', 'barrier', 'kra', 'log_D_new'])

    first_row = ['jumps']
    for atom in atoms:
        if atom.symbol == "Li":
            first_row.append(atom.tag)

    with open(f"displacments.csv", 'w') as file:
        writer = csv.writer(file)
        writer.writerow(first_row)

    initial_time0 = t.time()
    print(f"Making neighbors DataFrame for ensemble_{ensemble}...")
    all_neighbors_dict = {}
    for atom in atoms:
        if atom.symbol == 'Li' or atom.symbol == 'X':
            nearest_atom_object = find_nearest_neighbor(atoms, atom)
            nearest_atom_object.extend([None]*(6-len(nearest_atom_object)))
            all_neighbors_dict[atom.tag] = nearest_atom_object
    all_neighbors_DataFrane = pd.DataFrame(all_neighbors_dict)
    print("Time to get DataFrame: ", t.time() - initial_time0)
    # pd.set_option('max_columns', None)
    # pd.set_option('display.max_colwidth',500)
    # display(all_neighbors_DataFrane)

    system_changes = None
    new_energy = calc.get_energy_given_change(system_changes)
    # new_energy = calc.calculate()
    # new_energy = get_ce_energy(settings, atoms=atoms, eci=eci)

    # because this is the first time we wanna calculate the 'list_all_rates', all the 'atoms' object should be pass in the 'atoms_list'
    print("Making the intial list_all_rates...")
    initial_time = t.time()
    list_all_rates = all_rates(
        atoms, temp, new_energy, parallelization, atoms_list=atoms)
    print("Time to get the intial list_all_rates: ", t.time() - initial_time)

    # print("I am here 7")
    with open(f"original_list_all_rates_ensemble_{ensemble}.csv", 'w') as file:
        writer = csv.writer(file)
        writer.writerow(['rate', 'atom',  'atom.tag', 'target_end_point',
                         'target_end_point.tag', 'hopping_mechanism', 'barrier', 'kra'])
        for item in list_all_rates:
            writer.writerow(
                [float(item[0]), item[1], item[2], item[3], item[4], item[5], item[6], item[7]])

    ############################
    ############################
    ### NOTE: LOOP for jumps ###
    ############################
    ############################
    kmcs = 0
    time = 0
    jump_count = 0
    all_time_steps = []
    list_jumped_atoms = []
    count_jumped_type = {'monovacancy': 0, 'divacancy': 0}
    start_time_of_KMC_step = t.time()
    while True:
        one_jump = perform_jump(
            atoms, temp, list_all_rates, calc, parallelization)
        jump_count += 1
        new_energy = one_jump[0]
        time += one_jump[1]
        jumped_atom_tag = one_jump[2]
        # print("jumped_atom_tag", jumped_atom_tag)
        jump_type = one_jump[3]
        list_jumped_atoms.append(jumped_atom_tag)
        barrier = one_jump[4]
        kra = one_jump[5]
        count_jumped_type[jump_type] += 1
        # print("barrier: ", barrier)
        # print("new energy: ", new_energy)
        # print("list of jumped_atoms: ", list_jumped_atoms)
        # atoms.get_calculator().clear_history()
        # print("jumped_atoms: ", jumped_atoms)
        # print("energy", new_energy)
        # print(f"perfomed the {x}th step")
        # print("#################.....performed one jump.....#################")

        # for atom in atoms:
        #     if atom.tag == jumped_atom_tag:
        #         atoms[atom.index].symbol = 'Cs'
        #         write(filename=f'jump_{jump_count}.xsf', images=atoms)
        #         atoms[atom.index].symbol = 'Li'
        updating_time = t.time()
        new_list_all_rates = update_all_rates(
            atoms, temp, new_energy, list_jumped_atoms[-1], list_all_rates, parallelization)
        if jump_count < 11:
            timing = t.time() - updating_time
            print("Time to UPDATE list_all_rates (for developers): ", timing)
            all_time_steps.append(timing)
        if jump_count == 11:
            all_timing = np.array([timing])
            time_average = np.average(all_timing)
            print(
                "Average time for updating first 10 list_all_rates (for developers): ", time_average)
            print("Continue performing the jumps...")

        list_all_rates = new_list_all_rates
        # NOTE: This is to see how list_all_rate changes after some steps
        if jump_count % 100000 == 0:
            # if jump_count % 1 == 0:
            with open(f"list_all_rates_{jump_count}.csv", 'w') as file:
                writer = csv.writer(file)
                for item in list_all_rates:
                    writer.writerow(
                        [float(item[0]), item[1], item[2], item[3], item[4], item[5], item[6], item[7]])
            os.system(
                f"mv list_all_rates_{jump_count}.csv ensemble_{ensemble}")
        # ##### end of NOTE ######

        # # NOTE: calculate the list_all_rates from scratch
        # original_list_all_rates = all_rates(
        #     atoms, temp, new_energy, parallelization)
        # if jump_count % 1 == 0:
        #     with open(f"original_list_all_rates_{jump_count}.csv", 'w') as file:
        #         writer = csv.writer(file)
        #         for item in original_list_all_rates:
        #             writer.writerow(
        #                 [float(item[0]), item[1], item[2], item[3], item[4], item[5], item[6], item[7], item[8], item[9]])

        #     compare(list_all_rates, original_list_all_rates, jump_count)
        if jump_count % 5 == 0:

            log_D, cor_func1, cor_func2, log_D_new = displacement(
                atoms, time, set(list_jumped_atoms), jump_count, first_row, ensemble, number_of_Li)

            with open(f"results_ensemble_{ensemble}.csv", 'a') as file:
                writer = csv.writer(file)
                writer.writerow(
                    [kmcs, jump_count, time, new_energy, log_D, count_jumped_type, cor_func1, cor_func2, barrier, kra, log_D_new])
            # write(filename=f"structure_final_{ensemble}.xsf", images=atoms)

        if jump_count % number_of_Li == 0:
            print(
                f"Performed {kmcs} KMC steps for {t.time() - start_time_of_KMC_step} seconds")
            start_time_of_KMC_step = t.time()
            kmcs += 1

        if kmcs == max_num_KMC or jump_count == max_num_jumps:
            print("Reached maximum number of KMC/jumps!")
            print("kmcs:", kmcs)
            print("jump_count:", jump_count)
            print(" ")
            print(f"writing the final_ensemble_{ensemble}.csv and structure_final_{ensemble}.xsf files... ")
            arr_final = np.empty((0, 6), dtype=object)
            for atom in atoms:
                arr_final = np.append(
                    arr_final, [[atom.symbol, atom.position[0], atom.position[1], atom.position[2], atom.tag, atom.index]], axis=0)

            with open(f"final_ensemble_{ensemble}.csv", "w") as file:
                writer = csv.writer(file)
                for row in arr_final:
                    writer.writerow(row)

            write(filename=f'structure_final_{ensemble}.xsf', images=atoms)

            with open(f'1_completed_ensemble_{ensemble}.txt', 'w') as f:
                print(f"Done with {jump_count} jumps", file=f)

            if parallelization:
                os.system(f"rm energy.csv")
            os.system(f"mv structure_final_{ensemble}.xsf ensemble_{ensemble}")
            os.system(
                f"mv initial_ensemble_{ensemble}.csv ensemble_{ensemble}")
            os.system(
                f"mv stable_config_initial_ensemble_{ensemble}.xsf ensemble_{ensemble}")
            os.system(
                f"mv results_ensemble_{ensemble}.csv ensemble_{ensemble}")
            os.system(
                f"mv displacments.csv ensemble_{ensemble}")
            os.system(f"mv final_ensemble_{ensemble}.csv ensemble_{ensemble}")
            os.system(
                f"mv original_list_all_rates_ensemble_{ensemble}.csv ensemble_{ensemble}")

            print("final time (t parameter): ", time)
            with open('time.txt', 'w') as f:
                print(time, file=f)
            required_time = t.time() - start_time
            print(
                f"#################.....KMC of ensemble {ensemble} completed with {required_time} seconds.....#################")
            break


with open(f'All_COMPLETED.txt', 'w') as f:
    print(f"Now run get_results.py file to get all the results", file=f)
print("")
print("")
print("##########")
print("...Done...")
print("##########")
print("")
print("END OF CODE")

# complete2 = subprocess.run(["python", "mail.py"],
#                            capture_output=True, text=True)

# print("stdout", complete2.stdout)

# os.system("python mail.py")
