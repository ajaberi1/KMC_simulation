import numpy as np
import math
import csv
# from pathlib import Path
# from ase import atom
# from ase import atoms
# from ase.io import read, write


def displacement(atoms, time, jumped_atoms, jump_count, first_row, ensemble, number_of_Li):
    # NOTE: get the displacement of the displaced lithiums
    displaced_final_positions = {}
    n = 0
    # print("List jumped_atoms: ", jumped_atoms)
    for atom in atoms:
        if atom.symbol == 'Li':
            if atom.tag in jumped_atoms:
                n += 1
                # NOTE: the number of jumps (jump_count) does not necessarily equal to 'n' since one atom can be chosen for several jumps - n is the total number of diffusion ion
                displaced_final_positions[str(atom.tag)] = [
                    atom.position[0], atom.position[1], atom.position[2]]
                # print(row[0], row[1], row[2], row[3])
    # print(displaced_final_positions)
    all_initial_positions = {}
    with open(f"initial_ensemble_{ensemble}.csv") as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == 'Li':
                all_initial_positions[row[4]] = [
                    float(row[1]), float(row[2]), float(row[3])]
    # print(all_initial_positions)

    displaced_initial_positions = {}
    for tag in displaced_final_positions:
        displaced_initial_positions[tag] = all_initial_positions[tag]
    # print(displaced_initial_positions)

    length = int(len(displaced_final_positions))
    all_displacement = np.empty((0, length))
    tag_displacment = {}
    for tag in displaced_final_positions:
        displacement = math.sqrt((displaced_final_positions[tag][0] - displaced_initial_positions[tag][0])**2 + (
            displaced_final_positions[tag][1] - displaced_initial_positions[tag][1])**2+(displaced_final_positions[tag][2] - displaced_initial_positions[tag][2])**2)
        all_displacement = np.append(all_displacement, [[displacement]])
        tag_displacment[tag] = displacement

    # NOTE: print the displacement of each lithium atoms at each step time
    row = [jump_count]
    for i in range(len(first_row[1:])):
        row.append(0)

    for tag_atom, displace in tag_displacment.items():
        for i in range(len(row)-1):
            if int(tag_atom) == int(first_row[i+1]):
                row[i+1] = displace

    with open(f"displacments.csv", 'a') as file:
        writer = csv.writer(file)
        writer.writerow(row)

    # print("N: ", n)
    # NOTE: calculate the diffusion coefficient and the correlation factor
    squared_displacement = np.power(all_displacement, 2)
    cor_func1 = ((np.sum(squared_displacement))/n)/(jump_count*(2.84289827**2))
    cor_func2 = ((np.sum(squared_displacement))/n) / \
        ((jump_count/n)*(2.84289827**2))
    sum_disp = np.sum(all_displacement)
    # print("time", time)
    # print("n", n)
    # print("sum_disp", sum_disp)
    diff_coef = (1/(4*time*n))*(sum_disp*1e-8)**2

    # NOTE: At the very begining, there is a chance that after several jumps all the atoms return to their initial positions. That is why we need this if statement
    if diff_coef == 0:
        log_diff_coef = None
        print("Reached the initial positions...can NOT calculate diffusion coef")
    else:
        log_diff_coef = math.log10(diff_coef)

    diff_coef_new = (1/(4*time*number_of_Li))*(sum_disp*1e-8)**2
    if diff_coef_new == 0:
        log_diff_coef_new = None
    else:
        log_diff_coef_new = math.log10(diff_coef_new)

    return log_diff_coef, cor_func1, cor_func2, log_diff_coef_new
