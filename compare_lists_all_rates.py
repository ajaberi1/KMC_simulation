import numpy as np
from clease import NewStructures
from clease import CEBulk
from clease import CECrystal
from clease import Concentration
from ase.calculators.espresso import Espresso
from ase.db import connect
from clease.tools import update_db
from ase import Atom
from ase import atoms
from ase import atom
import math
from ase.io.trajectory import Trajectory
from clease.montecarlo.observers import Snapshot
from clease.montecarlo.observers import EnergyEvolution
from clease.montecarlo import Montecarlo
from clease.montecarlo.constraints import FixedElement
from clease.calculator import attach_calculator
import json
from ase.io import read, write
from pathlib import Path
import csv


def compare(list_all_rates, original_list_all_rates, jump_count):
    first_last_barriers = []
    # you just need two fors, , , no nned while

    for item in original_list_all_rates:
        for equivalent_item in list_all_rates:
            if item[2] == equivalent_item[2]:
                if item[4] == equivalent_item[4]:
                    first_last_barriers.append(
                        [item[2], item[4], item[6], equivalent_item[2], equivalent_item[4], equivalent_item[6]])
                    break

        # if counter2 < len()

    with open(f"original_vs_update_list_{jump_count}.csv", 'w') as file:
        writer = csv.writer(file)
        writer.writerow(['update barrier', 'update first', 'update last',
                         'original barrier', 'original first', 'original last', 'barrier_difference'])

    with open(f"original_vs_update_list_{jump_count}.csv", 'a') as file:
        writer = csv.writer(file)
        for item in first_last_barriers:
            writer.writerow(
                [item[0], item[1], item[2], item[3], item[4], item[5], abs(float(item[2]) - float(item[5]))])

    barrier_difference = []
    for item in first_last_barriers:
        difference = abs(float(item[2]) - float(item[5]))
        barrier_difference.append(difference)

    barrier_difference.sort(reverse=True)
    print("first five largest barrier_difference: ", barrier_difference[0:5])

    for i in barrier_difference:
        if i > 0.1:

            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("Large Diference found")
            print("Large Diference found")
            print("Large Diference found")
            print("Large Diference found")
            print("Large Diference found")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")
            print("#####################")

        break


# # NOTE: use this part to copare it manually
# jump_count = 1

# list_all_rates = np.empty((0, 10))
# with open(f"list_all_rates_{jump_count}.csv") as f:
#     reader = csv.reader(f)
#     for row in reader:
#         list_all_rates = np.append(
#             list_all_rates, np.array([[row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9]]]), axis=0)

# original_list_all_rates = np.empty((0, 10))
# with open(f"original_list_all_rates_{jump_count}.csv") as f:
#     reader = csv.reader(f)
#     for row in reader:
#         original_list_all_rates = np.append(
#             original_list_all_rates, np.array([[row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9]]]), axis=0)

# compare(list_all_rates, original_list_all_rates, jump_count)
