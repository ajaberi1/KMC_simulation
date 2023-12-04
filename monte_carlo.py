import csv
import json
# from ase.io.trajectory import Trajectory
# from clease.montecarlo.observers import Snapshot
from clease.montecarlo.observers import EnergyEvolution, LowestEnergyStructure
from clease.montecarlo import Montecarlo
from clease.montecarlo.constraints import FixedElement
from clease.calculator import attach_calculator
from clease.settings import CECrystal, Concentration
from ase.io import read, write
# from ase.visualize import view
# from clease.calculator.clease import Clease
from ase import build


######### NOTE: User input #########
target_li_conc = 
temp = 300
supercell_name = 'small_new.xsf'
# supercell_name = 'large.xsf'
####################################


# Specify CE settings

conc = Concentration(basis_elements=[['Li', 'X'], ['Co'], ['O']])
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
    db_name="LiCoO2.db",
    size=[3, 3, 3],
    max_cluster_size=4,
    max_cluster_dia=[7.0, 7.0, 6.0])

# Get the ecis
with open('eci_l2.json') as f:
    eci_file = json.load(f)
eci = eci_file

# NOTE: Method 1: If we do not the have the supercell structure: we have three options: a,b,c

# NOTE: a) make initial larger
# atoms = settings.atoms.copy()*(3, 3, 3)
# atoms = attach_calculator(settings, atoms=atoms, eci=eci)
# view(settings.atoms)

# NOTE: b) make the initial not only larger but also wider
# original_atoms = settings.atoms.copy()*(2, 2, 2)
# p = [[1, 1, -1], [-1, 1, 1], [1, -1, 1]]
# atoms = build.make_supercell(original_atoms, P=p)
# atoms = attach_calculator(settings, atoms=atoms, eci=eci)
# print(new_atoms.get_cell_lengths_and_angles())
# print(new_atoms.get_cell())
# view(settings.atoms)

# NOTE: c) make it from the primitive cell (suggested by alex). It gives the complete layered structure
# atoms = make_supercell(settings.prim_cell.copy(), [[9, -9, 0], [0, 9, -9], [3, 3, 3]])
# atoms = attach_calculator(settings, atoms=atoms, eci=eci)
# write(filename="read_by_clease.xsf", images=atoms)


# NOTE: d) Start from reading the supercell file
atoms = read(filename=supercell_name)
atoms = attach_calculator(settings, atoms=atoms, eci=eci)
write(filename="read_by_clease.xsf", images=atoms)
# # view(settings.atoms)


def count_atoms(atoms):
    """Count the number of each element."""
    atom_count = {key: 0 for key in atoms.symbols}
    for atom in atoms:
        atom_count[atom.symbol] += 1
    return atom_count


count = count_atoms(atoms)


approx_li = count['Li']*target_li_conc
num_li = round(approx_li)
num_x = count['Li']-num_li


def insert_X(num_X):
    counter = 0
    atom_index = 0
    while counter < num_X:
        if atoms[atom_index].symbol == 'Li':
            atoms[atom_index].symbol = 'X'
            counter += 1
        atom_index += 1


insert_X(num_x)
write(filename='random_config_initial.xsf', images=atoms)

# # NOTE: Method 2: if we have the initial supercell structure:
# atoms = read(filename='initial.xsf')
# calc = Clease(settings=settings, eci=eci)
# atoms.set_calculator(calc)

count2 = count_atoms(atoms)
li_sites = count2['Li']+count2['X']
final_li_conc = count2['Li']/(count2['Li']+count2['X'])
print("strable structure propoeries: ")
print("------------------------------")
print("N:  ", li_sites)
print("count X:  ", count2['X'])
print("count Li:  ", count2['Li'])
print("Li concentration", final_li_conc)

steps = 1000*count2['Li']

print("######")
print("######")
print("The number of Monte Carlo steps:", steps)
print("######")
print("######")

# Performing Monte Carlo calculation
# NOTE: For Montecarlo class we need an atoms object which its calculater attached to it.
mc = Montecarlo(atoms, temp)
cnst1 = FixedElement('Co')
mc.generator.add_constraint(cnst1)
cnst2 = FixedElement('O')
mc.generator.add_constraint(cnst2)
obs = EnergyEvolution(mc)
obs2 = LowestEnergyStructure(atoms, verbose=True)
mc.attach(obs, interval=200)
# snap = Snapshot(fname='snapshot', atoms=atoms)
# mc.attach(snap, interval=100)
mc.run(steps=steps)

write(
    filename=f'LowestEnergyStructure_{target_li_conc}.xsf', images=obs2.atoms)
print("Lowest energy: ", obs2.energy)
energies = obs.energies

with open("MC_energies.csv", 'a') as file:
    writer = csv.writer(file)
    writer.writerow(['Row', 'Energy/100 Step'])

with open("MC_energies.csv", 'a') as file:
    writer = csv.writer(file)
    for x in range(len(energies)):
        writer.writerow(
            [x, energies[x]])

# print('row', '  ', 'energies/100 step')
# for x in range(len(energies)):
#     print(x, '  ', energies[x])

thermo = mc.get_thermodynamic_quantities()
print(thermo)
write(filename='stable_config_initial_ensemble_0.xsf', images=atoms)
print("[- Monte Carlo simulation finished - stable structure generated -]")
