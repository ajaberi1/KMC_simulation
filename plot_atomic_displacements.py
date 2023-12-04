import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
from pathlib import Path
import os
import matplotlib.gridspec as gridspec

num_rows = 7

os.system(f"mkdir Li_dist_plots")
df = pd.read_csv(f"displacments.csv")
# print(df.head())
# print(df.columns[1:])
x = df["jumps"]
zeros = [0 for i in range(len(x))]
column_number = 0
while column_number < len(df.columns[1:]):
    all_ys = []
    group_count = 0
    for i in range(num_rows):
        all_ys.append(df[df.columns[column_number+1]])
        column_number += 1
        group_count += 1
        if column_number == len(df.columns[1:]):
            # print("break perfomed")
            break
    my_counter = group_count
    while my_counter != num_rows:
        all_ys.append(0)
        my_counter += 1
    counter = 0
    # print("group_count", group_count)
    fig, axs = plt.subplots(len(all_ys), 1, sharex=True,
                            constrained_layout=True, figsize=(8, 15), dpi=200)
    for i in range(group_count):
        if group_count != 1:
            axs[i].plot(x, all_ys[counter])
            axs[i].set_yticks(np.arange(0, 25, 10))
            axs[i].set_ylim(-1, 25)
            axs[i].set_title(
                f"{df.columns[column_number-group_count+i+1]}", y=1.0, pad=-14)
            counter += 1
        if group_count == 1:
            axs[0].plot(x, all_ys[counter])
            axs[0].set_yticks(np.arange(0, 25, 10))
            axs[0].set_ylim(-1, 25)
            axs[0].set_title(
                f"{df.columns[column_number-group_count+i+1]}", y=1.0, pad=-14)
            counter += 1
            break
    while counter != num_rows:
        axs[counter-1].plot(x, zeros)
        axs[counter-1].set_yticks(np.arange(0, 25, 10))
        axs[counter-1].set_ylim(-1, 25)
        axs[counter-1].set_title(
            f"{df.columns[column_number-group_count+i+1]}", y=1.0, pad=-14)
        counter += 1
    plt.savefig(f"displ_{column_number}.jpg")
    os.system(f"mv displ_{column_number}.jpg Li_dist_plots")
    # plt.show()

print("...DONE...")
