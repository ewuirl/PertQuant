# This script makes a bar plot showing the random split means and standard deviations of the MAE across the seven datasets. 
# Each dataset had five ANNs trained on a different train-val set.

# importing package
import matplotlib.pyplot as plt
import numpy as np

def figure_formatting(): # make the plot look xmgrace-esque
    font = {'family': 'sans-serif', 'weight': 'bold', 'size': 18}
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rc('font', **font)
    plt.rcParams['axes.linewidth'] = 2.5
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 2.5
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 2.5
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['mathtext.default'] = 'regular'
    plt.rcParams['legend.fancybox'] = False # no rounded legend box

figure_formatting()

# create data. Using test performance
# Order is: 2-2-a, 4-4-a, 10-10-a, 20-20-a, 2-2-sa, 4-4-sa, 10-10-sa 
C1_means = [9.72E-09, 6.06E-09, 2.70E-09, 1.41E-09, 8.00E-09, 5.94E-09, 7.10E-09]
C1_stdevs = [2.25052E-09, 4.98095E-09, 1.15211E-09, 6.73589E-10, 3.27101E-09, 1.33147E-09, 2.49932E-09]
C2_means = [1.67E-08, 5.70E-09, 2.68E-09, 2.66E-09, 7.90E-09, 5.12E-09, 6.13E-09]
C2_stdevs = [4.3908E-09, 2.16352E-09, 1.32326E-09, 1.23223E-09, 6.05472E-09, 7.20989E-10, 2.8668E-09]

tick_labels = ['2-2', '4-4', '10-10', '20-20', '2-2', '4-4', '10-10']

### Plotting C1
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))

x = np.arange(4)
ax1.bar(x, C1_means[:4], yerr=C1_stdevs[:4], capsize=3, color='#FF0000', edgecolor = "black", linewidth=2.5) # note to self: If linewidth >=3, can see a bit of bar under the axis
ax1.set_xticks(x)
ax1.set_xticklabels(tick_labels[:4])
ax1.set_xlabel("experiment type", fontweight='bold')
ax1.set_ylabel("mean absolute error (M/L)", fontweight='bold')
ax1.set_ylim(top=1.25E-08)
ax1.set_title('Case 1', fontweight='bold')

x = np.arange(3)
ax2.bar(x, C1_means[4:], yerr=C1_stdevs[4:], capsize=3, color='#FF0000', edgecolor = "black", linewidth=2.5) # note to self: If linewidth >=3, can see a bit of bar under the axis
ax2.set_xticks(x)
ax2.set_xticklabels(tick_labels[4:])
ax2.set_xlabel("experiment type", fontweight='bold')
ax2.set_ylabel("mean absolute error (M/L)", fontweight='bold')
ax2.set_ylim(top=1.25E-08)
ax2.set_title('Case 2', fontweight='bold')

plt.suptitle("C$_1$ MAE", fontweight='bold')
plt.savefig('MAE_bar_C1.eps', format='eps', bbox_inches='tight') # saving as a vector plot
plt.savefig('MAE_bar_C1.pdf', format='pdf', bbox_inches='tight') # saving as a vector plot


### Plotting C2
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))

x = np.arange(4)
ax1.bar(x, C2_means[:4], yerr=C2_stdevs[:4], capsize=3, color='#FF0000', edgecolor = "black", linewidth=2.5) # note to self: If linewidth >=3, can see a bit of bar under the axis
ax1.set_xticks(x)
ax1.set_xticklabels(tick_labels[:4])
ax1.set_xlabel("experiment type", fontweight='bold')
ax1.set_ylabel("mean absolute error (M/L)", fontweight='bold')
ax1.set_ylim(top=2.50E-08)
ax1.set_title('Case 1', fontweight='bold')

x = np.arange(3)
ax2.bar(x, C2_means[4:], yerr=C2_stdevs[4:], capsize=3, color='#FF0000', edgecolor = "black", linewidth=2.5) # note to self: If linewidth >=3, can see a bit of bar under the axis
ax2.set_xticks(x)
ax2.set_xticklabels(tick_labels[4:])
ax2.set_xlabel("experiment type", fontweight='bold')
ax2.set_ylabel("mean absolute error (M/L)", fontweight='bold')
ax2.set_ylim(top=2.50E-08)
ax2.set_title('Case 2', fontweight='bold')

plt.suptitle("C$_2$ MAE", fontweight='bold')
plt.savefig('MAE_bar_C2.eps', format='eps', bbox_inches='tight') # saving as a vector plot
plt.savefig('MAE_bar_C2.pdf', format='pdf', bbox_inches='tight') # saving as a vector plot