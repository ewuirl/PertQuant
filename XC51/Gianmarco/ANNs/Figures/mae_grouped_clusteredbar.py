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

# Squares in legend rather than rectangles
plt.rcParams['legend.handlelength'] = 1
plt.rcParams['legend.handleheight'] = 1.125
  
# create data. Using test performance
# Order is: 2-2-a, 4-4-a, 10-10-a, 20-20-a, 2-2-sa, 4-4-sa, 10-10-sa 
C1_means = [8.79E-09, 5.92E-09, 2.84E-09, 1.43E-09, 8.07E-09, 6.19E-09, 6.81E-09]
C2_means = [1.81E-08, 5.80E-09, 2.65E-09, 2.65E-09, 8.17E-09, 5.09E-09, 6.07E-09]
# For grouped, each tuple is (high, low)
C1_grouped = [(1.23E-07, 2.68E-08), (8.53E-08, 6.77E-08), (6.79E-10, 1.58E-08), 
(9.46E-10, 2.51E-08), (5.65E-08, 6.42E-08), (5.66E-08, 1.29E-07), (2.27E-08, 1.42E-08)]
C2_grouped = [(1.87E-07, 7.67E-08), (5.39E-08, 9.25E-09), (1.82E-08, 1.30E-08), 
(8.36E-09, 3.94E-09), (8.65E-08, 5.91E-08), (8.46E-08, 6.03E-08), (2.49E-08, 4.29E-08)]

C1_low = [item[1] for item in C1_grouped]
C1_high = [item[0] for item in C1_grouped]
C2_low = [item[1] for item in C2_grouped]
C2_high = [item[0] for item in C2_grouped]

tick_labels = ['2-2', '4-4', '10-10', '20-20', '2-2', '4-4', '10-10']

width = 0.20

plt.figure(figsize=(15, 8)) 

x = np.arange(7)

### Plotting C1
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))

x = np.arange(4)
C1_low_plot = ax1.bar(x-0.2, C1_low[:4], width, label='grouped low', color='#FF0000', edgecolor = "black", linewidth=2.5) # note to self: If linewidth >=3, can see a bit of bar under the axis
C1_plot = ax1.bar(x, C1_means[:4], width, label='random', color='#0000FF', edgecolor = "black", linewidth=2.5)
C1_high_plot = ax1.bar(x+0.2, C1_high[:4], width, label='grouped high', color='#FFFF00', edgecolor = "black", linewidth=2.5)
ax1.set_xticks(x)
ax1.set_xticklabels(tick_labels[:4])
ax1.set_xlabel("experiment type", fontweight='bold')
ax1.set_ylabel("mean absolute error (M/L)", fontweight='bold')
# ax1.set_ylim(top=1.25E-08)
ax1.set_title('Case 1', fontweight='bold')

legend = ax1.legend(handles=[C1_low_plot, C1_plot, C1_high_plot], loc='upper right', edgecolor='black') # 0.01 shifts the legend right slightly so it doesn't touch the y-axis tick marks
frame = legend.get_frame()
frame.set_linewidth(2.5) # setting the line width of the legend


x = np.arange(3)
C1_low_plot = ax2.bar(x-0.2, C1_low[4:], width, label='grouped low', color='#FF0000', edgecolor = "black", linewidth=2.5) # note to self: If linewidth >=3, can see a bit of bar under the axis
C1_plot = ax2.bar(x, C1_means[4:], width, label='random', color='#0000FF', edgecolor = "black", linewidth=2.5)
C1_high_plot = ax2.bar(x+0.2, C1_high[4:], width, label='grouped high', color='#FFFF00', edgecolor = "black", linewidth=2.5)
ax2.set_xticks(x)
ax2.set_xticklabels(tick_labels[4:])
ax2.set_xlabel("experiment type", fontweight='bold')
ax2.set_ylabel("mean absolute error (M/L)", fontweight='bold')
# ax1.set_ylim(top=1.25E-08)
ax2.set_title('Case 2', fontweight='bold')

legend = ax2.legend(handles=[C1_low_plot, C1_plot, C1_high_plot], loc='upper right', edgecolor='black') # 0.01 shifts the legend right slightly so it doesn't touch the y-axis tick marks
frame = legend.get_frame()
frame.set_linewidth(2.5) # setting the line width of the legend


plt.suptitle("C$_1$ MAE", fontweight='bold')
plt.savefig('MAE_grouped_C1.eps', format='eps', bbox_inches='tight') # saving as a vector plot
plt.savefig('MAE_grouped_C1.pdf', format='pdf', bbox_inches='tight') # saving as a vector plot


### Plotting C2
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))

x = np.arange(4)
C2_low_plot = ax1.bar(x-0.2, C2_low[:4], width, label='grouped low', color='#FF0000', edgecolor = "black", linewidth=2.5) # note to self: If linewidth >=3, can see a bit of bar under the axis
C2_plot = ax1.bar(x, C2_means[:4], width, label='random', color='#0000FF', edgecolor = "black", linewidth=2.5)
C2_high_plot = ax1.bar(x+0.2, C2_high[:4], width, label='grouped high', color='#FFFF00', edgecolor = "black", linewidth=2.5)
ax1.set_xticks(x)
ax1.set_xticklabels(tick_labels[:4])
ax1.set_xlabel("experiment type", fontweight='bold')
ax1.set_ylabel("mean absolute error (M/L)", fontweight='bold')
# ax1.set_ylim(top=1.25E-08)
ax1.set_title('Case 1', fontweight='bold')

legend = ax1.legend(handles=[C2_low_plot, C2_plot, C2_high_plot], loc='upper right', edgecolor='black') # 0.01 shifts the legend right slightly so it doesn't touch the y-axis tick marks
frame = legend.get_frame()
frame.set_linewidth(2.5) # setting the line width of the legend


x = np.arange(3)
C2_low_plot = ax2.bar(x-0.2, C2_low[4:], width, label='grouped low', color='#FF0000', edgecolor = "black", linewidth=2.5) # note to self: If linewidth >=3, can see a bit of bar under the axis
C2_plot = ax2.bar(x, C2_means[4:], width, label='random', color='#0000FF', edgecolor = "black", linewidth=2.5)
C2_high_plot = ax2.bar(x+0.2, C2_high[4:], width, label='grouped high', color='#FFFF00', edgecolor = "black", linewidth=2.5)
ax2.set_xticks(x)
ax2.set_xticklabels(tick_labels[4:])
ax2.set_xlabel("experiment type", fontweight='bold')
ax2.set_ylabel("mean absolute error (M/L)", fontweight='bold')
# ax1.set_ylim(top=1.25E-08)
ax2.set_title('Case 2', fontweight='bold')

legend = ax2.legend(handles=[C2_low_plot, C2_plot, C2_high_plot], loc='upper right', edgecolor='black') # 0.01 shifts the legend right slightly so it doesn't touch the y-axis tick marks
frame = legend.get_frame()
frame.set_linewidth(2.5) # setting the line width of the legend


plt.suptitle("C$_2$ MAE", fontweight='bold')
plt.savefig('MAE_grouped_C2.eps', format='eps', bbox_inches='tight') # saving as a vector plot
plt.savefig('MAE_grouped_C2.pdf', format='pdf', bbox_inches='tight') # saving as a vector plot
