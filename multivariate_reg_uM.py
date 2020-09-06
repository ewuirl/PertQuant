import matplotlib.pyplot as plt
import pandas as pd
import math

def read_eq_data_file(file_name):
	'''This function reads data from a file (named file_name) generated by 
	gen_eq_data and creates Ci_all_array and Am_array from it. It returns these
	arrays as (Ci_all_array, Am_array).'''
	with open(file_name, 'r') as file:
		# Read all the lines in the file
		lines = file.readlines()

		# Iterate through the first 6 lines in the file to figure out L, N, 
		# N_runs, Cmin, Cmax, and Ai
		for i in range(6):
			# Remove the newline character at the end of the line
			line = lines[i].rstrip('\n')
			# Remove the whitespace
			line_list = line.rsplit()
			# Figure out how many C strands there are
			if i == 0:
				L = int(line_list[2])
			# Figure out how many A strands there are
			elif i == 1:
				N = int(line_list[2])
			# Figure out how many runs there are
			elif i == 2:
				N_runs = int(line_list[2])
			elif i == 3:
				Cmin = float(line_list[2])
			elif i == 4:
				Cmax = float(line_list[2])
			else:
				Ai = float(line_list[2])

		# Make arrays to store the data in 
		Ci_all_array = np.zeros((N_runs, L))
		Am_array = np.zeros((N_runs, N))

		# Iterate through the lines to add each row to Ci_all_array and Am_array
		for j in range(len(lines)-6):
			i = j+6
			# Remove the newline character at the end of the line
			line = lines[i].rstrip('\n')
			# Remove the whitespace
			line_list = line.rsplit()
			# Write the Ci values
			for l in range(L):
				Ci = float(line_list[l])
				Ci_all_array[j][l] = Ci
			# Write the Am values
			for n in range(N):
				Am = float(line_list[L+n])
				Am_array[j][n] = Am

	return(Ci_all_array, Am_array, Cmin, Cmax, Ai)

def convert_np2df(array, strand_type):
	'''This function takes in a numpy array (either Ci_all_array or Am_array) 
	and strand_type, a string that specifies strand type (eg 'C', 'A'). 
	It converts the numpy array into a pandas data frame indexed by run and 
	labeled with strand numbers.'''
	# Figure out the dimensions of the array.
	N_runs, N_strands = np.shape(array)
	# Create the indices
	index = np.arange(N_runs)
	# Create a list to store the labels
	label_list = []
	# Create the labels
	for i in range(N_strands):
		label = strand_type + str(i)
		label_list.append(label)
	# Create the data frame
	df = pd.DataFrame(array, index=index, columns=label_list)

	return df

def get_stats(strand_df):
	'''This function takes a dataframe of strand concentrations and returns
	a transposed description of the statistics of the dataframe.'''
	strand_stats = strand_df.describe()
	strand_stats = strand_stats.transpose()
	return strand_stats

def Z_normalize_data(strand_df, strand_stats):
	'''This data takes a dataframe of strand concentrations and its statistics 
	and returns a dataframe of Z normalized data.'''
	return (strand_df - strand_stats['mean']) / strand_stats['std']

def prep_data(Am_array, Ci_all_array, train):
	'''This function takes in Am_array and Ci_all_array and train, a decimal 
	number between [0,1]. Am_array and Ci_all_array are arrays of data read in 
	using the function read_eq_data_file from gen_eq_data.py. train is the 
	fraction of the data that will be set aside as the training data. This 
	function returns a tuple of the following:
		Am_df = Am_array converted to a pandas dataframe,
		Ci_df = Ci_all_array converted to a pandas dataframe,
		train_Am_df = the subset of Am_df to be used for training,
		train_Ci_df = the subset of Ci_df to be used for training,
		test_Am_df = the subset of Am_df to be used for testing,
		test_Ci_df = the subset of Ci_df to be used for testing.'''
	# Convert array data to data frames
	Am_df = convert_np2df(Am_array, 'A')
	Ci_df = convert_np2df(Ci_all_array, 'C')

	# Pick out test and training Ci data
	train_Am_df_rand = Am_df.sample(frac=train,random_state=0)
	test_Am_df = Am_df.drop(train_Am_df_rand.index)
	train_Am_df = Am_df.drop(test_Am_df.index)

	# Separate test and training Am data
	test_Ci_df = Ci_df.drop(train_Am_df.index)
	train_Ci_df = Ci_df.drop(test_Am_df.index)

	return(Am_df, Ci_df, train_Am_df, train_Ci_df, test_Am_df, test_Ci_df)

def plot_data(Am_df, Ci_df, Cmax, **kwargs):
	'''This function can be used to plot the data in Am_df and Ci_df pairwise in
	an L x N grid. The keyword arguments are title, green_list, and y. title is 
	a string to use as the plot title. The default is 'Am vs Ci Data'. 
	green_list is a list of tuples that represent grid positions to plot green 
	instead of the standard blue. y is a decimal number specifying the y 
	location of the title in figure coordinates. The default is 0.92'''
	# Figure out how many strands we need to plot
	N = len(Am_df.keys())
	L = len(Ci_df.keys())
    
	if 'green_list' in kwargs:
		green_list = kwargs.get('green_list')
	else:
		green_list = []

	# Figure out the dimensions of the plot
	width = N*4.5
	height = L*3.25

	# Start plotting
	fig, ax = plt.subplots(L,N,figsize=(width, height))
	for l in range(L):
		Ci = Ci_df.keys()[l]
		# Label the y axes
		ax[l,0].set_ylabel(Ci)
		# Plot the data
		for n in range(N):
			Ai = Am_df.keys()[n]
			ax[l,n].set_ylim(0,Cmax)
			ax[l,n].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
			if l == n:
				ax[l,n].scatter(Am_df[Ai], Ci_df[Ci], c='g', marker='.')
			elif (l,n) in green_list:
				ax[l,n].scatter(Am_df[Ai], Ci_df[Ci], c='g', marker='.')
			else:
				ax[l,n].scatter(Am_df[Ai], Ci_df[Ci], marker='.')

			# Label the x axes
			if l == L-1:
				ax[l,n].set_xlabel(Ai)

	# Add the title and other labels
	if 'title' in kwargs:
		title = kwargs.get('title')
	else:
		title = 'Am vs Ci Data'

	if 'y' in kwargs:
		y = kwargs.get('y')
	else:
		y = 0.92
	fig.suptitle(title, y=y)


def plot_predict(predictions, true_vals, Cmax, **kwargs):
	'''This function takes in an array of predicted Ci values and an array of 
	true Ci values and plots them against each other in a grid. Each plot also 
	has a green y = x line. The function takes two optional keyword arguments, 
	'lims' and 'title'. 'lims' should be a list of end values to draw the y = x 
	line between. The default is [-0.25,2.25]. 'title' is a string to use as a 
	plot title. The default is 'Ci Predictions vs True Values'.'''
	L = len(true_vals.keys())

	# Figure size of plot
	if L <= 4:
		cols = L
	else:
		cols = 4
	rows = math.ceil(L/4)
	width = cols*4.5
	height = rows*4

	# Determine limits of y = x plot.
	if 'lims' in kwargs:
		lims = kwargs.get('lims')
	else:
		lims = [0,2e-6]

	# Plot!
	fig, ax = plt.subplots(rows,cols,figsize=(width, height))

	# Add the title
	if 'title' in kwargs:
		title = kwargs.get('title')
	else:
		title = 'Ci Predictions vs True Values'
	fig.suptitle(title)

	# Handle 1 row case
	if rows < 2:
		ax[0].set_ylabel('Predictions [Ci]')
		for l in range(L):
			Ci = true_vals.keys()[l]
			ax[l].scatter(true_vals[Ci], predictions[:,l], marker='.')
			ax[l].plot(lims, lims, 'g')
			ax[l].set_xlabel('True Values [{}]'.format(Ci))
			ax[l].set_ylim(0,Cmax)
			ax[l].set_xlim(0,Cmax)
			ax[l].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
			ax[l].ticklabel_format(axis="x", style="sci", scilimits=(0,0))

	# Handle multiple row case
	elif rows >= 2:
		for row in range(rows):
			ax[row].set_ylabel('Predictions [Ci]')

		for l in range(L):
			Ci = true_vals.keys()[l]
			row = math.floor(l/4)
			col = l % 4 - 1
			if col == -1:
				col = 3
			ax[row,col].scatter(true_vals[Ci], predictions[:,l], marker='.')
			ax[row,col].plot(lims, lims, 'g')
			ax[row,col].set_xlabel('True Values [{}]'.format(Ci))
			ax[row,col].set_ylim(0,Cmax)
			ax[row,col].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
			ax[row,col].set_xlim(0,Cmax)
			ax[row,col].ticklabel_format(axis="x", style="sci", scilimits=(0,0))

def plot_error_hist(error, **kwargs):
	'''This function plots a grid of histograms of the errors of the Ci 
	predictions. The keyword arguments are 'bins' and 'title'. 'bins' should be
	an integer to specify the number of bins the histograms have. The default is
	25 bins. 'title' is a string for the title. The default is 'Ci Prediction
	Error'.'''

	L = len(error.keys())

	# Figure size of plot
	if L <= 4:
		cols = L
	else:
		cols = 4
	rows = math.ceil(L/4)
	width = cols*4.5
	height = rows*4

	# Determine number of bins.
	if 'bins' in kwargs:
		bins = kwargs.get('bins')
	else:
		bins = 25
		
	# Plot!
	fig, ax = plt.subplots(rows,cols,figsize=(width, height))
	# Figure out the title
	if 'title' in kwargs:
		title = kwargs.get('title')
	else:
		title = 'Ci Prediction Error'
	fig.suptitle(title)

	# Handle 1 row case
	if rows < 2:
		ax[0].set_ylabel('Count')
		for l in range(L):
			Ci = error.keys()[l]
			ax[l].hist(error[Ci], bins = bins)
			ax[l].set_xlabel('Prediction Error [{}]'.format(Ci))

	# Handle multiple row case
	elif rows >= 2:
		for row in range(rows):
			ax[row].set_ylabel('Count')

		for l in range(L):
			Ci = true_vals.keys()[l]
			row = math.floor(l/4)
			col = l % 4 - 1
			if col == -1:
				col = 3
			ax[row,col].hist(error[Ci], bins = bins)
			ax[row,col].set_xlabel('True Values [{}]'.format(Ci))


def plot_true_v_error(true_vals, error, **kwargs):
	'''This function plots a grid of scatter plots of the true values vs the 
	errors of the Ci predictions. The keyword argument is 'title', which should
	be a string for the title. The default is 'True Values vs Error'.'''
	L = len(error.keys())

	# Figure size of plot
	if L <= 4:
		cols = L
	else:
		cols = 4
	rows = math.ceil(L/4)
	width = cols*4.5
	height = rows*4

	# Plot!
	fig, ax = plt.subplots(rows,cols,figsize=(width, height))
	# Figure out the title
	if 'title' in kwargs:
		title = kwargs.get('title')
	else:
		title = 'True Values vs Error'
	fig.suptitle(title)
		

	# Handle 1 row case
	if rows < 2:
		ax[0].set_ylabel('Prediction Error [Ci]')
		for l in range(L):
			Ci = error.keys()[l]
			ax[l].scatter(true_vals[Ci], error[Ci], marker='.')
			ax[l].set_xlabel('True Values [{}]'.format(Ci))
			ax[l].set_xlim(0,Cmax)
			ax[l].ticklabel_format(axis="x", style="sci", scilimits=(0,0))
			ax[l].ticklabel_format(axis="y", style="sci", scilimits=(0,0))

	# Handle multiple row case
	elif rows >= 2:
		for row in range(rows):
			ax[row].set_ylabel('Count')

		for l in range(L):
			Ci = true_vals.keys()[l]
			row = math.floor(l/4)
			col = l % 4 - 1
			if col == -1:
				col = 3
			ax[row,col].scatter(true_vals[Ci], error[Ci], marker='.')
			ax[row,col].set_xlabel('True Values [{}]'.format(Ci))
			ax[row,col].set_xlim(0,Cmax)
			ax[row,col].ticklabel_format(axis="x", style="sci", scilimits=(0,0))
			ax[row,col].ticklabel_format(axis="y", style="sci", scilimits=(0,0))


def plot_residuals(Am_df, error, **kwargs):
	'''This function is used to plot Am_df and oredicted Ci error pairwise in
	an L x N grid. he keyword arguments are title, green_list, and y. title is 
	a string to use as the plot title. The default is 'Residuals'. green_list is
	a list of tuples that represent grid positions to plot green instead of the 
	standard blue. y is a decimal number specifying the y location of the title 
	in figure coordinates. The default is 0.92'''
	# Figure out how many strands we need to plot
	N = len(Am_df.keys())
	L = len(error.keys())

	# Figure out the dimensions of the plot
	width = N*4.5
	height = L*3.25

	if 'green_list' in kwargs:
		green_list = kwargs.get('green_list')
	else:
		green_list = []

	# Start plotting
	fig, ax = plt.subplots(L,N,figsize=(width, height))
	for l in range(L):
		Ci = error.keys()[l]
		# Label the y axes
		ax[l,0].set_ylabel('Prediction Error [{}]'.format(Ci))
		# Plot the data
		for n in range(N):
			Ai = Am_df.keys()[n]
			ax[l,n].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
			if (l,n) in green_list:
				ax[l,n].scatter(Am_df[Ai], error[Ci], c='g', marker='.')
			else:
				ax[l,n].scatter(Am_df[Ai], error[Ci], marker='.')

			# Label the x axes
			if l == L-1:
				ax[l,n].set_xlabel(Ai)

	# Add the title and other labels
	if 'title' in kwargs:
		title = kwargs.get('title')
	else:
		title = 'Residuals'

	if 'y' in kwargs:
		y = kwargs.get('y')
	else:
		y = 0.92
	fig.suptitle(title, y=y)

def subset(Am_df, Ci_df, frac, train):
	'''This function takes in two data frames, Am_df (input) and Ci_df (output),
	and 2 real numbers between 0 and 1, frac (subset size) and train 
	(training subset size). It uses frac to pick out subsets of Am_df and Ci_df,
	and then divides those subsets into training sets of size train, and a
	testing sets of size 1-train. It also generates a large testing set, which 
	combines the testing set with the 1-frac subset. The function returns these
	data frames in a tuple: (train_Am_df, train_Ci_df, test_Am_df, test_Ci_df, 
	large_test_Am_df, large_test_Ci_df)'''
	# Pick out the subset of data to use.
	Am_df_subset_drop = Am_df.sample(frac=1-frac,random_state=0)
	Am_df_subset = Am_df.drop(Am_df_subset_drop.index)

	# Pick out test and training Am data
	train_Am_df_subset_rand = Am_df_subset.sample(frac=train,random_state=0)
	test_Am_df = Am_df_subset.drop(train_Am_df_subset_rand.index)
	train_Am_df = Am_df_subset.drop(test_Am_df.index)
	large_test_Am_df = Am_df.drop(train_Am_df.index)

	# Separate out the subset of Ci data to use.
	Ci_df_subset = Ci_df.drop(Am_df_subset_drop.index)
	# Pick out test and training Ci data
	test_Ci_df = Ci_df_subset.drop(train_Am_df.index)
	train_Ci_df = Ci_df_subset.drop(test_Am_df.index)
	large_test_Ci_df = Ci_df.drop(train_Am_df.index)
	return(train_Am_df, train_Ci_df, test_Am_df, test_Ci_df, large_test_Am_df, \
		large_test_Ci_df)



