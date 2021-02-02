
"""
"""
def read_line_FASTQ(file, counter):
	seq_ID = file.readline()
	seq = file.readline()
	file.readline()
	Q_score = file.readline()
