from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio import PDB
import pickle
import os
import sys
import warnings
import numpy
import subprocess
# import tmhmm # works on python2 on macs
import collections
if not sys.warnoptions:
	warnings.simplefilter("ignore")

from sklearn.neural_network import MLPClassifier
import numpy as np
from sklearn.model_selection import train_test_split
# from scipy import spatial as sp

# https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
def calc_residue_dist(residue_one, residue_two) :
	"""Computes and returns the distance between two residues, by comparing the position of their alpha carbons"""

	#TODO : return an integer representing the distance between the two residues, in Angstrom
	# print(residue_one["CA"].coord - residue_two["CA"].coord)
	try:
		distance  = residue_one["CA"] - residue_two["CA"]
		# distance = numpy.sqrt(numpy.sum(diff_vector * diff_vector))
	except:
		distance = 0
	# print(distance)
	return distance

# https://github.com/llrs/PYT-SBI/blob/master/contact_map.py#L63
def compute_distance_matrix(residues) :
	"""Computes a matrix of size len(Seq) * len(Seq) with distances between each residue."""

	#TODO : return a numpy 2D array of distances between each residue of the structure.
	#Tip : you might want to make sure you filter out consecutive residues at this step.

	# size = len(residues)
	
	# answer = np.zeros((size, size), np.float)
	# for row, residue_one in enumerate(residues):
	# 	for col, residue_two in enumerate(residues):
	# 		answer[row, col] = calc_residue_dist(residue_one, residue_two)
	# 		removeConsecutives(answer, row, residue_one, col, residue_two, residues)
	# print(answer)
	# return answer

	distance_matrix = []
	for residue1 in residues:
		distance_line = []
		for residue2 in residues:
			distance = abs(residue1["CA"] - residue2["CA"])
			distance_line.append(distance)
		distance_matrix.append(distance_line)
	return distance_matrix

# https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
def extract_residues(model):
	"""Returns a list of protein residues given a PDB model"""

	#TODO : return a list of protein residues given a PDB model
	residues = []
	for chain in model:
		for residue in chain:
			if PDB.is_aa(residue, standard=True): residues.append(residue)
	# print(residues)
	return residues

def get_dssp_info(PDB_file,model,dir):
	"""Runs DSSP on protein input"""

	#TODO : you can run DSSP through biopython. The output contains a lot of useful information.
	#Tip : make sure your secondary structure indexing matches the sequence order in the PDB file!

	return PDB.DSSP(model, dir + '/' + PDB_file, dssp = 'mkdssp')

def write_fasta(sequence,PDB_file):
	"""Writes sequence to fasta file named after the PDB the sequence comes from"""

	#TODO : implement the writing of a fasta file from the sequence obtained from the PDB file.

	#return the name of the file.
	name = PDB_file.split('.')
	# fp = open("FASTAs/" + 'all.fasta', "a")
	# # # print(type(sequence))
	# fp.write(">" + name[0] + "\n" + str(sequence) + "\n")
	# fp.close()
	return str(name[0])
	# return "FASTAs/" + str(name[0]) + '.fasta'

def run_tmhmm(filename, file):
	file = file.replace(".pdb", "")
	f = open("2019-04-12_tmhmm.txt", "r")
	cur_line = 0
	tmhmm = {}
	helices = []
	for line in f:
		templ = line.strip("# ")
		if (templ[0:4] == file):
			if line[0] != "#":
				protein_name = line[0:4]
				line = line.split()
				line.remove("TMHMM2.0")
				line.remove(protein_name)
				line[1] = [int(line[1]), int(line[2])]
				d = {}
				d[line[0]] = line[1]
				helices.append(d)
			else :
				line = line.replace("#", "")
				line = line.replace(" " + file + " ", "")
				line = line.rstrip()
				line = line.split(":")
				if (line[0] != "POSSIBLE N-term signal sequence"):
					line[1] = line[1].replace(" ", "")
					tmhmm[line[0]] = line[1]
		tmhmm['name'] = file
		tmhmm['helices'] = helices
	return tmhmm

def get_contact_numbers(contact_map):
	"""Returns the proportion of residues involved in intramolecular contacts"""

	#TODO : a utility function for the number of contacts in the protein, which you might use to make sure your output makes sense
	return numpy.count_nonzero(contact_map) / float(contact_map.size)

def split_by_n(seq, n):
	'''A generator to divide a sequence into chunks of n units.'''
	while seq:
		yield seq[:n]
		seq = seq[n:]

def generate_ML_dataset(sequence,dssp_ss,tm_ss,has_contact,DSSP_vector, TMHMM_vector, oracle):
# def generate_ML_dataset(sequence,dssp_ss,has_contact,DSSP_vector, TMHMM_vector, oracle):
	# sequence amino acids
	# match 
	# list of 9 tuples
	# each tuple a string 
	# intuition of formula
	# normalise
	# print("sequence", sequence)
	# print("len_sequence", len(sequence))
	# print("dssp_ss", dssp_ss)
	# print("len_dssp_ss", len(dssp_ss))
	# print("has_contact", has_contact)
	# print("DSSP_vector", DSSP_vector)

	# print("tm_ss", tm_ss)
	helices = tm_ss['helices']
	transmem = []
	inside = []
	outside = []

	for helix in helices:
		if(helix.keys()[0] == "inside"):
			start = helix.values()[0][0]
			end = helix.values()[0][1]
			i = start
			while i<=end:
				inside.append(i-1)
				i += 1
		elif (helix.keys()[0] == "outside"):
			start = helix.values()[0][0]
			end = helix.values()[0][1]
			i = start
			while i<=end:
				outside.append(i-1)
				i += 1
		elif (helix.keys()[0] == "TMhelix"):
			start = helix.values()[0][0]
			end = helix.values()[0][1]
			i = start
			while i<=end:
				transmem.append(i-1)
				i += 1

	vecT = []
	for i in enumerate(sequence):
		# print(i[0])
		if (int(i[0]) in inside or int(i[0]) in outside):
			vecT.append('C')
		elif int(i[0]) in transmem:
			vecT.append('H')

	# print(vecT)

	tm = ''.join(vecT)
	# print(tm)

	sequence = str(sequence)
	seq_nine = list(split_by_n(sequence, 9))
	# print("seq_nine", seq_nine)
	# print(type(dssp_ss))
	# dssp_ss = SeqRecord(dssp_ss, IUPAC.protein)
	tm_nine = list(split_by_n(tm, 9))
	# print("tm_nine", tm_nine)
	dssp_nine = list(split_by_n(dssp_ss, 9))
	# print("dssp_nine", dssp_nine)

	for i in enumerate(seq_nine):
		# print("i", i[0])
		# print("seq_nine[i[0]]", seq_nine[i[0]])
		if(len(seq_nine[i[0]]) < 9): continue
		DSSP_vector.append(zip(seq_nine[i[0]], dssp_nine[i[0]]))
		TMHMM_vector.append(zip(seq_nine[i[0]], tm_nine[i[0]]))
	# print("DSSP_vector", DSSP_vector)
	# print("TMHMM_vector", TMHMM_vector)

	i = 0
	while i<len(sequence):
		if i == 0:
			oracle.append(has_contact[i+4])
		else:
			oracle.append(has_contact[i])
		i = i + 9
	# print(oracle)
	"""generates vectors for machine learning based on sequence, DSSP and TMHMM outputs"""

	#TODO : generate machine learning dataset from PDB info, and append it to the three vectors in input.

	#This function takes as input the sequence, the two secondary structures, and the has_contact boolean array.
	#use the first 4 features to append elements defined as follows to the two vectors and the oracle.

	# the features are a list of 9 tuples ("AA","SS of AA"). The oracle states whether there is intramolecular contact (<5 Ang) at the center of this subsequence of length 9.
	#An example vector element is the following : [("Y","H),("S","C"),("A","C"),("S","C"),("A","C"),("S","C"),("A","C"),("Y","H)]
	#The ith list of tuples of DSSP_vector and TMHMM_vector should describe the same sequence.

	# the oracle is a list of booleans values, establishing presence of a contact for each vector.
	# the intramolecular contact label of the ith list of tuples in DSSP_vector and TMHMM is found at index i of the oracle.

	#return the three vectors after appending the new elements.

	return DSSP_vector, TMHMM_vector, oracle

def convert_info(dssp_ss):
	newString = ""
	for s in dssp_ss:
		if(s == 'G' or s == 'I' or s == 'H'):
			newString += 'H'
		elif(s == 'X'):
			newString += ''
		else:
			newString += 'C'
	# print(newString)
	return newString

def removeConsecutives(matrix):
	contact_map = []
	row_ind = 0
	for row in matrix:
		contact_list = []
		col_ind = 0
		for column in row:
			if (column < 5 and abs(col_ind - row_ind) >= 10):
				contact_list.append(True)
			else:
				contact_list.append(False)
			col_ind += 1
		row_ind += 1
		contact_map.append(contact_list)
	contact_map = np.array(contact_map)
	return(contact_map)

def get_PDB_info(dir):
	"""Extracts sequence, DSSP secondary structure, TMHMM secondary structure and contact information from PDB files in input directory"""

	#the three vectors you are required to fill.
	DSSP_vector, TMHMM_vector, oracle = [],[],[]

	print("There are",len(os.listdir(dir)),"PDB files to parse")


	#Assemble a machine learning dataset incrementally, for each PDB file in the directory
	for ind,PDB_file in enumerate(os.listdir(dir)):
		if ind%10==0:
			print("Working on structure",ind)
		
		if(str(PDB_file) == ".DS_Store"): continue
		# if(str(PDB_file) == "2dco.pdb"): break
		#Step 1 : parse your PDB file with biopython to obtain a model object
		p = PDB.PDBParser()
		structure = p.get_structure(PDB_file[:-4].upper(), dir + "/" + PDB_file)
		model = structure[0]

		#TODO : extract a list of residues from your model object
		residues = extract_residues(model)
		print("file", PDB_file, len(residues))
		# print("residue_size",len(residues))
		# if(len(residues) > 500): continue

		#TODO : compute a distance matrix of size len(sequence)*len(sequence) with the distance between each residue
		matrix = compute_distance_matrix(residues)
		# print("here")


		#TODO : contact map should be a boolean numpy array of the same size as the distance matrix.
		#if two amino acids are within 5 angstroms of each other in 3D, but distant of at least 10 in sequence, the table should have True, else False.
		

		contact_map = removeConsecutives(matrix)
		has_contact = [True if True in contact_map[residue] else False for residue in contact_map]

		#TODO : contact info should return the proportion of residues that have an intramolecular contact in your object.
		contact_info = get_contact_numbers(contact_map)
		# print(contact_info,"contacts")

		# TODO : obtain the secondary structure prediction of the PDB model with DSSP
		dssp_info = get_dssp_info(PDB_file,model,dir)

		#TODO : obtain the sequence of the PDB file in some way of your choice.
		sequence = ""
		ppb = PDB.PPBuilder()
		for pp in ppb.build_peptides(structure):
			sequence += pp.get_sequence()

		dssp_ss = "" #ss stands for secondary structure
		dssp_seq = ""

		dssp_keys = sorted(dssp_info.keys())
		for key in dssp_keys:
			curr_ss = dssp_info[key][2]
			dssp_ss += curr_ss
			dssp_seq += dssp_info[key][1]

		converted = convert_info(dssp_ss)
		# print(dssp_ss)
		#TODO : write the sequence to a fasta file to call TMHMM with it, or to use the webserver
		filename = write_fasta(sequence,PDB_file)

		#TODO : obtain secondary structure prediction for this FASTA file with TMHMM
		# run_tmhmm will now parse tmhmmm file
		
		# test_file = "6j20"

		tm_ss = run_tmhmm(filename,PDB_file)

		# if(len(sequence) != len(residues)): continue
		DSSP_vector, TMHMM_vector, oracle = generate_ML_dataset(sequence,converted,tm_ss,has_contact,DSSP_vector, TMHMM_vector, oracle)
		# DSSP_vector, TMHMM_vector, oracle = generate_ML_dataset(sequence,converted,has_contact,DSSP_vector, TMHMM_vector, oracle)
	return DSSP_vector, TMHMM_vector, oracle



def generate_dataset():
	"""Runs the PDB parsing utility"""
	# DSSP_vector, TMHMM_vector, oracle = get_PDB_info("PDBs/")
	DSSP_vector, TMHMM_vector, oracle = get_PDB_info("/Users/jakezhu/Documents/COMP564/W2019/PDBs/")

	#store a pickle of your results to avoid repeating get_PDB_info
	pickle.dump((DSSP_vector, TMHMM_vector, oracle),open("no_split_dataset.pickle","wb"))
	# pickle.dump((DSSP_vector, TMHMM_vector, oracle),open("ML_ready_dataset.pickle","wb"))
	return DSSP_vector, TMHMM_vector, oracle


def split_dataset(X, Y):
	"""Splits the dataset into training and testing"""
	# TODO : split the X and Y dataset into a reasonable training set and a test set. Your test set should have have 20% of the datapoints.

	# Tip : look up train_test_split with scikit-learn

	X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.20)
	new_Y_train = []
	new_X_train = []
	count_x = 0
	count_y = 0

	for(y_ind, y_val) in enumerate(Y_train):
		if(y_val == True):
			new_Y_train.append(y_val)
			new_X_train.append(X_train[y_ind])
			count_y += 1
	for(y_ind, y_val) in enumerate(Y_train):
		if(count_x < count_y and y_val == False):
			new_Y_train.append(y_val)
			new_X_train.append(X_train[y_ind])
			count_x += 1
	return X_train, X_test, Y_train, Y_test
	# return new_X_train, X_test, new_Y_train, Y_test


def format_simple_dataset(vector,solutions):
	"""takes as input a vector of sequence strings and a vector of booleans
	outputs a vector of size 9 vectors of tuples and an array of binary numbers"""

	AA = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T", "X"]
	contact_truth = [False, True]

	training_size = len(vector)
	formatted_X = []
	for i in vector[:training_size]:
		current_vec = np.zeros(189)
		for ind, j in enumerate(i):
			nuc = AA.index(j[0])
			current_vec[21 * ind + nuc] = 1
		current_vec = np.array(current_vec)
		formatted_X.append(current_vec)
	X = np.array(formatted_X)
	formatted_Y = np.zeros(training_size)
	for i in range(training_size):
		formatted_Y[i] = contact_truth.index(solutions[i])
	Y = np.array(formatted_Y)
	return X,Y

def format_one_hot_dataset(vector,solutions):
	"""Takes as input a vector of 9 sequence, ss tuples, outputs a vector of one-hot vectors of size 198 and a binary vector"""
	AA = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T", "X"]
	ss = ["H", "C"]
	contact_truth = [False, True]

	training_size = len(vector)
	formatted_X = []
	for i in vector[:training_size]:
		current_vec = np.zeros(198)
		for ind, j in enumerate(i):
			nuc = AA.index(j[0])
			sec = ss.index(j[1])
			current_vec[22 * ind + nuc] = 1
			current_vec[22 * ind + 21] = sec
		current_vec = np.array(current_vec)
		formatted_X.append(current_vec)
	X = np.array(formatted_X)
	formatted_Y = np.zeros(training_size)
	for i in range(training_size):
		formatted_Y[i] = contact_truth.index(solutions[i])
	Y = np.array(formatted_Y)
	return X,Y


def run_NN_on_sequence(vector, solutions):
	"""Trains a scikit-learn basic 4 layers neural network on a sequence vector, outputs results"""

	#X is a the input vector of ML-ready numpy arrays, Y is the corresponding oracle
	X, Y = format_simple_dataset(vector, solutions)

	#TODO : fill in split_dataset
	X_train, X_test, Y_train, Y_test = split_dataset(X,Y)


	clf = MLPClassifier(solver='sgd', alpha=1e-6, hidden_layer_sizes=(156), activation="logistic", shuffle=True,
						verbose=False, random_state=1, tol=1e-5, max_iter=350)
	clf.fit(X_train, Y_train)
	print("score, sequence only")
	print("TRAINING ACCURACY", clf.score(X_train, Y_train))
	print("TEST ACCURACY", clf.score(X_test, Y_test))

def run_NN_for_ss_type(vector,solutions):
	"""Trains a scikit-learn basic 4 layers neural network on a sequence-structure vector, outputs results"""


	#X is a the input vector of ML-ready numpy arrays, Y is the corresponding oracle
	X,Y = format_one_hot_dataset(vector,solutions)

	#TODO : fill in split_dataset
	X_train, X_test, Y_train, Y_test = split_dataset(X,Y)



	clf = MLPClassifier(solver='sgd', alpha=1e-4, hidden_layer_sizes=(156), activation="logistic", tol=1e-5,
						shuffle=True, verbose=False, random_state=1, max_iter=350)

	clf.fit(X_train,Y_train)
	print("score, ss and sequence")
	print("TRAINING ACCURACY",clf.score(X_train, Y_train))
	print("TEST ACCURACY",clf.score(X_test, Y_test))



def predict_intramolecular_contacts(dataset):
	"""Compares neural network results for DSSP and TMHMM secondary structures"""
	DSSP_vector, TM_vector, solutions = dataset
	print("SEQUENCE ONLY RESULTS")

	#use either vector for sequence only prediction
	run_NN_on_sequence(DSSP_vector,solutions)
	print("DSSP RESULTS")
	run_NN_for_ss_type(DSSP_vector,solutions)
	print("TM RESULTS")
	run_NN_for_ss_type(TM_vector,solutions)

	#TODO : analyze your results!

if __name__== "__main__":
	"""executes program. If you already have a dataset, you can load it via pickle"""

	# TODO : follow the instructions in get_PDB_info()
	dataset = generate_dataset()

	print("new_X,Y_train", 0)
	with open("no_split_dataset.pickle", "rb") as f:
		dataset0 = pickle.load(f)
		predict_intramolecular_contacts(dataset0)

	print("new_X,Y_train", 1)
	with open("ML_ready_dataset.pickle", "rb") as f:
		dataset1 = pickle.load(f)
		predict_intramolecular_contacts(dataset1)
	
