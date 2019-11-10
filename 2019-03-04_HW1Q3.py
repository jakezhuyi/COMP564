import sys, re
import math
import collections
# import numpy as np
# ============================================ Student info methods================================================
def get_student_name():
	# @TO_STUDENT: Write your name here
	student_name = "Jake Zhu"
	if not student_name:
		raise ValueError("Error: you forgot to add your name in get_student_name method.")
	return student_name

def get_student_id():
	# @TO_STUDENT: Write your student id here
	student_id = "260621270"
	if not student_id:
		raise ValueError("Error: you forgot to add your student id in get_student_id method.")
	return student_id
# =================================================================================================================

# =================================== Validate input/output methods================================================
def validate_Q3_1_input_format(subopt_result):
	if not isinstance(subopt_result, list) or [sr for sr in subopt_result if re.search("[^\(\)\.]", sr)]:
		raise ValueError("Error: your input should be list of strings (secondary structures, ie alphabet is '().')")

def validate_Q3_1_output_format(result):
	if not isinstance(result, list) or [sr for sr in result if not isinstance(sr, list)]:
		raise ValueError("Error: your output should be [ [i1, j1, freq_i1_j1 ], [i2, j2, freq_i2_j2 ], ...  ].")
	if [sr for sr in result if sr[0] >= sr[1]]:
		raise ValueError("Error: check your i and j values. They should satisfy: i > j.")
# =================================================================================================================
# ================================================== Helper methods================================================
def parse_subopt_result_file(filepath):
	'''
	Parsing of a standard txt file that contains result of
		"RNAsubopt -p __k__ < myFasta.fasta > subopt_result_file.txt"
		where __k__ is parameter. (Filename chosen randomly. Please, feel free to use your own names.)
	@args:
	filepath: (full or relative) path to the subopt_result_file.txt.
	(Filename chosen randomly. Please, feel free to use your own names.)

	@return: subopt_result: list of the strings (assumed to be secondary structures)
	'''
	subopt_result = []
	with open(filepath, 'r') as f:
		for i, line in enumerate(f):
			if i < 2:
				continue
			subopt_result.append(line.strip())
	return subopt_result

def parse_dot_ps_file(filepath):
	'''
	Parsing of a dot.ps file that contains result of RNAfold program
	@args:
	filepath: (full or relative) path to the dot.ps.
	@return:
	dot_ps_result: list f lists with i, j, freq_i_j
	'''
	dot_ps_result = []
	with open(filepath, 'r') as f:
		is_data = False
		for line in f:
			if not is_data and line.startswith('%start of base pair probability data'):
				is_data = True
				continue
			elif is_data and line.startswith('showpage'):
				break
			elif is_data:
				if line.find('ubox') > 0:
					# take only first 3 numbers
					data_line = line.split()[:3]
					dot_ps_result.append(
						[int(data_line[0]), int(data_line[1]), float(data_line[2])]
					)
	return dot_ps_result

# =================================================================================================================
def get_answer_Q3_1(subopt_result):
	'''
	This method should be implemented by student.
	@args:
	subopt_result: a list of the secondary structures produced by RNAsubopt -p <k> for particular input

	@return:
	result: list of lists (example is below) with indices and relevant frequency.
	example [ [0, 1, 0.10], [0, 2, 0.15], [0, 3, 0.16], ....... ]

	@note:
	check input/output as advised in code. Question will be marked as 0 in case of not following the formats.
	'''
	# basic check for the proper input
	validate_Q3_1_input_format(subopt_result)
	# print(subopt_result)
	# @TO_STUDENT: Write your code here
	# @TO_STUDENT: use result variable for results. below is an example of an expected format for result.

	xStack = []
	resultingPairs = []
	result = []

	# https://stackoverflow.com/questions/41066448/finding-out-rna-base-pairing-in-given-structure
	for seq in subopt_result:
		for i, x in enumerate(seq):
			if x == '(':
				xStack.append(i)
			elif x == ')':
				o = xStack.pop()
				tempPair = [o, i]
				resultingPairs.append(tempPair)
			else:
				continue
	
	tuples = [tuple(t) for t in resultingPairs]
	counter = collections.Counter(tuples)
	setTuples = set(tuples)
	# print resultingPairs
	for c in counter:
		a_list = list(c)
		a_list.append(float(counter[c])/len(subopt_result))
		# a_list.append(np.longdouble(counter[c])/len(subopt_result))
		result.append(a_list)

	# print result

	# result = [ [0, 1, 0.10], [0, 2, 0.15], [0, 3, 0.16] ]

	# @TO_STUDENT: output should be [ [i1, j1, freq_i1_j1 ], [i2, j2, freq_i2_j2 ], ...  ]
	# use validate_Q3_output_format(result) to validate the output
	validate_Q3_1_output_format(result)
	return result

def get_answer_Q3_2(q3_1_result, dot_ps_result):
	'''
	This method should be implemented by student.
	Compare output from RNAfold and result of question3_1 for the same sequence and return an error (see text assignment)
	result_error is expected to be numeric
	'''

	# print (q3_1_result)
	# print (dot_ps_resuclt)
	result_error = 0.0
	error_squared = 0.0
	# @TO_STUDENT: Write your code here (trust me, answer is not 0 :-) )

	for dot in dot_ps_result:
		for q3 in q3_1_result:
			# Offset due to Array Location starting at 0 versus starting at 1
			# For some reason I had an old version of the code and it was taking lbox data
			# Spent a week debugging 0.5-0.6 errors
			if(dot[0] == (q3[0]+1) and dot[1] == (q3[1]+1)):
			# if(dot[0] == (q3[0]+1) and dot[1] == (q3[1]+1) and dot[2] != 0.95)
				diff = (dot[2] ** 2 - q3[2]) ** 2
				# diff = (np.longdouble(dot[2]) ** 2 - q3[2]) ** 2
				error_squared += diff
				# break

	result_error = math.sqrt(error_squared)
	# print result_error
	return result_error

# @TO_STUDENT: You can test your methods below by calling methods. Workflow is given already (can be changed).
# @TO_STUDENT: Everything below this point will not be considered as a solution and will be deleted for grading.
# @TO_STUDENT: Advise: be careful with idents. Use only tabs, or only FOUR spaces. NEVER mix them.

print("This is a solution of %s, student_id is %s" % (get_student_name(), get_student_id()) )

for i in (10, 50, 100, 1000, 10000):
	i = str(i) 
	subopt_result_filepath = "subopt_" + i + ".txt"
	dot_ps_filepath = "HW1Q3/dot.ps"

	# parsing RNAsubopt result file
	subopt_result = parse_subopt_result_file(subopt_result_filepath)

	# solving quesion Q3_1
	q3_1_result = get_answer_Q3_1(subopt_result)

	# parsing dot.ps file
	dot_ps_result = parse_dot_ps_file(dot_ps_filepath)
	# print(dot_ps_result)

	# solving question Q3_2
	q3_2_result = get_answer_Q3_2(q3_1_result, dot_ps_result)
	print("%.7f" % q3_2_result, end=",")

print("")
