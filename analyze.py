#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division

from static import PATH
from boolfunc import PolynomialOverF2, BoolFunction, timer
from multiprocessing  import Process, Manager

import codecs

'''ANALYZE BOOL FUNCTION'''
def exception_checker(function):

	def wrapper(bool_function, *args):

		if not isinstance(bool_function, BoolFunction):
			raise TypeError("bool function must be an isinstance of BoolFunction class")

		if len(bool_function.truth_table) == 0:
			raise Exception("Truth table is empty, fill it from file or call bool_function.create_truth_table method")

		if len(bool_function.anf_table) == 0:
			raise Exception("ANF table is empty, fill it from file or call bool_function.create_anf_table method")

		if len(bool_function.walsh_spectrum_table[0]) == 0:
			raise Exception("Walsh spectrum table is empty, fill it from file or call function.create_anf_table method")

		return function(bool_function, *args)

	return wrapper

@exception_checker
@timer
def analyze_algebraic_degree(bool_function):

	def filter_false_value(enum_object):

		return bool(enum_object[1]) 

	anf_table = bool_function.anf_table
	n = bool_function.field_extension
	algebraic_degree = [ 0 for f in range(n) ]

	anf = filter(filter_false_value, enumerate(anf_table))

	for coef_index, value in anf:
		deg = sum([ bool(coef_index & (1 << i)) for i in range(n) ])
		for f in range(n):
			if value & (1 << f):
				algebraic_degree[f] = max(algebraic_degree[f], deg)

	return algebraic_degree

@exception_checker
@timer
def analyze_disbalance(bool_function):

	walsh_spectrum_table = bool_function.walsh_spectrum_table
	n = bool_function.field_extension
	disbalance = [ 0 for f in range(n) ]

	for f in range(n):
		disbalance[f] = walsh_spectrum_table[f][0]

	return disbalance

@exception_checker
@timer
def analyze_nonlinearity(bool_function):
	
	walsh_spectrum_table = bool_function.walsh_spectrum_table
	n = bool_function.field_extension
	nonlinearity = [ 0 for f in range(n) ]

	for f in range(n):
		nonlinearity[f] = pow(2, n - 1) - pow(2, -1) * max(map(abs, iter(walsh_spectrum_table[f])))

	return nonlinearity

@exception_checker
@timer
def analyze_correlation_immunity(bool_function):

	walsh_spectrum_table = bool_function.walsh_spectrum_table
	n = bool_function.field_extension
	k_balance = [ 0 for f in range(n) ]
	weight_table = {}

	for k in range(1, n + 1):
		weight_table[k] = []

	for e in xrange(1, pow(2, n)):
		k = sum([ bool(e & (1 << i)) for i in range(n) ])
		weight_table[k].append(e)

	for f in range(n):
		for k in range(1, n + 1):
			res = sum([ abs(walsh_spectrum_table[f][e]) for e in weight_table[k] ])
			if res == 0:
				k_balance[f] = k
			else:
				break

	return k_balance

@exception_checker
@timer
def analyze_rate_distribution_error(bool_function):

	truth_table = bool_function.truth_table
	n = bool_function.field_extension

	rate_distribution_error_multivariate = []
	for F in range(n):
		rate_distribution_error_multivariate.append(0)

	rate_distribution_error_univariate = []
	for f in range(n):
		rate_distribution_error_univariate.append([ 0 for i in range(n) ])

	for x in range(pow(2, n)):
		for i in range(n):
			res = truth_table[x] ^ truth_table[x ^ (1 << i)]
			for f in range(n):
				bit = int(bool(res & (1 << f)))
				rate_distribution_error_multivariate[i] += bit
				rate_distribution_error_univariate[f][i] += bit

	return dict(multi=rate_distribution_error_multivariate, uni=rate_distribution_error_univariate)

@exception_checker
@timer
def analyze_relative_deviation_of_rde(bool_function, rate_distribution_error_multivariate, rate_distribution_error_univariate):

	n = bool_function.field_extension
	average = pow(2, n-1)

	relative_deviation_univariate = []
	for f in range(n):
		relative_deviation_univariate.append([ 0 for i in range(n) ])

	relative_deviation_multivariate = []
	for i in range(n):
		relative_deviation_multivariate.append(0)

	for f in range(n):
		for i in range(n):
			relative_deviation_univariate[f][i] = (abs(rate_distribution_error_univariate[f][i] - average) / average) * 100


	for i in range(n):
		relative_deviation_multivariate[i] = (abs(rate_distribution_error_multivariate[i] - n * average) / (n * average)) * 100

	return dict(multi=relative_deviation_multivariate, uni=relative_deviation_univariate)
	
@exception_checker
@timer
def analyze_maximum_differential_probability(bool_function):

	manager = Manager()

	n = bool_function.field_extension

	mdp1 = manager.Value('d', 0.0)
	mdp2 = manager.Value('d', 0.0)
	mdp3 = manager.Value('d', 0.0)
	mdp4 = manager.Value('d', 0.0)

	process1 = Process(target=process_analyze_mdp, args=(bool_function.truth_table, n, 1, pow(2, n - 2),  mdp1))
	process2 = Process(target=process_analyze_mdp, args=(bool_function.truth_table, n, pow(2, n - 2), pow(2, n - 1), mdp2))
	process3 = Process(target=process_analyze_mdp, args=(bool_function.truth_table, n, pow(2, n - 1), pow(2, n - 1) + pow(2, n - 2),  mdp3))
	process4 = Process(target=process_analyze_mdp, args=(bool_function.truth_table, n, pow(2, n - 1) + pow(2, n - 2), pow(2, n), mdp4))

	process1.start()
	process2.start()
	process3.start()
	process4.start()
	
	process1.join()
	process2.join()
	process3.join()
	process4.join()

	return max(mdp1.value, mdp2.value, mdp3.value, mdp4.value)

def process_analyze_mdp(truth_table, field_extension, begin, end, mdp):

	MDP = 0.0
	field_power = pow(2, field_extension)
	field_elements = range(field_power)
	a_array = range(begin, end)

	for a in a_array:
		b = [ 0 for i in field_elements ]
		
		for x in field_elements:
			b[truth_table[x] ^ truth_table[x ^ a]] += 1

		MDP = max(MDP, max(b))

	mdp.value = MDP / field_power

'''WRITE ANALYZE INFO INTO FILE'''
def write_algebraic_degree_analyze_to(filename, algebraic_degree_list, time):

	with codecs.open(filename, 'w', 'utf-8-sig') as file:
		for f in range(len(algebraic_degree_list)):
			file.write("f%d" % (f + 1) + '\t' + str(algebraic_degree_list[f]) + '\n')
		
		file.write("\nF" + '\t' + str(min(algebraic_degree_list)) + '\n')
		file.write('time: %f' % time)

def write_disbalance_analyze_to(filename, disbalance_list, time):

	with codecs.open(filename, 'w', 'utf-8-sig') as file:
		for f in range(len(disbalance_list)):
			file.write("f%d" % (f + 1) + '\t' + str(disbalance_list[f]) + '\n')

		file.write('time: %f' % time)

def write_nonlinearity_analyze_to(filename, nonlinearity_list, time):

	with codecs.open(filename, 'w', 'utf-8-sig') as file:
		for f in range(len(nonlinearity_list)):
			file.write("f%d" % (f + 1) + '\t' + str(nonlinearity_list[f]) + '\n')

		file.write('time: %f' % time)

def write_correlation_immunity_analyze_to(filename, k_balance_list, time):

	with codecs.open(filename, 'w', 'utf-8-sig') as file:
		for f in range(len(k_balance_list)):
			file.write("f%d" % (f + 1) + '\t' + str(k_balance_list[f]) + '\n')

		file.write('time: %f' % time)

def write_mdp_analyze_to(filename, mdp, time):

	with codecs.open(filename, 'w', 'utf-8-sig') as file:
		file.write(str(mdp) + '\n')
		file.write('time: %f' % time)

def write_rde_analyze_to(filename, rde, time):

	rde_uni = rde['uni']
	rde_multi = rde['multi']

	with codecs.open(filename, 'w', 'utf-8-sig') as file:
		for f in range(len(rde_uni)):
			string_to_write = ''

			for i in range(len(rde_uni)):
				string_to_write += '\t' + str(rde_uni[f][i])

			file.write("f%d" % (f + 1) + string_to_write + '\n')

		string_to_write = ''
		for i in range(len(rde_multi)):
			string_to_write += '\t' + str(rde_multi[i]) 

		file.write('\nF' + string_to_write + '\n')
		file.write('time: %f' % time)

def write_relative_deviation_analyze_to(filename, relative_deviation, time):

	relative_deviation_uni = relative_deviation['uni']
	relative_deviation_multi = relative_deviation['multi']

	with codecs.open(filename, 'w', 'utf-8-sig') as file:
		for f in range(len(relative_deviation_uni)):
			string_to_write = ''

			for i in range(len(relative_deviation_uni)):
				string_to_write += '\t' + str(relative_deviation_uni[f][i])

			file.write("f%d" % (f + 1) + string_to_write + '\n')

		string_to_write = ''
		for i in range(len(relative_deviation_multi)):
			string_to_write += '\t' + str(relative_deviation_multi[i]) 

		file.write('\nF' + string_to_write + '\n')
		file.write('time: %f' % time)