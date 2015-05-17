#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division

from static import PATH
from boolfunc import PolynomialOverF2, BoolFunction, timer

from multiprocessing  import Process, Manager, Lock

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
def analyze_k_balance(bool_function):

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
			res = sum([ bool(walsh_spectrum_table[f][e]) for e in weight_table[k] ])
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
	lock = manager.Lock()

	truth_table = manager.list(bool_function.truth_table)
	n = bool_function.field_extension
	differential_probability = manager.dict()
	for a in range(1, pow(2, n)):
		differential_probability[a] = manager.dict()
	
	process1 = Process(target=process_analyze_mdp, args=(lock, truth_table, differential_probability, n, 1, pow(2, n - 1)))
	process2 = Process(target=process_analyze_mdp, args=(lock, truth_table, differential_probability, n, pow(2, n - 1), pow(2, n)))
	
	process1.start()
	process2.start()
	
	process1.join()
	process2.join()

	MDP = pow(2, -n) * max([ max(differential_probability[a].values()) for a in differential_probability ])

	return MDP
	
def process_analyze_mdp(lock, truth_table, differential_probability, field_extension, begin, end):

	for x in xrange(pow(2, field_extension)):
		for a in xrange(begin, end):
			b = truth_table[x] ^ truth_table[x ^ a]

			lock.acquire()
			try:
				differential_probability[a][b] += 1
			except KeyError:
				differential_probability[a][b] = 1
			lock.release()

def analyze_f():

	p = PolynomialOverF2(131081)
	N = 13

	f = BoolFunction(N, p)

	f.read_truth_table_from(PATH['f(x)']['Truth table'])
	f.read_anf_table_from(PATH['f(x)']['ANF table'])
	f.read_walsh_spectrum_table_from(PATH['f(x)']['Walsh spectrum table'])

	print analyze_maximum_differential_probability(f)
	
	'''
	f.create_truth_table()
	f.create_anf_table()
	f.create_walsh_spectrum_table()

	f.write_truth_table_to(PATH['f(x)']['Truth table'])
	f.write_anf_table_to(PATH['f(x)']['ANF table'])
	f.write_walsh_spectrum_table_to(PATH['f(x)']['Walsh spectrum table'])
	'''


def analyze_g():
	
	# p = PolynomialOverF2(69643)
	p = PolynomialOverF2(32771)
	M = 5

	g = BoolFunction(M, p)

	# g.read_truth_table_from('M14/truth_table.txt')
	# g.read_anf_table_from('M14/anf_table.txt')
	# g.read_walsh_spectrum_table_from('M14/walsh_spectrum_table.txt')

	g.create_truth_table()
	g.create_anf_table()
	g.create_walsh_spectrum_table()

	print analyze_algebraic_degree(g)
	print analyze_disbalance(g)
	print analyze_nonlinearity(g)
	print analyze_k_balance(g)

def main():

	analyze_f()

if __name__ == '__main__':
	main()
