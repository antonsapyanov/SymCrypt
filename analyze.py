#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division

from static import PATH
from boolfunc import PolynomialOverF2, BoolFunction, timer
from multiprocessing  import Process, Manager

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

	n = bool_function.field_extension

	mdp1 = manager.Value('d', 0.0)
	mdp2 = manager.Value('d', 0.0)

	process1 = Process(target=process_analyze_mdp, args=(bool_function.truth_table, n, 1, pow(2, n - 1),  mdp1))
	process2 = Process(target=process_analyze_mdp, args=(bool_function.truth_table, n, pow(2, n - 1), pow(2, n), mdp2))

	process1.start()
	process2.start()
	
	process1.join()
	process2.join()

	return max(mdp1.value, mdp2.value)

def process_analyze_mdp(truth_table, field_extension, begin, end, mdp):

	MDP = 0.0
	power = pow(2, field_extension)
	power_array = range(power)
	a_array = range(begin, end)

	for a in a_array:
		b = [ 0 for i in power_array ]
		
		for x in power_array:
			b[truth_table[x] ^ truth_table[x ^ a]] |= 1

		MDP = max(MDP, max(b))

	mdp.value = MDP / power

def analyze_f():

	p = PolynomialOverF2(131081)
	N = 257
	
	f = BoolFunction(N, p)

	f.read_truth_table_from(PATH['f(x)']['Truth table'])
	f.read_anf_table_from(PATH['f(x)']['ANF table'])
	f.read_walsh_spectrum_table_from(PATH['f(x)']['Walsh spectrum table'])

	print analyze_algebraic_degree(f)
	print analyze_disbalance(f)
	print analyze_nonlinearity(f)

	rate_distribution_error = analyze_rate_distribution_error(f)['result']
	print analyze_relative_deviation_of_rde(
		f,
		rate_distribution_error['multi'],
		rate_distribution_error['uni'])['result']['uni']

	'''
	mdp =  analyze_maximum_differential_probability(f)

	with codecs.open(PATH['f(x)']['mdp.txt'], 'w', 'utf-8-sig') as file:
		file.write(str(mdp))
	'''
	

def analyze_g():

	p = PolynomialOverF2(131081)
	M = 258

	g = BoolFunction(M, p)

	g.create_truth_table()
	g.write_truth_table_to(PATH['g(x)']['Truth table'])

	g.create_anf_table()
	g.write_anf_table_to(PATH['g(x)']['ANF table'])

	g.create_walsh_spectrum_table()
	g.write_walsh_spectrum_table_to(PATH['g(x)']['Walsh spectrum table'])

	print analyze_algebraic_degree(g)
	print analyze_disbalance(g)
	print analyze_nonlinearity(g)

	rate_distribution_error = analyze_rate_distribution_error(g)['result']
	print analyze_relative_deviation_of_rde(
		g,
		rate_distribution_error['multi'],
		rate_distribution_error['uni'])['result']['uni']

def analyze_prev_f():

	p = PolynomialOverF2(131081)
	N = 13
	
	f = BoolFunction(N, p)

	f.create_truth_table()
	f.create_anf_table()
	f.create_walsh_spectrum_table()

	# print analyze_algebraic_degree(f)
	# print analyze_disbalance(f)
	# print analyze_nonlinearity(f)

	rate_distribution_error = analyze_rate_distribution_error(f)['result']
	print analyze_relative_deviation_of_rde(
		f,
		rate_distribution_error['multi'],
		rate_distribution_error['uni'])['result']['uni']

def main():

	analyze_prev_f()
	# analyze_g()

if __name__ == '__main__':
	main()