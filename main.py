#!/usr/bin/python
# -*- coding: utf-8 -*-

from boolfunc import PolynomialOverF2, BoolFunction, timer
from analyze import *
from static import PATH

def main():
	''' g(x) '''
	p = PolynomialOverF2(131081)
	N = 258
	
	f = BoolFunction(N, p)
	
	f.read_truth_table_from(PATH['g(x)']['Truth table'])
	f.read_anf_table_from(PATH['g(x)']['ANF table'])
	f.read_walsh_spectrum_table_from(PATH['g(x)']['Walsh spectrum table'])

	res = analyze_algebraic_degree(f)
	write_algebraic_degree_analyze_to(PATH['g(x)']['Algebraic degree'], res['result'], res['time'])

	res = analyze_disbalance(f)
	write_disbalance_analyze_to(PATH['g(x)']['Disbalance'], res['result'], res['time'])

	res = analyze_nonlinearity(f)
	write_nonlinearity_analyze_to(PATH['g(x)']['Nonlinearity'], res['result'], res['time'])

	res = analyze_correlation_immunity(f)
	write_correlation_immunity_analyze_to(PATH['g(x)']['Correlation immunity'], res['result'], res['time'])

	res = analyze_rate_distribution_error(f)
	write_rde_analyze_to(PATH['g(x)']['RDE'], res['result'], res['time'])
	
	res = analyze_relative_deviation_of_rde(f, res['result']['multi'], res['result']['uni'])
	write_relative_deviation_analyze_to(PATH['g(x)']['Relative deviation of RDE'], res['result'], res['time'])

if __name__ == '__main__':
	main()