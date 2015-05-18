#!/usr/bin/python
# -*- coding: utf-8 -*-

from boolfunc import PolynomialOverF2, BoolFunction, timer
from analyze import *
from static import PATH

def main():

	p = PolynomialOverF2(131081)
	N = 257
	
	f = BoolFunction(N, p)

	f.read_truth_table_from(PATH['f(x)']['Truth table'])
	f.read_anf_table_from(PATH['f(x)']['ANF table'])
	f.read_walsh_spectrum_table_from(PATH['f(x)']['Walsh spectrum table'])

	res = analyze_k_balance(f)
	write_k_balance_analyze_to(PATH['f(x)']['k-balance'], res['result'], res['time'])

if __name__ == '__main__':
	main()