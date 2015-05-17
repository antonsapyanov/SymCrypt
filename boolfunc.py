#!/usr/bin/python
# -*- coding: utf-8 -*-

import time
import codecs

def timer(function):

	def wrapper(*args, **kwargs):

		print "Виклик функції %s (%d)..." % (function.__name__, id(function))
		t0 = time.time()
		result = function(*args, **kwargs)
		time_delta = time.time() - t0
		print "Час виконання фунції %s (%d): %f сек." % (function.__name__, id(function), time_delta)

		return dict(result=result, time=time_delta)

	return wrapper

class PolynomialOverF2(object):

	def __init__(self, number):
		
		if not isinstance(number, (int, long)):
			raise TypeError("number must be an isinstance of int or long")

		self._number = number
		self._deg = len(bin(number)[2:]) - 1

	def __str__(self):
		
		return '< ' + ', '.join([ str(int(bool(self._number & (1 << i)))) for i in range(self._deg + 1) ]) + ' >'

	def repr(self, n=0):
		
		string_coef_array = bin(self._number)[2:]
		k = n - len(string_coef_array)

		return tuple([ int(bool(self._number & (1 << i))) for i in range(self._deg + 1) ] + [ 0 for i in range(k)])

	'''to deal with as a number'''
	def __lshift__(self, value):
		
		if not isinstance(value, (int, long)):
			raise TypeError("value must be an isinstance of int or long")

		number = self._number << value
		
		return PolynomialOverF2(number)

	def __rshift__(self, value):
		
		if not isinstance(value, (int, long)):
			raise TypeError("value must be an isinstance of int or long")

		number = self._number >> value
		
		return PolynomialOverF2(number)

	def __eq__(self, other):

		if isinstance(other, PolynomialOverF2):
			return self._number == other._number
		else:
			raise TypeError("other must be an isinstance of PolynomialOverF2")

	def __ne__(self, other):
		
		return not eq(self, other)

	def is_(self, other):
		
		return eq(self, other)

	def is_not(self, other):
		
		return ne(self, other)

	def __add__(self, other):
		
		if isinstance(other, PolynomialOverF2):
			return PolynomialOverF2(self._number ^ other._number)
		else:
			raise TypeError("other must be an isinstance of PolynomialOverF2")

	def __sub__(self, other):
		
		return add(self, other)

	def __mod__(self, other):
		
		if isinstance(other, PolynomialOverF2):
			number = self._number
			k = self._deg
			n = other._deg

			while k >= n:
				number ^= other._number << (k - n)
				k = len(bin(number)[2:]) - 1

 			return PolynomialOverF2(number)
		else:
			raise TypeError("other must be an isinstance of PolynomialOverF2")

	def __mul__(self, other):
		
		if isinstance(other, PolynomialOverF2):
			hx = 0

			for i, coef in enumerate(other.repr()):
				if coef:
					hx ^= self._number << i

 			return PolynomialOverF2(hx)
		else:
			raise TypeError("other must be an isinstance of PolynomialOverF2")

	def scalar_mul(self, other):

		if isinstance(other, PolynomialOverF2):
			return sum([ int(bit) for bit in bin(self._number & other._number)[2:] ]) % 2
		else:
			raise TypeError("other must be an isinstance of PolynomialOverF2")

	def __div__(self, other):
		raise Exception("div is not implemented")

	def __pow__(self, value, modulo):
	
		if not isinstance(value, (int, long)):
			raise TypeError("value must be an isinstance of int or long")
		if not isinstance(modulo, PolynomialOverF2):
			raise TypeError("modulo must be an isinstance of PolynomialOverF2")

		# Horner's method

		b = PolynomialOverF2(1)
		c = self
		m = tuple([ int(bit) for bit in reversed(bin(value)[2:]) ])
		s = len(m)

		for i in range(s):
			if m[i] == 1:
				b = (b * c) % modulo
			c = (c * c) % modulo

		return b

	'''to deal with as an array'''
	def __getitem__(self, index):
		
		if not isinstance(key, (int, long)):
			raise TypeError("key must be an isinstance of int or long")

		return self._number & (1 << index)

	def __iter__(self):
		
		for coef in self.repr():
			yield coef

	def next(self):
		
		return iter(self).next()

	'''to deal with as a function'''
	def __call__(self, x):
		
		if x not in (0, 1):
			raise ValueError("x must a bit 0 or 1")

		return sum([ coef ^ x for coef in self.repr()[1:] ]) % 2 ^ self.repr()[0]

	'''Properties'''
	@property
	def number(self):

		return self._number

	@property
	def deg(self):

		return self._deg

	@property
	def weight(self):

		return sum(self.repr())

class BoolFunction(object):

	def __init__(self, power, generator):
		
		if not isinstance(power, (int, long)):
			raise TypeError("power must be an isinstance of int or long")
		if not isinstance(generator, PolynomialOverF2):
			raise TypeError("generator must be an isinstance of PolynomialOverF2 class")

		self._power = power
		self._generator = generator

		self._truth_table = []
		self._anf_table = []
		self._walsh_spectrum_table = []

		for i in range(self._generator.deg):
			self._walsh_spectrum_table.append([])

	def __call__(self, polynomial):
		if not isinstance(polynomial, PolynomialOverF2):
			raise TypeError("polynomial must be an isinstance of PolynomialOverF2 class")

		return pow(polynomial, self._power, self._generator)

	'''Truth table'''
	@timer
	def create_truth_table(self):

		for i in xrange(pow(2, self._generator.deg)):
			element = self(PolynomialOverF2(i))
			self._truth_table.append(element.number)

	def write_truth_table_to(self, filename):
		
		with codecs.open(filename, 'w', 'utf-8-sig') as file:
			for i in xrange(pow(2, self._generator.deg)):
				file.write(str(i) + '\t' + str(self._truth_table[i]) + '\n')

	def read_truth_table_from(self, filename):
		
		with codecs.open(filename, 'r', 'utf-8-sig') as file:
			for i, line in enumerate(file):
				self._truth_table.append(int(line.rsplit('\t')[1]))

	'''ANF table'''
	@timer
	def create_anf_table(self):

		if len(self._truth_table) == 0:
			raise Exception("truth table is empty")

		for u in self._truth_table:
			self._anf_table.append(u)

		for i in range(self._generator.deg):

			already_used = []

			for e in xrange(pow(2, self._generator.deg)):
				if e in already_used:
					already_used.remove(e)
					continue

				u0 = e
				u1 = e | (1 << i)
				
				self._anf_table[u1] = self._anf_table[u0] ^ self._anf_table[u1]

				already_used.append(u1)

	def write_anf_table_to(self, filename):

		with codecs.open(filename, 'w', 'utf-8-sig') as file:
			for i in xrange(pow(2, self._generator.deg)):
				file.write(str(i) + '\t' + str(self._anf_table[i]) + '\n')

	def read_anf_table_from(self, filename):
		
		with codecs.open(filename, 'r', 'utf-8-sig') as file:
			for i, line in enumerate(file):
				self._anf_table.append(int(line.rsplit('\t')[1]))

	'''Walsh spectrum table'''
	@timer
	def create_walsh_spectrum_table(self):

		if len(self._truth_table) == 0:
			raise Exception("truth table is empty")

		for u in self._truth_table:
			for i in range(self._generator.deg):
				self._walsh_spectrum_table[i].append(pow(-1, int(bool(u & (1 << i)))))

		for i in range(self._generator.deg):

			already_used = []

			for e in xrange(pow(2, self._generator.deg)):
				if e in already_used:
					already_used.remove(e)
					continue

				u0 = e
				u1 = e | (1 << i)
				
				for j in range(self._generator.deg):
					Cju0 = self._walsh_spectrum_table[j][u0]
					Cju1 = self._walsh_spectrum_table[j][u1]

					tmp0 = Cju0 + Cju1
					tmp1 = Cju0 - Cju1

					self._walsh_spectrum_table[j][u0] = tmp0
					self._walsh_spectrum_table[j][u1] = tmp1

				already_used.append(u1)

	def write_walsh_spectrum_table_to(self, filename):

		with codecs.open(filename, 'w', 'utf-8-sig') as file:
			for i in xrange(pow(2, self._generator.deg)):
				string_to_write = ''

				for j in range(self._generator.deg):
					string_to_write += '\t' + str(self._walsh_spectrum_table[j][i])

				file.write(str(i) + string_to_write + '\n')

	def read_walsh_spectrum_table_from(self, filename):
		
		with codecs.open(filename, 'r', 'utf-8-sig') as file:
			for i, line in enumerate(file):
				self._walsh_spectrum_table.append(int(line.rsplit()[1]))

		with codecs.open(filename, 'r', 'utf-8-sig') as file:
			for i, line in enumerate(file):
				spectra = line.rsplit('\t')[1:]

				for j in range(len(spectra)):
					self._walsh_spectrum_table[j].append(int(spectra[j]))

	'''Properties'''
	@property
	def truth_table(self):

		return self._truth_table[:]

	@property
	def anf_table(self):

		return self._anf_table[:]

	@property
	def walsh_spectrum_table(self):

		return self._walsh_spectrum_table[:]

	@property
	def  field_extension(self):

		return self._generator.deg