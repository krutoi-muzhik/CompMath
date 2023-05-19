import numpy as np
import matplotlib.pyplot as plt

x0 = 0.125
N = 10001
h = 1 / N


def k_left (x):
	return x + 1

def k_right (x):
	return x ** 2

def k (n):
	if n * h <= x0:
		return k_left(n * h)
	if n * h > x0:
		return k_right(n * h)

def q (n):
	return (n * h) ** 2

def f (n):
	return np.cos(n * h)


# Coefficients in equation
# a(n)*y(n-1) + b(n)*y(n) + c(n)y(n+1) = d(n)
# Given equation:
# 1/h*(k(n+1/2)*(y(n+1)-y(n))/h - k(n-1/2)*(y(n)-y(n-1))/h) - q(n)*y(n) = -f(n)

def a (n):
	return k(n - 1/2)

def b (n):
	return -(k(n + 1/2) + k(n - 1/2)) - q(n) * (h ** 2)

def c (n):
	return k(n + 1/2)

def d (n):
	return -f(n) * (h ** 2)


def main ():
	x = np.zeros(N)
	for i in range(x.size):
		x[i] = i * h

	u = np.zeros(N)
	u[0] = 1 										# Boundary
	u[N - 1] = 0 									# conditions

	alpha = np.zeros(N)
	betta = np.zeros(N)

	# Equation
	# a(n)*y(n-1) + b(n)*y(n) + c(n)*y(n+1) = d(n)
	# Starting conditions u[0] = 1 u[N-1] = 0 can be represented as equalities
	# a(1) + b(1)*y(1) + c(1)*y(2) = d(1)
	# a(N)*y(N-3) + b(N)y(N-2) = d(N)
	# This system can be represented as a matrix
	#  [b0	c0	0  ...	0 ]
	#  [a1	b1	c1 ...	0 ]
	#  [0	a2	b2  ...	0 ]
	#  [...	... ... ...   ]
	#  [0 	0 	0 	... bN]
	
	n0 = int(x0 * N)

	# [0------------x0------------------1]
	#
	# Deviding segment [0, N-1] to 2 smaller segments by n0
	#
	# In the left segment
	# Using this system every Yn can be represented by the next Yn+1, considering the previous Yn-1 as a known one
	#
	# y(n) = alpha[n]*y(n+1) + betta[n]
	# y(n) = -c(n)/b(n)*y(n+1) + (d(n) - a(n)*y(n-1))/b(n)
	# y(1) = -c(1)/b(1) + d(n)/b(n)
	# 
	# After substitution the starting equations transforms to
	# c(n)*y(n+1) + b(n)*y(n) + a(n)*(y(n)*alpha(n-1) + betta(n-1)) = d(n)
	# y(n) = -c(n)/(b(n) + a(n)*alpha(n-1))*y(n+1) + (d(n) - a(n)*betta(n-1))/(b(n) + a(n)*alpha(n-1))
	# Alpha and betta are eeasily represented by the previous alpha and betta
	#
	# In the right segment
	# Using this system every Yn can be represented by the previous Yn-1, considering the previous Yn+1 as a known one
	#
	# y(n) = alpha[n]*y(n-1) + betta[n]
	# y(n) = -a(n)/b(n)*y(n-1) + (d(n) - c(n)*y(n+1))/b(n)
	# y(N-2) = -a(N-2)/b(N-2) + d(N-2)/b(N-2)
	
	alpha[0] = 0 										# No
	betta[0] = u[0] 									# need
	alpha[N-1] = 0 										# in
	betta[N-1] = u[N-1]									# them


	for n in range(1, n0):
		alpha[n] = c(n) / (-b(n) - alpha[n-1] * a(n))
		betta[n] = (d(n) - a(n) * betta[n-1]) / (b(n) + a(n) * alpha[n-1])

	for n in range(N-2, n0+1, -1):
		alpha[n] = a(n) / (-b(n) - alpha[n+1] * c(n))
		betta[n] = (d(n) - c(n) * betta[n-1]) / (b(n) + c(n) * alpha[n+1])
 

	# u[n0] = betta[n0+2] / (1 - alpha[n0+2])
	# u[n0+1] = u[n0]
	# u[n0+2] = u[n0]
	# u[n0-1] = alpha[n0-1] * u[n0] + betta[n0-1]

	u[n0] = (k(n0+1)*betta[n0+2] + k(n0)*betta[n0-1])/(k(n0)*(1-alpha[n0-1]) + k(n0+1)*(1-alpha[n0+2]))
	u[n0+1] = u[n0]
	u[n0-1] = alpha[n0-1]*u[n0] + betta[n0-1]
	u[n0+2] = alpha[n0+2]*u[n0+1] + betta[n0+2]


	for i in range(n0-2, 0, -1):
		u[i] = alpha[i] * u[i+1] + betta[i]

	for i in range(n0+3, N-1):
		u[i] = alpha[i] * u[i-1] + betta[i]


	plt.plot(x, u)
	plt.show()


	# plt.plot(x, alpha)
	# plt.show()


	# plt.plot(x, betta)
	# plt.show()



main()
