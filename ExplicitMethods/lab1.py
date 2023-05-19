import numpy as np
import matplotlib.pyplot as plt
import math

# y`` + e*(y^2-1)*y` + y = 0
#
# x` = z
# z` = e*(1-x^2)*z - x
# x(0) = 2
# z(0) = 0
# e > 0
# 0 < t <= 100

e = 10
tl = 0 			# left t boundary
tr = 100 		# right t boundary
z0 = 0
x0 = 2


def func1 (x, z):
	return z, e * (1 - x * x) * z - x

def methodRungeKutta4 (f, h):
	N = int((tr - tl) / h)
	t = np.zeros(N)
	for i in range(N):
		t[i] = tl + (tr - tl) * i / h

	z = np.zeros(N)											# 2 dimensional system
	x = np.zeros(N)

	z[0] = z0
	x[0] = x0

	# Butcher tableau for Runge-Kutta method of 4th order
	# [0.5 0.5 	0 	0 	0 ]
	# [0.5	0  0.5	0 	0 ]
	# [	1 	0 	0 	1 	0 ]
	# [	   1/6 2/6 2/6 1/6]

	# From Butcher tableau
	# y(n+1) = y(n) + h*(k1 + 2*k2 + 2*k3 + k4)/6
	# k1 = f(t(n), y(n))
	# k2 = f(t(n) + h/2, y(n) + h*k1/2)
	# k3 = f(t(n) + h/2, y(n) + h*k2/2)
	# k4 = f(t(n) + h, y(n) + h*k3)
	# In this case f(y) has no dependence on t


	for i in range(1, N - 1):
		k1x, k1z = f(x[i-1], z[i-1])
		k2x, k2z = f(x[i-1] + h*k1x/2, z[i-1] + h*k1z/2)
		k3x, k3z = f(x[i-1] + h*k2x/2, z[i-1] + h*k2z/2)
		k4x, k4z = f(x[i-1] + h*k3x, z[i-1] + h*k3z)
		x[i] = x[i-1] + h*(k1x + 2*k2x + 2*k3x + k4x)/6
		z[i] = z[i-1] + h*(k1z + 2*k2z + 2*k3z + k4z)/6

	plt.plot(t, z)
	plt.show()

	return x, z


def methodAdams4 (f, h, initx, initz):
	N = int((tr - tl) / h)
	t = np.zeros(N)
	for i in range(N):
		t[i] = tl + (tr - tl) * i / h

	x = np.zeros(N)
	z = np.zeros(N)

	for i in range(4):				# cnt of dots is got from RungeKutta
		x[i] = initx[i]
		z[i] = initz[i]

	# Adams method
	# In the next expression f(i) is f(t(i), y(i))
	# In this case f(i) is f(y(i)) as t is not needed
	# y(n+4) = h*(55/24*f(n+3) - 59/24*f(n+2) + 37/24*f(n+1) - 3/8f(n))
	# To get the sollution 4 startingg dots are needed
	# They can be got from RungeKutta method

	for i in range(4, N - 1):
		k1x, k1z = f(x[i-1], z[i-1])
		k2x, k2z = f(x[i-2], z[i-2])
		k3x, k3z = f(x[i-3], z[i-3])
		k4x, k4z = f(x[i-4], z[i-4])
		x[i] = x[i-1] + h * (55/24 * k1x - 59/24 * k2x + 37/24 * k3x - 3/8 * k4x)
		z[i] = z[i-1] + h * (55/24 * k1z - 59/24 * k2z + 37/24 * k3z - 3/8 * k4z)

	plt.plot(t, z)
	plt.show()

	return x, z


def main ():
	h = 0.1
	x1, z1 = methodRungeKutta4(func1, h)

	n = 4
	x2, z2 = methodAdams4(func1, h, x1[0:4], z1[0:4])


main()