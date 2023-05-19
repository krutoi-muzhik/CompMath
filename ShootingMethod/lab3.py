import numpy as np
import matplotlib.pyplot as plt
import math


# y`` + 1/2 * 1/(1 - 1/2*y) * y`^2 = 0
# 0 < x <= 1
# y(0) = y0
# y(1) = 0
# y0 = 0,25; 0,5; 1; 1,5; 1,8; 1,9; 1,95;

L = 1
y0_arr = np.array([0.25, 0.5, 1, 1.5, 1.8, 1.9, 1.95])

left = 0
right = L

# u` = f(x, y, u)
# y` = u

def func (y, u):
	return u, -1/2 * u**2 / (1 - 1/2 * y)



def methodRungeKutta4 (f, h, y0, u0):
	N = int((right - left) / h)
	t = np.zeros(N)
	for i in range(N):
		t[i] = left + (right - left) * i / h

	y = np.zeros(N)
	u = np.zeros(N)

	y[0] = y0
	u[0] = u0

	for i in range(1, N):
		k1y, k1u = f(y[i-1], u[i-1])
		k2y, k2u = f(y[i-1] + h*k1y/2, u[i-1] + h*k1u/2)
		k3y, k3u = f(y[i-1] + h*k2y/2, u[i-1] + h*k2u/2)
		k4y, k4u = f(y[i-1] + h*k3y, u[i-1] + h*k3u)
		y[i] = y[i-1] + h*(k1y + 2*k2y + 2*k3y + k4y)/6
		u[i] = u[i-1] + h*(k1u + 2*k2u + 2*k3u + k4u)/6

	return t, y, u


def searchBinary (f, y0, yL, h, epsilon, al, ar):
	yLfinal = yL + 10*epsilon

	while (abs(yLfinal - yL) > epsilon):
		alpha = (al + ar) / 2
		# print("left = ", al, "right = ", ar)
		# print("alpha = ", alpha)
		u0 = alpha
		t, y, u = methodRungeKutta4(f, h, y0, u0)

		yLfinal = y[y.size - 1]

		#print("yLfinal = ", yLfinal)

		if yLfinal - yL < 0:
			print("1")
			al = alpha
		else:
			print("2")
			ar = alpha

	print("y0 = ", y0)
	print("alpha = ", alpha)
	print("yLfinal = ", yLfinal)

	plt.plot(t, u)
	plt.show()

def main ():
	h = 0.001
	eps = 0.001

	for i in range(y0_arr.size):
		searchBinary(func, y0_arr[i], 0, h, eps, -10, 10)

main()