# gamma(n + 1) = factorial(n)
# gamma(n + 1) = integration(e^-x * x^n),where n is any positive real number.

import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np
import scipy.special

def integrand(x,n):
	return np.exp(-x)*(x**n)
	
def graph(list_1,arr,gam):
	plt.plot(arr,list_1,label = 'Integrand')
	plt.scatter(arr,gam,c = 'r',label = 'Inbuilt factorial function')
	plt.legend()
	plt.show()

def factorial(arr):
	fact = []
	for n in arr:
		y = scipy.special.factorial(n)
		fact.append(y)
	return fact

if __name__ == "__main__":
	arr = np.arange(1,7.5,0.1)
	list_1 = []
	for n in arr:
		I = quad(integrand,0,10000,args = (n))
		list_1.append(I[0])
	print(arr)
	print(list_1)
	fact = factorial(arr)
	graph(list_1,arr,fact)
	
