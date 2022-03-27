import math as m

def expo(x,n):
	ex = 0
	for i in range(n+1):
		ex = ex + (x**i) / m.factorial(i)
	return ex

if __name__ == "__main__":
	n = int(input("Enter the number upto which you want to calculate series: "))
	x = float(input("Enter the value: "))
	print("Value of exponential %.3f"%expo(x,n))
	
		
		
