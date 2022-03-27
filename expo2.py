import matplotlib.pyplot as plt
import math
n = [1,2,3,4,5,6,7,8,9,10]
def exp(x,m):
  e=0
  for i in range(0,m): 
    e=e+(x**i)/math.factorial(i)
  return(e)

x = float(input("Enter No: "))
y =[]
for i in n:
  yy = exp(x,i)
  y.append(yy)

  print ('%.3f'%yy)
plt.plot(n,y,'r*')
plt.show()
