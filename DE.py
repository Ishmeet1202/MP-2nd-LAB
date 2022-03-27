import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable

def nodal(a,b,h):
	N = (b - a)/h
	return N

def slope(c,tau):
	f = -(c/tau)
	return f

def analytic_sol(y_o,t,tau):
	anl_sol = []
	for i in range(len(t)):
		y = y_o * np.exp(-t[i]/tau)
		anl_sol.append(y)
	return anl_sol

def euler(h,c,tau,N):
	x = c
	y = [x]
	for i in range(0,N):
		x = x + h * slope(x,tau)
		y.append(x)
	return y
	
def RungeKutta_2(h,c,tau,N):
	x = c
	y = [x]
	for i in range(0,N):
		k1 = h * slope(x,tau)
		k2 = h * slope((x + k1),tau)
		k = (k1 + k2)/2
		x = x + k
		y.append(x)
	return y

def absolute_error(arr1,arr2):
	err = np.array(arr1) - np.array(arr2)
	return err

def graph(t_r1,t_rc,t_st,e1,e2,e3,r1,r2,r3,e_1a,e_1b,e_2a,e_2b,e_3a,e_3b):
	fig,ax =  plt.subplots(1,2)
	fig1,ax1 = plt.subplots(1,2)
	fig2,ax2 = plt.subplots(1,2) 
	ax[0].scatter(t_r1,e1,label = "Euler method")
	ax[0].scatter(t_r1,r1,c = "orange",label = "Runge Kutta order 2",s = 50,marker = "*")
	ax[1].scatter(t_r1,e_1a,label = "Absolute Error in euler")
	ax[1].scatter(t_r1,e_1b,c = "orange",label = "Absolute Error in RK2")
	ax1[0].scatter(t_rc,e2,label = "Euler method")
	ax1[0].scatter(t_rc,r2,c = "orange",label = "Runge Kutta order 2",s = 50,marker = "*")
	ax1[1].scatter(t_rc,e_2a,label = "Absolute Error in euler")
	ax1[1].scatter(t_rc,e_2b,c = "orange",label = "Absolute Error in RK2")
	ax2[0].scatter(t_st,e3,label = "Euler method")
	ax2[0].scatter(t_st,r3,c = "orange",label = "Runge Kutta order 2",s = 50,marker = "*")
	ax2[1].scatter(t_st,e_3a,label = "Absolute Error in euler")
	ax2[1].scatter(t_st,e_3b,c = "orange",label = "Absolute Error in RK2")
	for i in range(2):
		if i == 0:
			ax[i].set(xlabel = "Time / Half life",ylabel = "No. of Nuclei (N)",title = "Method's Graph")
			ax1[i].set(xlabel = "Time (second)",ylabel = "Voltage (V)",title = "Method's Graph",xlim = (-0.0005,0.0055))
			ax2[i].set(xlabel = "Time (second)",ylabel = "Velocity (v)",title = "Method's Graph")
		else:
			ax[i].set(xlabel = "Time / Half life",ylabel = "Absolute error",title = "Absolute Error Graph")
			ax1[i].set(xlabel = "Time (second)",ylabel = "Absolute error",title = "Absolute Error Graph",xlim = (-0.0005,0.0055))
			ax2[i].set(xlabel = "Time (second)",ylabel = "Absolute error",title = "Absolute Error Graph")
		ax[i].grid(ls = "--")
		ax1[i].grid(ls = "--")
		ax2[i].grid(ls = "--")
		ax[i].legend()
		ax1[i].legend()
		ax2[i].legend()
	fig.suptitle("RADIOACTIVE DECAY")
	fig1.suptitle("RC CIRCUIT")
	fig2.suptitle("STOKE'S LAW")
	plt.show()

def table(t_r1,t_rc,t_st,e1,e2,e3,r1,r2,r3,e_1a,e_1b,e_2a,e_2b,e_3a,e_3b,a1,a2,a3):
	table1 = PrettyTable(['Sr No.','Time / Half life','No. of nuclei (Analytical)','No. of nuclei (Euler)','No. of nuclei (Runge kutta order 2)','Abs Error(Euler)','Abs Error(RK2)'])
	table2 = PrettyTable(['Sr No.','Time (second)','Voltage (Analytical)','Voltage (Euler)','Voltage (Runge kutta order 2)','Abs Error(Euler)','Abs Error(RK2)'])
	table3 = PrettyTable(['Sr No.','Time (second)','Velocity (Analytical)','Velocity (Euler)','Velocity (Runge kutta order 2)','Abs Error(Euler)','Abs Error(RK2)'])
	for i in range(len(t_r)):
		table1.add_row([i+1,t_r1[i],a1[i],e1[i],r1[i],e_1a[i],e_1b[i]])
		table2.add_row([i+1,t_rc[i],a2[i],e2[i],r2[i],e_2a[i],e_2b[i]])
		table3.add_row([i+1,t_st[i],a3[i],e3[i],r3[i],e_3a[i],e_3b[i]])
	print("\nRADIOACTIVE DECAY:\n",table1,"\n")
	print("RC CIRCUIT:\n",table2,"\n")
	print("STOKE'S LAW:\n",table3,"\n")
	
if __name__ == "__main__":
	# INITIAL CONDITION
	half_t = 4                                         # q1
	R = 1000 ; C = 10**-6					     	   # q2	
	m = 50 ;r = 4 *10**-2 ;neta = 20		     	   # q3
	tau = np.array([1/(0.693/half_t),R*C,m/(6*np.pi*neta*r)]) 
	c = [20000,10,5]
	# CALCULATIONS FOR STEP SIZE AND TIME AXIS
	a = 0
	b = np.array([5 * tau[0],5 * tau[1],5 * tau[2]])
	h = np.array([tau[0]/10,tau[1]/10,tau[2]/10])		   
	N_r = int(nodal(a,b[0],h[0])) ; N_rc = int(nodal(a,b[1],h[1])) ; N_st = int(nodal(a,b[2],h[2]))
	t_r = np.linspace(0,b[0],51)
	t_r1 = (np.linspace(0,b[0],51))/half_t
	t_rc = (np.linspace(0,b[1],51))
	t_st = (np.linspace(0,b[2],51)) 
	# CALLING FUNCTIONS FOR CALCULATIONS
	e1 = euler(h[0],c[0],tau[0],N_r)
	e2 = euler(h[1],c[1],tau[1],N_rc)
	e3 = euler(h[2],c[2],tau[2],N_st)
	r1 = RungeKutta_2(h[0],c[0],tau[0],N_r)
	r2 = RungeKutta_2(h[1],c[1],tau[1],N_rc)
	r3 = RungeKutta_2(h[2],c[2],tau[2],N_st)
	a1 = analytic_sol(c[0],t_r,tau[0])
	a2 = analytic_sol(c[1],t_rc,tau[1])
	a3 = analytic_sol(c[2],t_st,tau[2])
	e_1a = absolute_error(a1,e1) ; e_1b = absolute_error(a1,r1)
	e_2a = absolute_error(a2,e2) ; e_2b = absolute_error(a2,r2)
	e_3a = absolute_error(a3,e3) ; e_3b = absolute_error(a3,r3)
	table(t_r1,t_rc,t_st,e1,e2,e3,r1,r2,r3,e_1a,e_1b,e_2a,e_2b,e_3a,e_3b,a1,a2,a3)
	graph(t_r1,t_rc,t_st,e1,e2,e3,r1,r2,r3,e_1a,e_1b,e_2a,e_2b,e_3a,e_3b)
