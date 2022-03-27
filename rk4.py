import numpy as np
import matplotlib.pyplot as plt

def f1(x,y):
    dx = y + x - x**3
    return dx

def f2(x):
    dy = -x
    return dy

def euler(h,x_o,y_o,N):
    x = x_o ; y = y_o
    x_list = [x_o] ; y_list = [y_o] ; all = []
    for i in range(0,N):
        x = x + h * f1(x,y)
        y = y + h * f2(x)
        x_list.append(x) ; y_list.append(y)
    all.extend([x_list,y_list])
    return all

def rk4(h,x_o,y_o,N):
    x = x_o ; y = y_o
    x_list = [x_o] ; y_list = [y_o] ; all = []
    for i in range(0,N):
        k1 = h * f1(x,y)
        l1 = h * f2(x)

        k2 = h * f1(x + (k1/2),y + (k1/2))
        l2 = h * f2(x + (l1/2))

        k3 = h * f1(x + (k2/2),y + (k2/2))
        l3 = h * f2(x + (l2/2))

        k4 = h * f1(x + k3 ,y + k3)
        l4 = h * f2(x + l3)

        x = x + 1/6 * (k1 + 2*(k2 + k3) + k4)
        y = y + 1/6 * (l1 + 2*(l2 + l3) + l4)
        x_list.append(x) ; y_list.append(y)
    all.extend([x_list,y_list])
    return all

def graph(sol1,sol2,sol3,sol4,t):
    fig,ax =  plt.subplots(2,2)
    fig1,ax1 = plt.subplots(2,2)
    fig2,ax2 = plt.subplots(2,2)
    ax[0,0].scatter(sol1[0],sol1[1])
    ax[0,1].scatter(sol2[0],sol2[1])
    ax[1,0].scatter(sol3[0],sol3[1])
    ax[1,1].scatter(sol4[0],sol4[1])

    ax1[0,0].scatter(t,sol1[0])
    ax1[0,1].scatter(t,sol2[0])
    ax1[1,0].scatter(t,sol3[0])
    ax1[1,1].scatter(t,sol4[0])

    ax2[0,0].scatter(t,sol1[1])
    ax2[0,1].scatter(t,sol2[1])
    ax2[1,0].scatter(t,sol3[1])
    ax2[1,1].scatter(t,sol4[1])
    plt.show()

if __name__ == "__main__":
    # 1st Part
    x_o = 0
    y_o = [-1,-2,-3,-4]
    N1 = int(input("Enter the value of nodal points: "))
    h1 = (15 - 0)/N1
    t1 = np.linspace(0,15,N1+1)
    # 2nd Part
    # CALLING FUNCTION FOR CALCULATIONS
    sol1 = rk4(h1,x_o,y_o[0],N1)
    sol2 = rk4(h1,x_o,y_o[1],N1)
    sol3 = rk4(h1,x_o,y_o[2],N1)
    sol4 = rk4(h1,x_o,y_o[3],N1)
    euler1 = euler(h1,x_o,y_o[0],N1)
    #graph(sol1,sol2,sol3,sol4,t1)

    plt.plot(sol1[0],sol1[1],c = "r")
    plt.plot(euler1[0],euler1[1])
    plt.show()
'''
def rk4_1(x_o,v_o,h,N,b,n):
    x = x_o ; v = v_o ; all = []
    dis = [x] ; vel = [v] 
    for i in range(0,N):
        X = [x,v]
        k1 = h * q1(X,b,n)
        X = [x + k1[0]/2,v + k1[1]/2]
        k2 = h * q1(X,b,n)
        X = [x + k2[0]/2,v + k2[1]/2]
        k3 = h * q1(X,b,n)
        X = [x + k3[0],v + k3[1]]
        k4 = h * q1(X,b,n)
        
        x = x + 1/6 * (k1[0] + 2*(k2[0] + k3[0]) + k4[0])
        v = v + 1/6 * (k1[1] + 2*(k2[1] + k3[1]) + k4[1])
        dis.append(x) ; vel.append(v)
    all.extend([dis,vel])
    return all
'''

    