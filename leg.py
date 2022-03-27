import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre
import scipy.special
from prettytable import PrettyTable

def gamma(num):
    if num == 0.5:
        return np.sqrt(np.pi)
    elif num == 1:
        return num
    else:
        return (num-1)*gamma(num-1)

def kron(m,n):
    if m == n:
        return 1
    else:
        return 0

def legendre_inbuilt(x,n):
    leg = legendre(n)
    P_n = leg(x)
    return P_n
 
def legendre_m(n,x):
    p = 0
    if n%2 == 0:
        for i in np.arange((n/2)+1):
            p = p + (((-1)**i * scipy.special.factorial(2*n - 2*i))/(2**n * scipy.special.factorial(i) * scipy.special.factorial(n - i) * scipy.special.factorial(n - 2*i))) * x**(n - 2*i)
    else:
        for i in np.arange(((n-1)/2)+1):
            p = p + (((-1)**i * scipy.special.factorial(2*n - 2*i))/(2**n * scipy.special.factorial(i) * scipy.special.factorial(n - i) * scipy.special.factorial(n - 2*i))) * x**(n - 2*i)
    return p

def derivative(x,n):
    p_n = legendre_m(n,x)
    p_m = legendre_m(n-1,x)
    der_p = (n * (p_m - (x*p_n)))/ (1 - x**2)
    return der_p

def recursion_relation1(n,x):
    list_x = []
    list_n = []
    list_n1 = []
    data3 = []
    P_n = []
    P_d = []
    P_d1 = []
    l = []
    r = []
    all1 = []
    for j in x:
        data3.extend([j,n,(n-1),legendre_m(n,j),derivative(j,n),derivative(j,n-1)])
        list_x.append(j)
        list_n.append(n)
        list_n1.append(n-1)
        P_n.append(legendre_m(n,j))
        P_d.append(derivative(j,n))
        P_d1.append(derivative(j,n-1))
        r1_l = n * legendre_m(n,j)
        r1_r = ((j * derivative(j,n)) - derivative(j,n-1))
        l.append(r1_l)
        r.append(r1_r)
    all1.extend([list_x,list_n,list_n1,P_n,P_d,P_d1,l,r])
    data_3 =np.reshape(data3,(99,6))
    file3 = np.savetxt("/home/ishmeet/BSc physics/MP 2nd/Lab_Codes/leg02.txt",(data_3),delimiter= '\t\t',fmt = '%.8e')
    return all1

def recursion_relation2(n,x):
    list_x = []
    list_n = []
    list_nmin1 = []
    list_nplus1 = []
    data4 = []
    P_n = []
    P_nmin = []
    P_nplus = []
    l = []
    r = []
    all2 = []
    for j in x:
        data4.extend([j,n,(n-1),(n+1),legendre_m(n,j),legendre_m(n-1,j),legendre_m(n+1,j)])
        list_x.append(j)
        list_n.append(n)
        list_nmin1.append(n-1)
        list_nplus1.append(n+1)
        P_n.append(legendre_m(n,j))
        P_nmin.append(legendre_m(n-1,j))
        P_nplus.append(legendre_m(n+1,j))
        r2_l = (2*n+1)*j*legendre_m(n,j)
        r2_r = (((n +1)*legendre_m(n+1,j)) + (n*legendre_m(n-1,j)))
        l.append(r2_l)
        r.append(r2_r)
    all2.extend([list_x,list_n,list_nmin1,list_nplus1,P_n,P_nmin,P_nplus,l,r])
    data_4 = np.reshape(data4,(99,7))
    file4 = np.savetxt("/home/ishmeet/BSc physics/MP 2nd/Lab_Codes/leg03.txt",(data_4),delimiter= '\t\t',fmt = '%.8e')
    return all2

def recursion_relation3(n,x):
    list_x = []
    list_n = []
    list_nmin1 = []
    list_nplus1 = []
    data5 = []
    P_n = []
    P_nmin = []
    P_nplus = []
    l = []
    r = []
    all3 = []
    for j in x:
        r3_l = n*legendre_m(n,j)
        r3_r = ((2*n - 1)*j*legendre_m(n-1,j)) - ((n - 1)*legendre_m(n-2,j))
        data5.extend([j,n,(n-1),(n+1),legendre_m(n,j),legendre_m(n-1,j),legendre_m(n+1,j)])
        list_x.append(j)
        list_n.append(n)
        list_nmin1.append(n-1)
        list_nplus1.append(n+1)
        P_n.append(legendre_m(n,j))
        P_nmin.append(legendre_m(n-1,j))
        P_nplus.append(legendre_m(n+1,j))
        r2_l = (2*n+1)*j*legendre_m(n,j)
        r2_r = (((n +1)*legendre_m(n+1,j)) + (n*legendre_m(n-1,j)))
        l.append(r3_l)
        r.append(r3_r)
    all3.extend([list_x,list_n,list_nmin1,list_nplus1,P_n,P_nmin,P_nplus,l,r])
    data_5 = np.reshape(data5,(99,7))
    file4 = np.savetxt("/home/ishmeet/BSc physics/MP 2nd/Lab_Codes/leg04.txt",(data_5),delimiter= '\t\t',fmt = '%.8e')
    return all3

def table(x,p_o,p_1,p_2,der_po,der_p1,der_p2,der_p3,all1,all2,all3):
    table1 = PrettyTable(['X','P_o','P_1','P_2'])
    table2 = PrettyTable(['X',"P_o'","P_1'","P_2'","P_3'"])
    table3 = PrettyTable(['X','n','n-1','P_n',"P_n'","P_n-1'",'LHS','RHS'])
    table4 = PrettyTable(['X','n','n-1','n+1','P_n','P_n-1','P_n+1','LHS','RHS'])
    table5 = PrettyTable(['X','n','n-1','n+1','P_n','P_n-1','P_n+1','LHS','RHS'])
    for i in range(len(x)):
        table1.add_row([x[i],p_o[i],p_1[i],p_2[i]])
        table2.add_row([x[i],der_po[i],der_p1[i],der_p2[i],der_p3[i]])
        table3.add_row([all1[0][i],all1[1][i],all1[2][i],all1[3][i],all1[4][i],all1[5][i],all1[6][i],all1[7][i]])
        table4.add_row([all2[0][i],all2[1][i],all2[2][i],all2[3][i],all2[4][i],all2[5][i],all2[6][i],all2[7][i],all2[8][i]])
        table5.add_row([all3[0][i],all3[1][i],all3[2][i],all3[3][i],all3[4][i],all3[5][i],all3[6][i],all3[7][i],all3[8][i]])
    print(table1)
    print(table2)
    print('Recursion relation 1: \n',table3)
    print('Recursion relation 2: \n',table4)
    print('Recursion relation 3: \n',table5)

def orthogonality(m,n1,x,kron_del):
    I = 0
    inte = []
    inte1 = []
    main = []
    for j in x:
        I = I + legendre_m(m,j)*legendre_m(n1,j)*0.02
    inte.append(I)
    r = (2/(2*n1 + 1)) * kron_del
    inte1.append(r)
    main.extend([inte,inte1])
    return main

    
def graph(x,p_o,p_1,p_2,der_po,der_p1,der_p2,der_p3,func1,func2,func3):
    fig, ax = plt.subplots()
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    ax.plot(x,p_o,label = '$P_o$')
    ax.plot(x,p_1,label = '$P_1$')
    ax.plot(x,p_2,label = '$P_2$')
    ax.set(title = 'Legendre Polynomial',xlabel = 'Values of x',ylabel = 'Value of polynomial P(x)')
    ax1.plot(x,der_po,label = 'Dervatiove of $P_o$')
    ax1.plot(x,der_p1,label = 'Derivative of $P_1$')
    ax1.plot(x,der_p2,label = 'Derivative of $P_2$')
    ax1.plot(x,der_p3,label = 'Derivative of $P_3$')
    ax1.set(title = 'Derivative of Legendre Polynomial',xlabel = 'Values of x',ylabel = 'Value of derivative of polynomial P(x)')
    ax2.plot(x,func1,label = '$P_o$ via inbuilt function')
    ax2.plot(x,func2,label = '$P_1$ via inbuilt function')
    ax2.plot(x,func3,label = '$P_2$ via inbuilt function')
    ax2.set(title = 'Legendre Polynomial(Inbuilt)',xlabel = 'Values of x',ylabel = 'Value of polynomial (Inbuilt) P(x)')
    ax.grid(ls = '--')
    ax1.grid(ls = '--')
    ax2.grid(ls = '--')
    ax.legend()
    ax1.legend()
    ax2.legend()
    plt.legend()
    plt.show()

if __name__ == "__main__":
    num = float(input("Enter the number to calculate its gamma function or (num - 1) factorial: \n"))
    x = np.arange(-0.98,1,0.02)
    p_o = []
    p_1 = []
    p_2 = []
    data1 = []
    der_po = []
    der_p1 = []
    der_p2 = []
    der_p3 = []
    data2 = []
    oi1 = [] ; oi2 = []
    for j in x:
        data1.append(j)
        data2.append(j)
        for n in range(4):
            data1.append(legendre_m(n,j))
            data2.append(derivative(j,n))
            if n == 0:
                p_o.append(legendre_m(n,j))
            elif n == 1:
                p_1.append(legendre_m(n,j))
            elif n == 2:
                p_2.append(legendre_m(n,j))
            if n == 0:
                der_po.append(derivative(j,n))    
            elif n == 1:
                der_p1.append(derivative(j,n))
            elif n == 2:
                der_p2.append(derivative(j,n))
            elif n == 3:
                der_p3.append(derivative(j,n))
    data_1 = np.reshape(data1,(99,5))
    data_2 = np.reshape(data2,(99,5))
    file1 = np.savetxt("/home/ishmeet/BSc physics/MP 2nd/Lab_Codes/leg00.txt",(data_1),delimiter= '\t\t')
    file2 = np.savetxt("/home/ishmeet/BSc physics/MP 2nd/Lab_Codes/leg01.txt",(data_2),delimiter= '\t\t')
    func1 = legendre_inbuilt(x,0)
    func2 = legendre_inbuilt(x,1)
    func3 = legendre_inbuilt(x,2)
    all1 = recursion_relation1(2,x)
    all2 = recursion_relation2(2,x)
    all3 = recursion_relation3(3,x)
    for m in range(0,3):
        for n1 in range(0,3):
            kron_del = kron(m,n1)
            o = orthogonality(m,n1,x,kron_del)
            oi1.append(o[0][0])
            oi2.append(o[1][0])
    table(x,p_o,p_1,p_2,der_po,der_p1,der_p2,der_p3,all1,all2,all3)
    print('\nThe value of gamma funtion of',num,'and value of',num -1,'factorial is:',gamma(num),'\n')
    print('LHS: \n',np.reshape(oi1,(3,3)))
    print('RHS: \n',np.reshape(oi2,(3,3)))
    graph(x,p_o,p_1,p_2,der_po,der_p1,der_p2,der_p3,func1,func2,func3)        
