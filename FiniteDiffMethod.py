from math import *
from numpy import linspace
import matplotlib.pyplot as plt
def p_x(x, e):
    return ((sin(x))**2+cos(x))/e

def q_x(x, e):
    return (x**2+x)/e

def f_x(x, e):
    return (-x**3+x**2-x)/e

def Solve_3diag(a, b, c, d):
    n = len(b)
    a.insert(0, 0)
    P = [0]*(n-1)
    Q = [0]*n
    P[0] = -c[0]/b[0]
    Q[0] = d[0]/b[0]
    for i in range(1, n-1):
        P[i] = c[i]/(-b[i] - a[i]*P[i-1])
        Q[i] = (a[i]*Q[i-1] - d[i])/(-b[i] - a[i]*P[i-1])
    x = [0]*n
    x[n-1] = (a[n-1]*Q[n-2] - d[n-1])/(-b[n-1] - a[n-1]*P[n-2])
    for i in range(n-2, -1, -1):
        x[i] = P[i]*x[i+1] + Q[i]
    return x

def Fin_Diff_Method(x, h, e, a, b):
    n = round((b - a)/h)
    c1, c2, c3, c4 = [], [], [], []
    c2.append(-1/h)
    c3.append(1/h)
    c4.append(-3)
    for i in range(1, n):
        p = p_x(x[i], e)
        q = q_x(x[i], e)
        f = f_x(x[i], e)
        c1.append(1/(h**2)-p/(2*h))
        c2.append(-2/(h**2)+q)
        c3.append(1/(h**2)+p/(2*h))
        c4.append(f)
    c1.append(-1/h)
    c2.append(1/h)
    c4.append(6)
    y = Solve_3diag(c1, c2, c3, c4)
    return y

def Delta(y1, y2):
    n = len(y1)
    dif = []
    for i in range(n):
        dif.append(abs(y1[i]-y2[2*i]))
    return max(dif)

def Error1(h, e, a, b): #Погрешность
    n = round((b - a)/h) #Целый шаг
    x = linspace(a, b, n+1)
    y11 = Fin_Diff_Method(x, h, e, a, b)
    h /= 2 #Половинный шаг
    n = round((b - a)/h)
    x = linspace(a, b, n+1)
    y12 = Fin_Diff_Method(x, h, e, a, b)
    h /= 2 #Четверть шага
    n = round((b - a)/h)
    x = linspace(a, b, n+1)
    y13 = Fin_Diff_Method(x, h, e, a, b)
    h /= 2 #Восьмая шага
    n = round((b - a)/h)
    x = linspace(a, b, n+1)
    y14 = Fin_Diff_Method(x, h, e, a, b)
    delta1 = Delta(y11, y12)
    delta2 = Delta(y12, y13)
    delta3 = Delta(y13, y14)
    print('e = ', e)
    print(delta1, delta2, delta3, sep='\n')
    print()

def Demonstrate(h, a=0, b=1):
    n = round((b - a)/h)
    x = linspace(a, b, n+1)
    y1 = Fin_Diff_Method(x, h, 1, a, b)
    y2 = Fin_Diff_Method(x, h, 0.1, a, b)
    y3 = Fin_Diff_Method(x, h, 0.01, a, b)
    y4 = Fin_Diff_Method(x, h, 0.001, a, b)
    Error1(h, 1, a, b)
    Error1(h, 0.1, a, b)
    Error1(h, 0.01, a, b)
    Error1(h, 0.001, a, b)
    plt.subplot(2, 3, 1)
    plt.plot(x, y1, color='red')
    plt.plot(x, y2, color='blue')
    plt.plot(x, y3, color='green')
    plt.plot(x, y4, color='orange')
    plt.legend(['eps = 1', 'eps = 0.1', 'eps = 0.01', 'eps = 0.001'])
    plt.subplot(2, 3, 2)
    plt.plot(x, y1)
    plt.title('eps = 1')
    plt.subplot(2, 3, 3)
    plt.plot(x, y2)
    plt.title('eps = 0.1')
    plt.subplot(2, 3, 4)
    plt.plot(x, y3)
    plt.title('eps = 0.01')
    plt.subplot(2, 3, 5)
    plt.plot(x, y4)
    plt.title('eps = 0.001')
    plt.show()
