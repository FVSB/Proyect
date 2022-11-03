import numpy as np
import matplotlib.pyplot as plt
import math


def f(x, y): return x+y**1/3
def f2(t, v): return -0.04*v - 9.8
def Solf2(t): return (294*math.e**(-t/25)) - 245


def Numeric_Methods(X0, Y0, h, f, method, range_end):
    xn = np.arange(X0, range_end, h)
    yn = [Y0]
    dic = {}
    for i in range(0, len(xn)-1):
        q = xn[i]
        a = round(method(q, yn[i], h, f), 2)
        dic[q] = a
        yn.append(a)
    return xn, yn


def Method_euler(Xn, Yn, h, f):
    yn1 = Yn + (h * f(Xn, Yn))
    return yn1


def Method_eulerMejorado(Xn, Yn, h, f):
    k1 = f(Xn, Yn)
    un1 = Yn + (h * k1)
    k2 = f(Xn+h, un1)
    yn1 = Yn + (h * (1/2)*(k1 + k2))
    return yn1


def Runge_Kutta(Xn, Yn, h, f):
    k1 = f(Xn, Yn)
    k2 = f(Xn + (1/2)*h, Yn+(1/2)*h*k1)
    k3 = f(Xn + (1/2)*h, Yn+(1/2)*h*k2)
    k4 = f(Xn+h, Yn + h*k3)
    Yn1 = Yn + (h/6 * (k1 + 2*k2 + 2*k3 + k4))
    return Yn1


def SolutionValues():
    xn = np.arange(0, 10, 0.1)
    yn = []
    for i in range(0, len(xn)):
        yn.append(round(Solf2(xn[i]), 2))
    return xn, yn


ptosEulerM = Numeric_Methods(0, 49, 0.2, f2, Method_eulerMejorado, 10)
ptosEulerM2 = Numeric_Methods(0, 49, 0.1, f2, Runge_Kutta, 10)


ptos = SolutionValues()


print(ptosEulerM)
print("------------")
# print(ptosEulerM2)
print("------------")
print(ptos)


def Plot():
    plt.plot(ptosEulerM[0], ptosEulerM[1])
    plt.plot(ptosEulerM2[0], ptosEulerM2[1])
    plt.plot(ptos[0], ptos[1])
    plt.show()


def Test():
    a = ptosEulerM[1]
    b = ptos[1]
    for i in range(len(ptos)):
        c = abs(a[i]-b[i])
        if (c > 0.9):
            print("Error")
            print([a[i]])
            print(b[i])

            break


Test()

a = {}

a["d"] = 1
print(a["d"])
