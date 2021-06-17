# import from modules
import numpy as num
from scipy.integrate import odeint
import matplotlib.pyplot as pyplot

# define fixed parameters
a2 = 0.025
b1 = 15
b2 = 0.8
d = 5E-5
g =  2.5E-7
K1 = 3000
K2 = 750
n = 2

def diff_eq(p, t, a1):
    """Function returning the differential equation with given parameters"""
    x = p[0]
    y = p[1]

    dxdt = a1 + ((b1*(x**n)/((K1**n+x**n))))- g*x*y - d*x
    dydt = a2 + ((b2*(x**n))/((K2**n) + (x**n))) - d*y

    return[dxdt, dydt]

# set initial state parameters
p0 = [0, 100]
t = num.linspace(0, 2000000, 1000)

# Examine different values for a1
a1 = [0.005, 0.007, 0.01]

for a in a1:
    p1 = odeint(diff_eq, p0, t, args=(a,))
    x = p1[:,0]
    y = p1[:,1]
    # plot results using pyplot
    pyplot.title("Activator-Repressor Model with a1 = %.3f" %a)
    pyplot.xlabel("time")
    pyplot.ylabel("concentration")
    pyplot.plot(t,x,label="activator")
    pyplot.plot(t,y, label="repressor")
    pyplot.legend()
    pyplot.show()
