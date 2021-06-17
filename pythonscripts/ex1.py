# Import the required modules
import numpy as num
from scipy.integrate import odeint
import matplotlib.pyplot as pyplot

# define fixed parameters
n = 2       # Hill coeff
b = 1       # max exp rate
K = 100     # act threshold
g = 0.001   # degrad rate

def diff_eq(x, t):
    """Function returning differential expression"""
    dxdt = ((b * x**n) / (K**n + x**n)) - g * x
    return dxdt

#
t = num.linspace(0, 10000, 1000)
x0 = [10,10.1,10.2,10.5,11]
#x0 = [5, 10, 50, 100, 1000]

pyplot.title('Direct positive feedback model')
pyplot.xlabel('Time')
pyplot.ylabel('x(t)')
#pyplot.legend(["x = %.0f" %x])

for x in x0:
    p0 = odeint(diff_eq, x, t)

    # get threshold value
    len_p0 = len(p0)
    print("The function for x =", x,"converges towards a value of", round(float(p0[len_p0-1]),4))
    p1 = round(float(p0 [len_p0 - 1]), 3)

    # plot results using pyplot
    pyplot.plot(t, p0, label=["x = %.1f" %x])
pyplot.legend()
pyplot.show()
