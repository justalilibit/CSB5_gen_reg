# Import the required modules
import numpy as num
from scipy.integrate import odeint
import matplotlib.pyplot as pyplot

#define given parameters
a_k = 0.875
a_s = 0.0004
b_k = 7.5
b_s = 0.06
g = 0.001
d = 0.0001
T_k = 25000
T_s = 20
k_k = 5000
k_s = 833
n = 2
p = 5

def diff_eq(p0, t, a_k, a_s, b_s):
    """Function returning the differential equation"""
    S = p0[1]
    K = p0[0]

    dKdt = (a_k+ ((b_k * K**n) / (k_k**n + K**n))) - ((g*K) / (1+(K/T_k) + (S/T_s))) - d*K
    dSdt = (a_s+ (b_s / (1 + (K/k_s)**p)))-((g*S) / (1+(K/T_k) + (S/T_s))) - d*S
    return [dKdt, dSdt]

# set initial state parameters
p0 = [1, 0.1]
t = num.linspace(0, 1000000, 1000)

# solve model
p1 = odeint(diff_eq, p0, t, args=(a_k, a_s, b_s))
dK = p1[:,0]
dS = p1[:,1]

# plot results using pyplot
pyplot.title(f'Model of genetic competence \n with a_k = {a_k}, a_s = {a_s}, b_s= {b_s}')
pyplot.xlabel('Time in s')
pyplot.ylabel('Concentration in nM')
pyplot.plot(t, dK)
pyplot.plot(t, dS)
pyplot.legend(["ComK", "ComS"])
pyplot.show()

# get values for K and S
len_dK = len(dK)
len_dS = len(dS)
print("[ComK] converges towards a value of ", round(float(dK[len_dK - 1]), 3))
print("[ComS] converges towards a value of ", round(float(dS[len_dS - 1]), 3))
