import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from sympy import *
from sympy.solvers.ode.systems import dsolve_system

t = symbols('t')
a1 = Function('a1')(t)
a2 = Function('a2')(t)
q = Matrix([[a1],[a2]])
dq = diff(q,t)
eqns = Matrix([[0],[0]])
m1 = 0.5
m2 = 0.15
l1 = 0.1
l2 = 0.1
g = 9.8

x1 = l1*sin(a1)
y1 = -l1*cos(a1)
x2 = l1*sin(a1) + l2*sin(a2)
y2 = -l1*cos(a1) - l2*cos(a2)

T = 0.5*(m1*((diff(x1,t))**2 + (diff(y1,t))**2) + m2*((diff(x2,t))**2 + (diff(y2,t))**2))
V = m1*g*y1 + m2*g*y2
L = T-V

for i in range(0,2):
    eqns[i] = diff(diff(L,dq[i]),t) - diff(L,q[i])

print(eqns[1])
print(eqns[0])
# ics = {a1(0):np.pi/2,a2(0):np.pi/3}
funcs = [a1,a2]
# sol = dsolve_system(eqns,funcs=funcs,t=t)
print(classify_ode(eqns[0]))


