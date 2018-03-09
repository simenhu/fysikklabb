import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import csv
import pandas as pd

# Input a tab seperated list of (t, x, y) values.
# Produce coefficients corresponding to a function of degree 15.
def iptrack(filename):
    data = np.loadtxt(filename, skiprows=2)
    polyfit = np.polyfit(data[:, 1], data[:, 2], 15)
    return polyfit

# Input coefficients from iptrack.
# Return a list of [y-coordinates, velocities, accelerations, alphas, radius]
def trvalues(p, x):
    y = np.polyval(p, x)
    dp = np.polyder(p)
    dydx = np.polyval(dp, x)
    ddp = np.polyder(dp)
    d2ydx2 = np.polyval(ddp, x)
    alpha = np.arctan(-dydx)
    R = (1.0 + dydx ** 2) ** 1.5 / d2ydx2
    return [y, dydx, d2ydx2, alpha, R]

data = pd.read_csv('M1Data.txt', delimiter='\t', skiprows=1)
data = data.sort_values('t')

# Constants
m = 0.0027  # Mass of ball
g = 9.82  # Constant of gravitation
c = 2/3
start_x = 0.0 # The start position in x-axis
end_x = 2.0 # The end position in x-axis
num_points = 10000 # Number of points
h = (data['t'].iloc[-1]-data['t'].iloc[0])/num_points # Stepsize

# Find the coefficients of a polynomial of degree 15 from the (t, x, y) values.
coefficients = iptrack("data_tabs.txt")

gca = plt.gca()
gca.set_ylim([-0.01, 0.05])
gca.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.005))
gca.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.05))

# This is a implemented eulers method
s_now = 0
v_now = 0
x_now = 0
y_now = data['y'].iloc[0]
s = []                     # List of positions traveled in x direction
v = []                     # List of velocity values
x = []                     # List of x values
y = []                     # List of y(x) values
a = []                     # List of accelerations
alpha_list = []            # List of alpha values (angles)
n = []                     # List of normal force values

for i in range(num_points):
    if (x_now > 1.21): break
    alpha = trvalues(coefficients, x_now)[3]
    s_next = s_now + v_now*h # h is the time derivative
    a_next = ((g*np.sin(alpha))/(1 + c))
    v_next = v_now + a_next*h
    x_next = x_now + (s_next - s_now)*np.cos(alpha)
    y_next = y_now - (s_next - s_now)*np.sin(alpha)

    n.append(np.cos(alpha) * m * g)
    alpha_list.append(alpha)
    a.append(a_next)
    x.append(x_next)
    y.append(y_next)
    s.append(s_next)
    v.append(v_next)

    x_now = x_next
    y_now = y_next
    s_now = s_next
    v_now = v_next

plt.subplot(231)
plt.plot(x, s)
plt.ylabel("posisjon x [m/s]")
plt.xlabel("tid t [s]")
data_cut = data[data['x']<=1.2]
y_interpolant = trvalues(coefficients, data_cut['x'])[0]

plt.subplot(232)
plt.plot(x, v)
plt.ylabel("hastighet v [m/s]")
plt.xlabel("tid t [s]")
plt.legend()

plt.subplot(233)
plt.plot(s, a)
plt.ylabel("akselerasjon a [m/s^2]")
plt.xlabel("posisjon x [m]")
plt.legend()
plt.plot()

plt.subplot(234)
plt.plot(x, n)
plt.ylabel("normalkraft N [mN]")
plt.xlabel("posisjon x [m]")
plt.legend()

plt.subplot(235)
f = np.subtract(0.03*(9.82*np.sin(alpha_list)), a)
plt.plot(s, f)
plt.ylabel("friksjonskraft f [N]")
plt.xlabel("posisjon x [m]")
plt.legend()
plt.plot()

plt.show()
