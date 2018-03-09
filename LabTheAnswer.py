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
    print(polyfit)
    return polyfit


# Input coefficients from iptrack.
# Return a list of [y-coordinates, velocities, accelerations, alphas, radius']
def trvalues(p, x):
    y = np.polyval(p, x)
    dp = np.polyder(p)
    dydx = np.polyval(dp, x)
    ddp = np.polyder(dp)
    d2ydx2 = np.polyval(ddp, x)
    alpha = np.arctan(-dydx)
    R = (1.0 + dydx ** 2) ** 1.5 / d2ydx2
    return [y, dydx, d2ydx2, alpha, R]


data = pd.read_csv('data_tabs.txt', delimiter='\t', skiprows=1)
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

# Domain space from 0 to 1.20, 10 000 points for precision.
# x = np.linspace(start_x, end_x, num_points)

# Values are a list of [y(x), dy/dx, d^2y/d^2x, alphas, radius']
#values = trvalues(coefficients, x)

# Acceleration values of the rolling ball.
#acceleration_values = values[2]

# alpha_values are the angles with respect to the horizontal plane.
#alpha_values = values[3]

# The alpha_values transformed by the function sin(x).
#sinus_alpha_values = np.sin(alpha_values)

#plt.subplot(311)
# This will plot the friction force with respect to x.
#plt.plot(x, m * g * sinus_alpha_values - m * acceleration_values)

gca = plt.gca()
gca.set_ylim([-0.01, 0.05])
gca.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.005))
gca.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.05))


# This is a implemented eulers method
s_now = 0
v_now = 0
x_now = 0
y_now = data['y'].iloc[0]
s = []
v = []
x = []
y = []

for i in range(num_points):
    alpha = trvalues(coefficients, x_now)[3]
    s_next = s_now + v_now*h # h is the time derivative
    v_next = v_now + ((g*np.sin(alpha))/(1 + c))*h
    x_next = x_now + (s_next - s_now)*np.cos(alpha)
    y_next = y_now - (s_next - s_now)*np.sin(alpha)
    x.append(x_next)
    y.append(y_next)
    s.append(s_next)
    v.append(v_next)
    x_now = x_next
    y_now = y_next
    s_now = s_next
    v_now = v_next


plt.subplot(221)
plt.plot(x, s, label='S(x)')
plt.legend()
plt.subplot(222)
plt.plot(x, v, label='v(x)')
plt.legend()
plt.subplot(223)
plt.plot(x, y)
plt.plot()

plt.show()


print(v[:100])