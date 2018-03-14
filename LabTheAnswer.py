import numpy as np
import matplotlib
from matplotlib import pyplot as plt
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


dataFile = "M3Data.txt"
data = pd.read_csv('%s' % dataFile, delimiter='\t', skiprows=1)
data = data.sort_values('t')

# Constants
m = 0.0027  # Mass of ball
g = 9.82  # Constant of gravitation
c = 2 / 3
start_x = 0.0  # The start position in x-axis
end_x = 2.0  # The end position in x-axis
num_points = 10000  # Number of points
h = (data['t'].iloc[-1] - data['t'].iloc[0]) / num_points  # Stepsize

# Find the coefficients of a polynomial of degree 15 from the (t, x, y) values.
coefficients = iptrack(dataFile)

gca = plt.gca()
gca.set_ylim([-0.01, 0.05])
gca.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.005))
gca.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.05))

# This is a implemented eulers method
position_now = 0
velocity_now = 0
time_now = 0
y_now = data['y'].iloc[0]
position = []  # List of positions traveled in x direction
velocity = []  # List of velocity values
time = []  # List of x values
y = []  # List of y(x) values
acceleration = []  # List of accelerations
alpha_list = []  # List of alpha values (angles)
n = []  # List of normal force values

for i in range(num_points):
    alpha = trvalues(coefficients, time_now)[3]
    position_next = position_now + velocity_now * h  # h is the time derivative
    acceleration_next = ((g * np.sin(alpha)) / (1 + c))
    velocity_next = velocity_now + acceleration_next * h
    time_next = time_now + h
    y_next = y_now - (position_next - position_now) * np.sin(alpha)

    # Add calculated values to lists
    n.append(np.cos(alpha) * m * g)
    alpha_list.append(alpha)
    acceleration.append(acceleration_next)
    time.append(time_next)
    y.append(y_next)
    position.append(position_next)
    velocity.append(velocity_next)

    time_now = time_next
    y_now = y_next
    position_now = position_next
    velocity_now = velocity_next

plt.style.use("bmh")

plt.subplot(231)
plt.plot(time, position)
plt.ylabel("posisjon x [m/s]")
plt.xlabel("tid t [s]")
plt.title("Numerisk posisjon/tid")


plt.subplot(232)
plt.plot(time, velocity)
plt.ylabel("hastighet v [m/s]")
plt.xlabel("tid t [s]")
plt.title("Numerisk hastighet/tid")

plt.subplot(233)
plt.plot(time, acceleration)
plt.ylabel("akselerasjon a [m/s^2]")
plt.xlabel("tid t [s]")
plt.title("Numerisk akselerasjon/tid")


#plt.subplot(234)
#plt.plot(position, n)
#plt.ylabel("normalkraft N [mN]")
#plt.xlabel("posisjon x [m]")
#plt.title("Numerisk pos")
#plt.legend()

plt.subplot(234)
f = np.subtract(np.multiply(m * g, np.sin(alpha_list)), np.multiply(m, acceleration))
plt.plot(time, f)
plt.ylabel("friksjonskraft f [N]")
plt.xlabel("tid t [s]")
plt.title("Numerisk friksjonskraft/tid")


plt.show()