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


dataFile = "M2Data.txt"
data = pd.read_csv('%s' % dataFile, delimiter='\t', skiprows=1)
data = data.sort_values('t')

dataFile2 = "M3Velocity.txt"
data2 = pd.read_csv('%s' % dataFile, delimiter='\t', skiprows=1)
data2 = data2.sort_values('t')

# Constants
m = 0.0027  # Mass of ball
g = 9.82  # Constant of gravitation
c = 2 / 3
start_x = 0.0  # The start position in x-axis
end_x = 2.0  # The end position in x-axis
num_points = 100000  # Number of points
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
x_now = data['x'].iloc[0] # Simen: Startposisjonen i x retning
position = []  # List of positions traveled in x direction
velocity = []  # List of velocity values
time = []  # List of x values
y = []  # List of y(x) values
x = [] # List of x values
acceleration = []  # List of accelerations
alpha_list = []  # List of alpha values (angles)
n = []  # List of normal force values

print(data)

for i in range(num_points):
    alpha = trvalues(coefficients, x_now)[3] # Burde være x_now ikke time_now
    position_next = position_now + velocity_now * h  # h is the time derivative
    acceleration_next = ((g * np.sin(alpha)) / (1 + c)) # Burde kalles acceleration now
    velocity_next = velocity_now + acceleration_next * h
    time_next = time_now + h
    y_now = y_now - (position_next - position_now) * np.sin(alpha) # The next y value to use
    x_now = x_now + (position_next - position_now)*np.cos(alpha) # The next x value to use

    # Add calculated values to lists
    n.append(np.cos(alpha) * m * g)
    alpha_list.append(alpha)
    acceleration.append(acceleration_next)
    time.append(time_next)
    y.append(y_now)
    position.append(position_next)
    velocity.append(velocity_next)

    time_now = time_next
    position_now = position_next
    velocity_now = velocity_next

plt.style.use("bmh")



def plotXandYExperimental():
    plt.subplot(211)
    plt.plot(data['t'], data['x'], label="Eksperimentell")

    plt.legend()
def plotNumericalPosition():
    global gca
    plt.plot(time, position, label="Numerisk")
    plt.ylabel("posisjon x [m/s]")
    plt.xlabel("tid t [s]")
    plt.legend()
    gca = plt.gca()
    gca.set_xlim([0, 1.5])

def plotNumericalAcceleration():
    global gca
    plt.subplot(211)
    plt.plot(time, acceleration, label ="Numerisk")
    plt.ylabel("akselerasjon a [m/s^2]")
    plt.legend()
    plt.xlabel("tid t [s]")
    gca = plt.gca()
    gca.set_xlim([0, 1.5])


def plotFandN():
    global gca
    plt.subplot(211)
    plt.plot(position, n)
    plt.ylabel("normalkraft N [mN]")
    plt.xlabel("posisjon x [m]")

    gca = plt.gca()
    gca.set_xlim([0, 1.371])

    plt.subplot(212)
    f = np.subtract(np.multiply(m * g, np.sin(alpha_list)), np.multiply(m, acceleration))
    plt.plot(position, f)

    plt.ylabel("friksjonskraft f [mN]")
    plt.xlabel("posisjon x [m]")
    gca = plt.gca()
    gca.set_xlim([0, 1.371])


def plotExperimentalVelocity():
    global gca
    plt.subplot(212)

    # Data2['x'] is actually Data2['v']
    plt.plot(data2['t'], data2['x'], label="Eksperimentell")
    plt.ylabel("Hastighet v [m/s]")
    plt.legend()
    plt.xlabel("tid t [s]")

def plotNumericalVelocity():
    global gca
    plt.subplot(212)
    plt.plot(time, velocity, label="Numerisk")
    plt.legend()
    plt.ylabel("hastighet v [m/s]")
    plt.xlabel("tid t [s]")
    gca = plt.gca()

    gca.set_xlim([0, 1.5])



def plotF():

    plt.subplot(211)
    f = np.subtract(np.multiply(m * g, np.sin(alpha_list)), np.multiply(m, acceleration))
    plt.plot(position, f, label="Numerisk")
    plt.legend()
    plt.ylabel("friksjonskraft f [mN]")
    plt.xlabel("posisjon x [m]")
    gca = plt.gca()
    gca.set_xlim([0, 1.371])


plotNumericalAcceleration()

plt.show()
