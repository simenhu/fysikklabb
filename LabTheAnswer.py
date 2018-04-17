import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

SMALL_SIZE = 9
MEDIUM_SIZE = 30
BIGGER_SIZE = 30

plt.rc('font', size=BIGGER_SIZE)            # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)       # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)       # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)      # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)      # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)      # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)     # fontsize of the figure title
left = 0.09                                 # the left side of the subplots of the figure
right = 0.98                                # the right side of the subplots of the figure
bottom = 0.09                               # the bottom of the subplots of the figure
top = 0.97                                  # the top of the subplots of the figure
wspace = 0.21                               # width between subplots
hspace = 0.31                               # height between subplots

plt.subplots_adjust(left, bottom, right, top, wspace, hspace)


def iptrack(filename):
    """
    Input a tab seperated list of (t, x, y) values.
    Produce coefficients corresponding to a function of degree 15.
    :param filename: name of file
    :return:
    """
    data = np.loadtxt(filename, skiprows=2)
    polyfit = np.polyfit(data[:, 1], data[:, 2], 15)
    return polyfit


def trvalues(p, x):
    """
    Input coefficients from iptrack.
    Return a list of [y-coordinates, velocities, accelerations, alphas, radius]
    :param p:
    :param x:
    :return:
    """
    y = np.polyval(p, x)
    dp = np.polyder(p)
    dydx = np.polyval(dp, x)
    ddp = np.polyder(dp)
    d2ydx2 = np.polyval(ddp, x)
    alpha = np.arctan(-dydx)
    R = (1.0 + dydx ** 2) ** 1.5 / d2ydx2
    return [y, dydx, d2ydx2, alpha, R]


# Import data and sort it
fileName = "M2DataWithVelocity.txt"
dataFile = pd.read_csv('%s' % fileName, delimiter='\t', skiprows=1)
dataFile = dataFile.sort_values('t')

# Constants
m = 0.0027                                  # Mass of ball
g = 9.82                                    # Constant of gravitation
c = 2 / 3                                   # Torque constant for ball
start_x = 0.0                               # The start position in x-axis
end_x = 2.0                                 # The end position in x-axis
num_points = 196                            # Number of points
h = dataFile['t'].iloc[-1] / num_points     # Time step size

# Find the coefficients of a polynomial of degree 15 from the (t, x, y) values.
coefficients = iptrack(fileName)

# Constants for and implementation of Eulers method starts here
position_now = 0
velocity_now = 0                                            # Position the traveled direction
time_now = 0
y_now = dataFile['y'].iloc[0]
x_now = dataFile['x'].iloc[0]
position = []                                               # List of positions traveled in the track
velocity = []                                               # List of velocity values
time = np.linspace(0, dataFile['t'].iloc[-1], num_points)   # list of time samples to use when plotting by time
y = []                                                      # List of y(x) values
x = []                                                      # LIst of x values
acceleration = []                                           # List of accelerations
alpha_list = []                                             # List of alpha values (angles)
n = []                                                      # List of normal force values

for i in range(num_points):
    alpha = trvalues(coefficients, position_now)[3]
    position_next = position_now + velocity_now * h
    acceleration_next = ((g * np.sin(alpha)) / (1 + c))
    velocity_next = velocity_now + acceleration_next * h
    time_next = time_now + h
    y_next = y_now - (position_next - position_now) * np.sin(alpha)
    x_next = x_now + (position_next - position_now) * np.cos(alpha)

    n.append(np.cos(alpha) * m * g)
    alpha_list.append(alpha)
    acceleration.append(acceleration_next)
    y.append(y_next)
    x.append(x_next)
    position.append(position_next)
    velocity.append(velocity_next)

    y_now = y_next
    x_now = x_next
    position_now = position_next
    velocity_now = velocity_next

plt.style.use("bmh")

f = np.subtract(np.multiply(m * g, np.sin(alpha_list)), np.multiply(m, acceleration))


def plot_experimental_and_numerical_position_and_velocity():
    """
    Compares the experimental and numerical values of position and value.
    """
    plt.subplot(211)
    plt.plot(dataFile['t'], dataFile['y'], label="Eksperimentell")
    plt.ylabel("posisjon y [m]")
    plt.xlabel("tid t [s]")
    plt.legend()
    gca = plt.gca()
    gca.set_xlim([0, 1.5])

    plt.plot(time, y, label="Numerisk")
    plt.ylabel("posisjon x [m]")
    plt.xlabel("tid t [s]")
    plt.legend()
    gca = plt.gca()
    gca.set_xlim([0, 1.5])
    gca.set_ylim([0.45, 0.65])

    plt.subplot(212)
    plt.plot(dataFile['t'], dataFile['v'], label="Eksperimentell")
    plt.ylabel("hastighet v [m/s]")
    plt.xlabel("tid t [s]")
    plt.legend()
    gca = plt.gca()
    gca.set_xlim([0, 1.5])
    gca.set_ylim([0, 2])

    plt.plot(time, velocity, label="Numerisk")
    plt.ylabel("hastighet v [m/s]")
    plt.xlabel("tid t [s]")
    plt.legend()


def plot_friction_force_and_acceleration():
    """
    Plot the numerical acceleration and friction force.
    """
    ax1 = plt.subplot(211)
    plt.plot(x, acceleration, label="Numerisk")
    plt.ylabel("akselerasjon a [m/s²]")
    plt.legend()
    plt.xlabel("posisjon x [m]")
    gca = plt.gca()
    gca.set_xlim([0, 1.371])
    plt.axvline(x=0.48, color='r', linestyle='dashed')
    plt.axvline(x=1.016, color='r', linestyle='dashed')

    plt.subplot(212, sharex=ax1)

    plt.plot(x, f)
    plt.ylabel("friksjonskraft f [mN]")
    plt.xlabel("posisjon x [m]")
    plt.axvline(x=0.48, color='r', linestyle='dashed')
    plt.axvline(x=1.016, color='r', linestyle='dashed')


def plot_friction_and_normal_force():
    """
    Plot friction and normal force.
    """
    # This plots the normal force N(x).
    ax1 = plt.subplot(211)
    plt.plot(x, n)
    plt.ylabel("normalkraft N [mN]")
    plt.xlabel("posisjon x [m]")
    gca = plt.gca()
    gca.set_xlim([0, 1.371])
    gca.set_ylim([0.019, 0.027])
    plt.axvline(x=0.48, color='r', linestyle='dashed')
    plt.axvline(x=1.016, color='r', linestyle='dashed')

    # This plots the friction force f(x).
    plt.subplot(212, sharex=ax1)
    plt.plot(x, f)
    plt.ylabel("friksjonskraft f [mN]")
    plt.xlabel("posisjon x [m]")
    plt.axvline(x=0.48, color='r', linestyle='dashed')
    plt.axvline(x=1.016, color='r', linestyle='dashed')


def plot_interpolated_curve_with_experimental_curve():
    plt.subplot(211)

    x = np.linspace(0, 1.353721644E0, len(dataFile['x']))

    y = np.polyval(coefficients, x)
    plt.plot(x, y, label="Interpolant")
    plt.plot(dataFile['x'], dataFile['y'], label="Eksperimentell")
    plt.ylabel("posisjon y [m]")
    plt.xlabel("posisjon x [m]")
    plt.legend()

    diff = y - dataFile['y']
    plt.subplot(212)
    plt.plot(x, diff, label="Differanse")
    plt.ylabel("posisjon y [m]")
    plt.xlabel("posisjon x [m]")
    plt.legend()


# plot_experimental_and_numerical_position_and_velocity()
# plot_friction_and_normal_force()
# plot_friction_force_and_acceleration()
# plot_interpolated_curve_with_experimental_curve()
plt.show()
