#!/usr/bin/env python
#
# Piezo generator model:
#
#       ,-->-- R -- L -- C ----+-->--o
#      +|  Iq           +UC-   |  I  +
#      Uq                      Cp    U
#      -|                      |     -
#       '----------------------+-----o
#
# Interface principle: Synchronous electrical charge extraction (SECE)
#
#       (1) Leave the generator in open circuit.
#       (2) When the generator voltage U has a extremum, then immediately
#           discharge the piezo capacitor Cp. (Idealy, this discharge doesn't
#           take any time.)

from hyint import hyint
from vector import vector as array
from numpy import pi, sin, cos, sign

# Piezo model parameters
#       Mechanical
a  = 5              # m/s^2
m  = 0.0018         # kg
D  = 0.0426
k  = 2076.5
ga = 0.00107
#       Electrical
Cp  = 23e-9         # F
#       Calculate the circuit parameters
F  = m * a
UqPeak = F / ga
R  = D / ga**2
L  = m / ga**2
C  = ga**2 / k

# Simulation parameters
t0 = 0.0
freq = 172.7
t1 = 10/freq
dt = 1.0 / freq / 100
eps = 1e-15

def Uq(t):
    """Source voltage"""
    return UqPeak*sin(2*pi*freq * t)

def f(t, x, y):
    """Vector field from the ODE system"""
    Iq, UC, U = x                               # time continous state vector
    dIq = (Uq(t) - R*Iq - UC - U) / L
    dUC = Iq / C
    dU  = Iq / Cp
    return array([dIq, dUC, dU])

def ev_harvest(t, x, y):
    """The harvesting event occurs when the voltage U has an extremum. This is
    equivalent with I_Cp = Iq = 0. To distinguish the positive and negative
    zero crossing, the event function becomes Iq * sign(U) <= 0."""
    Iq, UC, U = x
    return Iq * sign(U)

def CHARGE(t, x, y):
    """The SECE principle could modeled with only one FSM state CHARGE. After
    every extremum of the output voltage U is reseted to zero. Also the
    current Iq is set to zero avoiding numerical problems. Additional the
    time discrete state vector saves the energy E transfered through the
    output and every extremum of the voltage U."""
    Iq, UC, U = x
    E, Umax = y                                 # time discrete state vector
    E_new = E + Cp/2 * U**2
    return array([0.0, UC, 0.0]), array([E_new, U])

SECEgraph = {                           # FSM graph
    CHARGE : {ev_harvest : CHARGE},
}
z0 = CHARGE                             # Init FSM state
x0 = array([0, 0, 0])                   # Init values for Iq, UC, U
y0 = array([0, 0])                      # Init values for E, Umax

if __name__ == '__main__':
    # Run the simulation
    t, x, y = hyint(f, x0, t0, t1, dt, SECEgraph, z0, eps, y0)
    # Transpose the results
    Iq, UC, U = zip(*x)
    E, Umax = zip(*y)

    # If the system gets stationary, the last Umax is
    # related to the averaged power
    Pmean = Cp/2 * Umax[-1]**2 * 2*freq
    print 'Pmean =', Pmean
    # The maximum output power in the case of impedance matching
    Pmax = UqPeak**2 / (8*R)
    print 'Pmax  =', Pmax
    print 'Pmean / Pmax =', Pmean / Pmax

    # Plot the results
    from pylab import *

    inch = 2.54
    figure(figsize=(16/inch, 8/inch))
    subplots_adjust(bottom=0.2)
    plot(t, U)
    ylabel('Upiezo / V')
    xlabel('t / s')

    savefig('piezo_sece.png')
    show()
