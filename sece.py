from hyint import hyint
from vector import vector as array
from numpy import pi, sin, cos, sign

a  = 5              # m/s^2
m  = 0.0018         # kg
D  = 0.0426
k  = 2076.5
ga = 0.00107

Cp  = 23e-9

t0 = 0.0
freq = 172.7
t1 = 50/freq
dt = 1.0 / freq / 100

F  = m * a
Upeak = F / ga
R  = D / ga**2
L  = m / ga**2
C  = ga**2 / k

def Uq(t):
    return Upeak*sin(2*pi*freq * t)

def f(t, x, y):
    Iq, UC, U = x
    dIq = (Uq(t) - R*Iq - UC - U) / L
    dUC = Iq / C
    dU  = Iq / Cp
    return array([dIq, dUC, dU])

def CHARGE(t, x, y):
    Iq, UC, U = x
    E, Umax = y
    E_new = E + Cp/2 * U**2
    return array([0.0, UC, 0.0]), array([E_new, U])

def ev_harvest(t, x):
    Iq, UC, U = x
    iCp = Iq * sign(U)
    return iCp

graph = {CHARGE : {ev_harvest : CHARGE}}
x0 = array([0, 0, 0])
y0 = array([0, 0])
t, x, y = hyint(f, x0, t0, t1, dt, graph, CHARGE, 1e-15, y0)
Iq, UC, U = zip(*x)
E, Umax = zip(*y)

from pylab import subplot, plot, ylabel, gcf, setp, show

ax1 = subplot(411)
plot(t, Iq)
ylabel('Iq')

subplot(412, sharex=ax1)
plot(t, UC)
ylabel('UC')

subplot(413, sharex=ax1)
plot(t, U)
ylabel('U')

subplot(414, sharex=ax1)
plot(t, E)
ylabel('E')

for ax in gcf().axes:
    ax.grid(True)
    setp(ax.get_yaxis().label, rotation='horizontal')
show()

Pmean = Cp/2 * Umax[-1]**2 * 2*freq
print 'Mittlere Leistung  Pmean =', Pmean
Pmax = Upeak**2 / (8*R)
print 'Maximale Leistung  Pmax  =', Pmax
print 'Verhaeltnis Pmean / Pmax =', Pmean / Pmax

def p(*args):
    return 1.0 / sum(1.0/z for z in args)

j = 1j
omega = 2*pi*freq
s = j*omega

Z = p(1/(s*Cp), R + s*L + 1/(s*C))
U0 = Upeak * 1/(s*Cp) / (R + s*L + 1/(s*C) + 1/(s*Cp))
P_max = abs(U0)**2 / (8*Z.real)
print 'Maximale Leistung (komplexe Rechnung): ', P_max

U_ =   Z.conjugate() / (Z + Z.conjugate()) * U0
