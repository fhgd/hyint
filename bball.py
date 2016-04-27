from hyint import hyint
from vector import vector as array

def f(t, x, y):
    """vector field from the ODE system"""
    h, v = x
    dv = -9.81
    dh = v
    return array([dh, dv])

def ev_bottom(t, x, y):
    """event for detecting if the ball hits the ground"""
    h, v = x
    return h

def RISE(t, x, y):
    """the ball touched the ground and now bounces back"""
    h, v = x
    v_new = -v * 0.8
    return array([h, v_new]), y

def FALL(t, x, y):
    """do nothing"""
    return x, y

def ev_top(t, x, y):
    """event for detecting if the ball reaches the peak"""
    h, v = x
    return v


graph = {                       # FSM graph
    FALL : {ev_bottom : RISE},
    RISE : {ev_top    : FALL},
}
z0 = FALL                       # Init FSM state
x0 = array([5, 0])              # Init values for h, v
y0 = None                       # no discrete states

t0 = 0.0
t1 = 8.0
dt = 0.02
eps = 1e-3


if __name__ == '__main__':
    # Run the simulation
    t, x, y = hyint(f, x0, t0, t1, dt, graph, z0, eps, y0)

    # Transpose the results
    h, v = zip(*x)

    #Plot the results
    from pylab import *

    inch = 2.54
    figure(figsize=(16/inch, 8/inch))
    subplots_adjust(bottom=0.2)
    plot(t, h)
    ylabel('h / m')
    xlabel('t / s')

    savefig('bball.png')
    show()
