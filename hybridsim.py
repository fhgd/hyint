#!/usr/bin/env python
# -*- coding: utf-8 -*-

def fsolve(g, x0, x1, eps):
    """Zero finding with the bisect methode.

    Where x0 < x1 and g(x0) * g(x1) <= 0. The algorithm is from

    Hans Petter Langtangen: A Primer on Scientific Programming with Python.
    Page 156. http://books.google.com/books?id=cVof07z_rA4C
    """

    assert x0 < x1, 'x0 >= x1'
    g0 = g(x0)
    assert g0*g(x1) <= 0, 'g does not change sign in [x0, x1]'
    while x1 - x0 > eps:
        xm = (x0 + x1) / 2.0
        gm = g(xm)
        if g0*gm <= 0:
            x1 = xm     # zero is in left half
        else:
            x0 = xm     # zero is in right half
            g0 = gm
    return x0, x1

def hyint(f, x0, t0, t1, dt, graph, z0, eps, y0):
    """Integration of a hybrid system (with a first order ode system)"""

    # start values
    t = [t0]            # time
    x = [x0]            # time continuous states
    y = [y0]            # time discrete states
    z = z0              # state from the finite state machine (FSM)

    # integrate until time is over
    while t[-1] <= t1:
        # test all event functions
        events = [ev for ev in graph[z].iterkeys() if ev(t[-1], x[-1]) < 0]
        if events:
            print t[-1], events
            # at least one event is detected, so find the first event
            k_min = 1.0
            event = None
            for ev in events:
                ev_local = lambda k: ev(t[-2] + k*dt, x[-2] + (x[-1] - x[-2])*k)
                k0, k1 = fsolve(ev_local, 0.0, 1.0, eps)
                # Use k1 which terminates the actual process
                if k1 <= k_min:
                    k_min = k1
                    event = ev
            assert event, 'No zero was found in [0, 1]'
            # correct the last integration step, which was too far
            x[-1] = x[-2] + (x[-1] - x[-2])*k_min
            t[-1] = t[-2] + k_min*dt
            # transition and action of the FSM until all events are off
            while 1:
                z = graph[z][event]
                print z.__name__
                x_new, y_new = z(t[-1], x[-1], y[-1])
                x.append(x_new)
                y.append(y_new)
                t.append(t[-1])
                # test all event functions
                events = [ev for ev in graph[z].iterkeys() if ev(t[-1], x[-1]) < 0]
                if events:
                    event = events[0]
                    assert len(events) == 1, \
                        'More the one active event after a FSM transaction.'
                else:
                    break
        # integrate with an euler step
        x.append(x[-1] + f(t[-1], x[-1]) * dt)
        t.append(t[-1] + dt)
        # constant extension of the time discrete vaules
        y.append(y[-1])
    return t, x, y

if __name__ == '__main__':

    Upeak = 2.5
    freq = 750.
    t1 = 5./freq

    R  = 5.
    L  = 800e-6
    C  = 10e-6
    deltaU = 0.1

    from numpy import pi, sin, cos
    from vector import vector

    def Uq(t):
        return Upeak*sin(2*pi*freq * t)

    def iUq(t):
        return Upeak/(2*pi*freq) * (1 - cos(2*pi*freq * t))

    def dUq(t):
        return 2*pi*freq * Upeak*cos(2*pi*freq * t)

    def d2Uq(t):
        return -(2*pi*freq)**2 * Upeak*sin(2*pi*freq * t)

    def Uref(t):
        return 0.5 * (Uq(t) - L/R*dUq(t) - 1/(R*C)*iUq(t))

    def f(x, t):
        I, UC, E = x
        dI  = (Uq(t) - R*I - UC - Uref(t)) / L
        dUC = I / C
        dE  = Uref(t) * I
        return vector([dI, dUC, dE])

    from scipy.integrate import odeint
    from numpy import linspace

    t = linspace(0, t1, 500)
    z = odeint(f, [0, 0, 0], t)
    I, UC, E = z.T
    Pmean = E[-1] / t1
    print 'Mittlere Leistung  P =', Pmean

    #~ def HARVEST(t, x, E):
        #~ UC, IM, UK = x
        #~ E_new = E + C/2*(UC**2 - UC_ref(t, x)**2)
        #~ return vector([UC_ref(t, x), IM, UK]), E_new

    #~ def ev_too_low(t, x):
        #~ UC, IM, UK = x
        #~ return UC - (UC_ref(t, x) - deltaUC)

    #~ def ev_too_high(t, x):
        #~ UC, IM, UK = x
        #~ return (UC_ref(t, x) + deltaUC) - UC

    #~ graph = {HARVEST : {ev_too_low : HARVEST, ev_too_high : HARVEST}}
    #~ x0 = vector([0, 0, 0])
    #~ t, x, y = hyint(f, x0, 0.0, 3*math.pi, 0.001, graph, HARVEST, 1e-15, 0.0)

    #~ UC, IM, UK = [vector(_) for _ in zip(*x)]
    #~ ULC = vector(map(UF, t)) - RD*IM - UC
