#!/usr/bin/env python
# -*- coding: utf-8 -*-

def fsolve(g, x0, x1, eps):
    """Zero finding with the secant methode"""

    while abs(g(x1)) > eps:
        x0, x1 = x1, x1 - g(x1) * (x1 - x0) / (g(x1) - g(x0))
    return x1

def odeint(f, x0, t0, t1, dt, graph, z0, eps, y0):
    """Integration of a hybrid system (with a first order ode system)"""

    # start values
    t = [t0]            # time
    x = [x0]            # time continuous states
    y = [y0]            # time discrete states
    z = z0              # state from the finite state machine (FSM)

    while t[-1] < t1:

        # integrate until at least one event is detected
        while 1:
            # integrate with an euler step
            x.append(x[-1] + f(t[-1], x[-1]) * dt)
            t.append(t[-1] + dt)
            # constant extention of the time discrete vaules
            y.append(y[-1])
            # test all event functions
            events = [ev for ev in graph[z].iterkeys() if ev(t[-1], x[-1]) < 0]
            if events:
                # at least one event is detected
                break

        # find the first event
        dt_min = dt
        x_local = lambda t_: x[-2] + f(t[-2] + t_, x[-2]) * t_
        for ev in events:
            ev__x_local = lambda t_: ev(t[-2] + t_, x_local(t_))
            dt_ = fsolve(ev__x_local, 0, dt, eps)
            if dt_ <= dt_min:
                dt_min = dt_
                event = ev

        # correct the last integration step, which was too far
        x[-1] = x_local(dt_min)
        t[-1] = t[-2] + dt_min

        # transition and action of the FSM
        z = graph[z][event]
        x_new, y_new = z(t[-1], x[-1], y[-1])
        x.append(x_new)
        y.append(y_new)
        t.append(t[-1])
    return t, x, y


if __name__ == '__main__':

    R  = 0.5
    C  = 1.0
    CK = 1.0
    LM = 1.0
    RD = 2.0
    deltaUC = 0.01

    import math
    from vector import vector

    def UF(t):
        return math.sin(t) # + 0.8*math.sin(2*t)

    def UC_ref(t, x):
        UC, IM, UK = x
        return 0.5 * R/(R + RD) * UF(t)

    def f(t, x):
        UC, IM, UK = x
        dUC = (IM - UC/R) / C
        dIM = (UF(t) - RD*IM - UK - UC) / LM
        dUK = IM / CK
        return vector([dUC, dIM, dUK])

    def HARVEST(t, x, E):
        UC, IM, UK = x
        E_new = E + C/2*(UC**2 - UC_ref(t, x)**2)
        return vector([UC_ref(t, x), IM, UK]), E_new

    def ev_too_low(t, x):
        UC, IM, UK = x
        return UC - (UC_ref(t, x) - deltaUC)

    def ev_too_high(t, x):
        UC, IM, UK = x
        return (UC_ref(t, x) + deltaUC) - UC

    graph = {HARVEST : {ev_too_low : HARVEST, ev_too_high : HARVEST}}
    x0 = vector([0, 0, 0])
    t, x, y = odeint(f, x0, 0.0, 3*math.pi, 0.001, graph, HARVEST, 1e-15, 0.0)

    UC, IM, UK = [vector(_) for _ in zip(*x)]
    ULC = vector(map(UF, t)) - RD*IM - UC
