#!/usr/bin/env python
# -*- coding: utf-8 -*-

def fsolve(g, x0, x1, eps):
    """Nullstellensuche mit dem Sekantenverfahren."""

    while abs(g(x1)) > eps:
        x0, x1 = x1, x1 - g(x1) * (x1 - x0) / (g(x1) - g(x0))
    return x1

def odeint(f, x0, t0, t1, dt, graph, z0, eps):
    """Integration eines hybriden Systems (erster Ordnung)"""

    # start values
    t = [t0]
    x = [x0]
    z = z0

    while t[-1] < t1:

        # integrate until at least one event is detected
        while 1:
            # integrate with an euler step
            x.append(x[-1] + f(t[-1], x[-1]) * dt)
            t.append(t[-1] + dt)
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
        x.append(z(t[-1], x[-1]))
        t.append(t[-1])
    return t, x


if __name__ == '__main__':

    R = 0.5
    C = 1.0

    import math

    def IP(t):
        return math.sin(t)

    def f(t, UC):
        IC = IP(t) - UC/R
        return IC/C

    def event(t, UC):
        IC = abs(IP(t)) - abs(UC)/R
        return IC

    def action(t, UC):
        return 0.0

    graph = {action : {event : action}}
    t, x = odeint(f, 0.0, 0.0, 3*math.pi, 0.001, graph, action, 1e-15)