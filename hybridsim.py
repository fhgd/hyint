#!/usr/bin/env python
# -*- coding: utf-8 -*-

def fsolve(g, x0, x1, eps):
    """Nullstellensuche mit dem Sekantenverfahren."""

    while abs(g(x1)) > eps:
        x0, x1 = x1, x1 - g(x1) * (x1 - x0) / (g(x1) - g(x0))
    return x1

def odeint(f, events, actions, x0, t0, t1, dt, eps):
    """Integration eines hybriden Systems (erster Ordnung)"""

    # start values
    t = [t0]
    x = [x0]

    while t[-1] < t1:
        # integrate until at least one event is detected
        while 1:
            # integrate with a euler step
            x.append(x[-1] + f(t[-1], x[-1]) * dt)
            t.append(t[-1] + dt)
            # test all event functions
            idx = [i for i, ev in enumerate(events) if ev(t[-1], x[-1]) < 0]
            if idx:
                # at least one event is detected
                break

        # find the first event
        x_local = lambda t_: x[-2] + f(t[-2] + t_, x[-2]) * t_
        dt_all = []
        for i in idx:
            ev = events[i]
            ev__x_local = lambda t_: ev(t[-2] + t_, x_local(t_))
            dt_all.append(fsolve(ev__x_local, 0, dt, eps))
        # the first event has the smallest dt
        dt_, i = min((dt_, i) for dt_, i in zip(dt_all, idx))

        # correct the last integration step, which was too big
        x[-1] = x_local(dt_)
        t[-1] = t[-2] + dt_

        # discrete step in the continious state x
        x.append(actions[i](x[-1]))
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

    def action(UC):
        return 0.0

    t, x = odeint(f, [event], [action],
        x0=0.0, t0=0.0, t1=3*math.pi, dt=0.001, eps=1e-15)
