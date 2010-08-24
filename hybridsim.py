#!/usr/bin/env python
# -*- coding: utf-8 -*-

def fsolve(g, x0, x1, eps):
    """Nullstellensuche mit dem Sekantenverfahren."""

    while abs(g(x1)) > eps:
        x0, x1 = x1, x1 - g(x1) * (x1 - x0) / (g(x1) - g(x0))
    return x1

def odeint(f, g, h, x0, t0, t1, dt, eps):
    """Integration eines hybriden Systems (erste Ordnung)"""

    # start values
    t = [t0]
    x = [x0]

    while t[-1] < t1:

        # integrate
        while g(x[-1]) >=0:
            x.append(x[-1] + f(x[-1])*dt)
            t.append(t[-1] + dt)

        # correct the last step, which was too big
        x_local = lambda t_: x[-2] + f(x[-2])*t_
        g__x_local = lambda t_: g(x_local(t_))
        dt_ = fsolve(g__x_local, 0, dt, eps)
        x[-1] = x_local(dt_)
        t[-1] = t[-2] + dt_

        # discrete step in the continious state x
        x.append(h(x[-1]))
        t.append(t[-1])

    return t, x


if __name__ == '__main__':

    def f(x):
        return -x

    def g(x):
        return x - 1

    def h(x):
        return 10.0

    t, x = odeint(f, g, h, x0=10.0, t0=0.0, t1=7.0, dt=0.1, eps=1e-6)
