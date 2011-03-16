"""A pure python module for integration of a hybrid system."""

def fsolve(g, x0, x1, eps):
    """Zero finding of function g with the bisect methode.

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
    """Integration of a hybrid system within the time interval [t0, t1].

    The hybrid system consists of a first order ordinary differential equation
    (ODE) system

        x' = f(t, x, y)    with x(t0) = x0

    and a finite-state machine (FSM)

        z_new = graph(z, event)    with z(t0) = z0.

    Where t is the (global) time, x is the continuous state (vector), y is
    the time discrete state (vector) and z is the FSM state. The state
    transition function of the FSM is given as a dictionary like following:

        graph = {z_1 : {event_a : z_1, event_b : z_2},
                 z_2 : {event_b : z_1, event_b : z_2},
                }

    The finite states z_i are functions which are called after a state
    transition and defines the new init value for the ODE system and the time
    discrete state transition of y:

        x_new, y_new = z(t, x, y)    with y(t0) = y0.

    The events are also functions which must be negative to activate the
    event:

        ev(t, x) < 0.

    It is important to note that only one event of the current state can be
    active. Otherwise the next FSM state would be ambiguous.

    The integration of the ODE system is done with the Runge-Kutta methode of
    fourth order and a stepsize dt until the next event occurs. The exact
    event time t_ev for

        event(t_ev, x(t_ev), y) = 0

    is determined within the tolerance eps. It is also important to note that
    events are only detected if the event function change the sign within one
    time step:

        event(t, x, y) * event(t + dt, x, y) <= 0.

    But the correct event detection is not the only reason for a small step
    size. The used Runge-Kutte methode for integration is not a-stabil. So the
    integration error becomes unbounded if the step size is too large.
    """
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
                        'More then one active event after a FSM transaction.'
                else:
                    break
        # integrate with Runge-Kutta fourth order
        t_, x_, y_ = t[-1], x[-1], y[-1]
        k1 = f(t_, x_, y_)
        k2 = f(t_ + dt/2.0, x_ + k1*dt/2.0, y_)
        k3 = f(t_ + dt/2.0, x_ + k2*dt/2.0, y_)
        k4 = f(t_ + dt, x_ + k3*dt, y_)
        x.append(x_ + (k1 + 2*(k2 + k3) + k4)*dt/6.0)
        t.append(t_ + dt)
        # constant extension of the time discrete vaules
        y.append(y[-1])
    return t, x, y

if __name__ == '__main__':

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

    from numpy import pi, sin, cos, sign
    from vector import vector as array

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