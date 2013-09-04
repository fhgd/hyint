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
        if g0*gm < 0:
            x1 = xm     # zero is in left half
        else:
            x0 = xm     # zero is in right half
            g0 = gm
    return x0, x1

def odestep(f, t, dt, x, y):
    """One time step dt with Runge-Kutta fourth order"""
    dt_2 = dt*0.5
    t_dt_2 = t + dt_2
    k1 = f(t, x, y)
    k2 = f(t_dt_2, x + k1*dt_2, y)
    k3 = f(t_dt_2, x + k2*dt_2, y)
    k4 = f(t + dt, x + k3*dt, y)
    return x + (k1 + 2*(k2 + k3) + k4)*dt*0.16666666666666666

def hyint(f, x0, t0, t1, dt, graph, z0, eps, y0, debug=False):
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

        ev(t, x, y) < 0.

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
        events = [ev for ev in graph[z].iterkeys() if ev(t[-1], x[-1], y[-1]) < 0]
        if events:
            # at least one event is detected, so find the first event
            if debug:
                ev_names = ', '.join([ev.__name__ for ev in events])
                print 't = %.15f : detected %s' % (t[-1], ev_names)
            k_min = 1.0
            event = None
            for ev in events:
                def ev_local(k):
                    x_21 = odestep(f, t[-2], k*dt, x[-2], y[-2])
                    return ev(t[-2] + k*dt, x_21, y[-2])
                k0, k1 = fsolve(ev_local, 0.0, 1.0, eps)
                # Use k1 which terminates the actual process
                if k1 <= k_min:
                    k_min = k1
                    event = ev
            assert event, 'No zero was found in [0, 1]'
            # correct the last integration step, which was too far
            x[-1] = odestep(f, t[-2], k_min*dt, x[-2], y[-2])
            t[-1] = t[-2] + k_min*dt
            if debug:
                print 't = %.15f : located  %s' % (t[-1], event.__name__)
            # transition and action of the FSM until all events are off
            while 1:
                z = graph[z][event]
                if debug:
                    print z.__name__
                x_new, y_new = z(t[-1], x[-1], y[-1])
                x.append(x_new)
                y.append(y_new)
                t.append(t[-1])
                # test all event functions
                events = [ev for ev in graph[z].iterkeys() if ev(t[-1], x[-1], y[-1]) < 0]
                if events:
                    event = events[0]
                    assert len(events) == 1, \
                        'More then one active event after a FSM transaction.'
                else:
                    break
        # integrate one time step
        x.append(odestep(f, t[-1], dt, x[-1], y[-1]))
        t.append(t[-1] + dt)
        # constant extension of the time discrete vaules
        y.append(y[-1])
    return t, x, y
