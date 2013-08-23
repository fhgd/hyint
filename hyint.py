"""A pure python module for integration of a hybrid system."""

"""
ToDo:

* store not the state x, but only output y = Cx due to memory issue

* Do the saving outside, remaining if crashed

* Add warning, when the ode has only one time step between the FSM transitions

* Replace graph[z].iterkeys() and graph[z][event] with

    znew, events = fsm(z, ev)

    def fsm(z, ev):
        znew = graph[z][ev]
        events = graph[znew].iterkeys()
        return znew, events

* Parallel FSM:

    def make_fsm(graphs):
        def fsm(Z, ev):
            '''Return the new state and event list for all FSMs in a dict

                Znew = {z1 : events, z2 : events, ...}
            '''
            Znew = dict((g[z][ev], g[znew].iterkeys()) for g, z in zip(graphs, Z))
            return Znew
        return fsm


* Problem-1: The events from different FSMs could be true at the same time

* Problem-2: How to deal with the setting of new init values?

    x0 = z1(x)
    x0 = z2(x)
    ...

  Each fsm should only set his own components in x! Then the succesive
  setting of the new init values is (should be) no problem. But the different
  x0 must be merge in a good way. Maybe

    x0 = z1(x), z2(x), ...

  should be good enough. Or with a dict:

    z1(x) = {x1 : ?, x2 : ?}
    z2(x) = {x3 : ?, x4 : ?}

    x0 = x.update(z1(x), z2(x))

  Maybe use namedtuple as vectors?

  Or each FSM should have his own ODE which could be empty?
  (And each ODE shoulf have his own FSM which could be empty?)

  In general: Each FSM could read all variables but should have his own set
  of writeable variables (inital conditions and time discrete vars).

    So, maybe for now, only the firste FSM can change the initial values
    of the ODE. The other FSM acting on theire own time discrete variables:

        x0, y1 = z1(x, y1, ...)
        y2 = z2(x, y2, ...)

    or
        x0 = z1(x, y1, ...)
        y2 = z2(x, y2, ...)


* General structure of functions?

    def STATE(x, params):
        return xnew


"""

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
            N = 0
            while 1:
                N = N + 1
                assert N < 10, 'More the 10 FSM tansitions at the same time.'
                z = graph[z][event]
                if debug:
                    print event.__name__, '=>', z.__name__
                x_new, y_new = z(t[-1], x[-1], y[-1])
                x.append(x_new)
                y.append(y_new)
                t.append(t[-1])
                # test all event functions
                events = [ev for ev in graph[z].iterkeys() if ev(t[-1], x[-1], y[-1]) < 0]
                if events:
                    event = events[0]
                    assert len(events) == 1, \
                        'More then one active event after a FSM transaction:\n    %s' \
                        % ', '.join([f.__name__ for f in events])
                else:
                    break
        # integrate one time step
        x.append(odestep(f, t[-1], dt, x[-1], y[-1]))
        t.append(t[-1] + dt)
        # constant extension of the time discrete vaules
        y.append(y[-1])
    return t, x, y
