# hyint

hyint is a pure python module for integration of a hybrid system (aka Analog-Mixed-Signal system).

----

# Example: Piezo Harvester with the SECE Interface

SECE stands for: Synchronous electrical charge extraction which do two things:

    1. Leave the harvester in open circuit
    2. When the generator voltage U has an extremum, then immediately
    discharge the piezo capacitor Cp. (Idealy, this discharge doesn't
    take any time.)

The SECE interface is modeled with a simple finite state machine (FSM)
and the piezo harvester is modeled with a circuit model

```text
    ,-->-- R -- L -- C ----+-->--o
   +|  Iq           +UC-   |  I  +
   Uq                      Cp    U
   -|                      |     -
    '----------------------+-----o
```

The simulation result shows the periodic steady state of the total system


----

# Theory

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
