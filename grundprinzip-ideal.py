#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from scipy.integrate import cumtrapz

from pylab import *
rc('font',   size=10)
rc('axes',   labelsize=14)
rc('xtick',  labelsize=8)
rc('ytick',  labelsize=8)
rc('legend', fontsize=12)
rc('axes', grid=True)

""" Parameter """
UC0 = 0
R = 36e3
C = 12e-9
I = 100e-6
f = 250

deltaUC = R*I/5

sw = True

def f_IC(IP, UC, R):
    return IP - UC/R

def integrate(IP, t, UC0, R, C, deltaUC):
    N = len(t)
    UC  = zeros(N)                          # saved UC
    E   = zeros(N)
    UC[0] = UC0                             # Init UC
    for n in xrange(N-1):
        UC_ref = R/2.*IP[n]
        if  sw and UC_ref - UC[n] > deltaUC:       # event too low
            UC[n+1] = UC_ref
            E[n+1] = E[n] + C/2*(UC[n]**2-UC[n+1]**2)
        elif sw and UC_ref - UC[n] < -deltaUC:     # event too high
            UC[n+1] = UC_ref
            E[n+1] = E[n] + C/2*(UC[n]**2-UC[n+1]**2)
        else:                               # normal integration
            IC = IP[n] - UC[n]/R
            m_UC = IC/C                     # calc ascent of UC
            dt = t[n+1] - t[n]
            UC[n+1] = UC[n] + m_UC*dt       # calc new UC (Euler)
            E[n+1] = E[n]
    return UC, E

def calc(R, deltaUC):
    t = linspace(0, 1/f, 1e4)
    IP = I*sin(2*pi*f*t)
    UC, E = integrate(IP, t, UC0, R, C, deltaUC)
    return UC, IP, E, t

def plot_result(UC, IP, t, deltaUC):
    #f = figure(figsize=[7/2.54, 6/2.54])
    #ax1 = axes([0.2, 0.75, 0.7, 0.15])
    subplots_adjust(right=0.85)
    ax1 = subplot(311)
    plot(t, IP)
    ylabel(r'$I_\mathrm{P}$')
    yticks([])
    ylim(-1.2*I, 1.2*I)

    #ax2 = axes([0.2, 0.1, 0.7, 0.6], sharex=ax1)
    ax2 = subplot(312, sharex=ax1)
    plot(t, R/2.*IP, '-.', label='$U^*$', color='blue')
    plot(t, R/2.*IP + deltaUC, ':', color='0.6')
    plot(t, R/2.*IP - deltaUC, ':', color='0.6')
    plot(t, UC, label='$U$', color='green')
    ylabel(r'$U$')
    yticks([])
    legend(loc=(1.02, 0.4))
    ylim(-R*I, R*I)

    ax3 = subplot(313, sharex=ax1)
    plot(t[:-1], R/4*cumtrapz(IP**2, t), '-.', label='$E_\mathrm{max}$')
    plot(t, E, label='$E_\mathrm{out}$')
    legend(loc=(1.02, 0.4))
    ylabel(r'$E$')
    yticks([])
    ylim(-2e-8, 2e-7)
    xticks([1/f/4, 1/f/2, 3/f/4, 1/f], ['$T/4$', '$T/2$','$3T/4$', '$T$'])
    xlabel('$t$', x=0.5, ha='center')
    
    """ Layout """
    setp(ax1.get_xticklabels(), visible=False)
    setp(ax2.get_xticklabels(), visible=False)
    for ax in ax1, ax2, ax3:
        ax.set_autoscale_on(False)
        ax.grid(True)
        #ax.set_xlim(ax.dataLim.intervalx)
        #ax.set_ylim(ax.dataLim.intervaly*1.2)
        ax.axhline(linestyle='-', color='k')
        setp(ax.get_yaxis().label, y=0.5, va='center', rotation='horizontal')
        # y-Label bündig mit 0V tick von ax2
        #ax.yaxis.set_label_coords(get_tick_position(ax2.yaxis)[0], 1)
        #postPlot(ax)
    return ax3

if __name__ == '__main__':
    UC, IP, E, t = calc(R, deltaUC)
    plot_result(UC, IP, t, deltaUC)
    show()
