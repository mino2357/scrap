#!/usr/bin/env python3
"""Plot thermodynamic properties from NASA polynomial data in therm.dat.

The script reads the thermodynamic coefficients, evaluates specific heat,
enthalpy and entropy over the valid temperature range and writes three PNG
files (`cp.png`, `h.png`, `s.png`).
"""
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = 8.3144621  # J/mol/K

def load_thermo(fname):
    with open(fname) as f:
        lines = [line.rstrip('\n') for line in f]
    # locate header with temperature limits
    i = 0
    while i < len(lines):
        if lines[i].strip().startswith('THERMO'):
            i += 1
            t_low, t_mid, t_high = map(float, lines[i].split()[:3])
            i += 1
            break
        i += 1
    data = {}
    while i < len(lines):
        line = lines[i]
        if line.strip().startswith('END'):
            break
        name = line[:16].strip()
        l2 = lines[i+1][:-2]
        l3 = lines[i+2][:-2]
        l4 = lines[i+3][:-2]
        def parse(line):
            fields = []
            for j in range(0, len(line), 15):
                chunk = line[j:j+15].strip()
                if chunk:
                    fields.append(float(chunk))
            return fields
        coeffs = parse(l2) + parse(l3) + parse(l4)
        data[name] = {
            't_mid': t_mid,
            'high': coeffs[:7],
            'low': coeffs[7:]
        }
        i += 4
    return t_low, t_high, data

def thermo_props(T, coeffs):
    c = coeffs['high'] if T >= coeffs['t_mid'] else coeffs['low']
    t = T
    t2, t3, t4 = t*t, t*t*t, t*t*t*t
    cp = R*(c[0] + c[1]*t + c[2]*t2 + c[3]*t3 + c[4]*t4)
    h = R*t*(c[0] + c[1]*t/2 + c[2]*t2/3 + c[3]*t3/4 + c[4]*t4/5 + c[5]/t)
    s = R*(c[0]*math.log(t) + c[1]*t + c[2]*t2/2 + c[3]*t3/3 + c[4]*t4/4 + c[6])
    return cp, h, s

def main():
    base = Path(__file__).resolve().parent
    t_low, t_high, thermo = load_thermo(base / 'therm.dat')
    T = np.linspace(t_low, t_high, 400)
    figs = {
        'cp': plt.figure(),
        'h': plt.figure(),
        's': plt.figure(),
    }
    axes = {k: figs[k].gca() for k in figs}
    for name, coeffs in thermo.items():
        cp_vals = []
        h_vals = []
        s_vals = []
        for temp in T:
            cp_i, h_i, s_i = thermo_props(temp, coeffs)
            cp_vals.append(cp_i)
            h_vals.append(h_i)
            s_vals.append(s_i)
        axes['cp'].plot(T, cp_vals, label=name)
        axes['h'].plot(T, h_vals, label=name)
        axes['s'].plot(T, s_vals, label=name)
    axes['cp'].set_xlabel('Temperature [K]')
    axes['cp'].set_ylabel('cp [J/mol/K]')
    axes['h'].set_xlabel('Temperature [K]')
    axes['h'].set_ylabel('h [J/mol]')
    axes['s'].set_xlabel('Temperature [K]')
    axes['s'].set_ylabel('s [J/mol/K]')
    for ax in axes.values():
        ax.legend()
        ax.grid(True)
    figs['cp'].tight_layout()
    figs['h'].tight_layout()
    figs['s'].tight_layout()
    figs['cp'].savefig(base / 'cp.png')
    figs['h'].savefig(base / 'h.png')
    figs['s'].savefig(base / 's.png')

if __name__ == '__main__':
    main()
