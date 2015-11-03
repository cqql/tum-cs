#!/usr/bin/env python3

from math import sin, cos

# Parameters
M_s = 250.0
M_p = 70.0
R = 2.5
I = M_p * R**2
G = 9.81

# Initial conditions
s = 0.0
v = 0.0
phi = 1.25
omega = 0.0

# Time step
h = 0.001
time = 0.0

while time < 5:
    ds_dt = v
    dphi_dt = omega
    dv_dt = (M_p / M_s * (omega**2 / R * sin(phi) - G * sin(phi) * cos(phi))) / (1 + ((cos(phi))**2 + 1) * M_p / M_s)
    domega_dt = (G * sin(phi) - dv_dt * cos(phi)) / R

    s += h * ds_dt
    phi += h * dphi_dt
    v += h * dv_dt
    omega += h * domega_dt
    
    time += h
    print("{}\t{}\t{}".format(time, v, omega))

