# Import math module
import math as m

# --- Planetary orbit data (Values in km) ---
planets = {
	'Mercury':              {'r_p': 46001200, 'r_a': 69816900},
	'Venus':                {'r_p': 107477000, 'r_a': 108939000},
	'Earth':                {'r_p': 147098290, 'r_a': 152100000},
	'Mars':                 {'r_p': 206700000, 'r_a': 249200000},
    'The Asteroid Belt':    {'r_p': 375000000, 'r_a': 375000000},
	'Jupiter':              {'r_p': 740520000, 'r_a': 816620000},
	'Saturn':               {'r_p': 1353572956, 'r_a': 1513325783},
	'Uranus':               {'r_p': 2734998229, 'r_a': 3006318143},
	'Neptune':              {'r_p': 4452940833, 'r_a': 4553946490},
}

# Standard gravitational parameter for the Sun (mu = GM), km^3/s^2
mu = 1.327e11

# --- User selection ---
# Set the names of the departure and arrival planets
departure_planet = 'Earth'
arrival_planet = 'Mars'

# Input Angles
theta_A_deg = 0
theta_B_deg = 150

# --- End of user selection ---

# --- Orbit parameters for each planet ---
r_p_A = planets[departure_planet]['r_p']
r_a_A = planets[departure_planet]['r_a']
r_p_B = planets[arrival_planet]['r_p']
r_a_B = planets[arrival_planet]['r_a']

# Convert true anomalies to radians
theta_A = m.radians(theta_A_deg)
theta_B = m.radians(theta_B_deg)

# Calculate orbital elements for each planet
a_A = (r_p_A + r_a_A) / 2
e_A = (r_a_A - r_p_A) / (r_a_A + r_p_A)
p_A = a_A * (1 - e_A**2)
h_A = m.sqrt(p_A * mu)

a_B = (r_p_B + r_a_B) / 2
e_B = (r_a_B - r_p_B) / (r_a_B + r_p_B)
p_B = a_B * (1 - e_B**2)
h_B = m.sqrt(p_B * mu)

# Position (radius) at selected true anomaly for each planet
r_A = p_A / (1 + e_A * m.cos(theta_A))
r_B = p_B / (1 + e_B * m.cos(theta_B))

# --- Transfer orbit parameters (assuming common apse line) ---
e_t = (r_B - r_A) / (r_A * m.cos(theta_A) - r_B * m.cos(theta_B))
#if e_t < 0:
    #raise ValueError(f"Transfer orbit eccentricity is negative (e_t={e_t:.6f}). No valid solution for these input angles.")
if e_t > 1:
    raise ValueError(f"Transfer orbit eccentricity is >1. No valid solution for these input angles.")
p_t = r_A * r_B * (m.cos(theta_A) - m.cos(theta_B)) / (r_A * m.cos(theta_A) - r_B * m.cos(theta_B))
a_t = p_t / (1 - e_t**2)
h_t = m.sqrt(p_t * mu)

# --- Calculate true anomaly of departure and arrive for transfer orbit ---
arg_A = (p_t / r_A - 1) / e_t
arg_B = (p_t / r_B - 1) / e_t
tol = 1e-3

if abs(arg_A-1) < tol:  
    arg_A = 1.0
elif abs(arg_A+1) < tol:
    arg_A = -1.0

if abs(arg_B-1) < tol:
    arg_B = 1.0
elif abs(arg_B+1) < tol:
    arg_B = -1.0


if theta_A_deg <= 180:
    theta_t_A = m.acos(arg_A)
else:
    theta_t_A = 2*m.pi - m.acos(arg_A)

    
if theta_B_deg <= 180:
    theta_t_B = m.acos(arg_B)
else:
    theta_t_B = 2*m.pi - m.acos(arg_B)

# --- Convert true anomaly to eccentric anomaly ---
E_t_A = 2 * m.atan2(m.tan(theta_t_A / 2) * m.sqrt((1 - e_t) / (1 + e_t)), 1)
E_t_B = 2 * m.atan2(m.tan(theta_t_B / 2) * m.sqrt((1 - e_t) / (1 + e_t)), 1)

# --- Convert eccentric anomaly to mean anomaly ---
M_t_A = E_t_A - e_t * m.sin(E_t_A)
M_t_B = E_t_B - e_t * m.sin(E_t_B)

# --- Calculate orbital period of transfer orbit ---
T_t = 2 * m.pi * m.sqrt(a_t**3 / mu)  # seconds

# --- Calculate transfer time between anomalies ---
if M_t_B >= M_t_A:
    delta_t = ((M_t_B - M_t_A) / (2 * m.pi)) * T_t
else:
    delta_t = ((M_t_B - M_t_A + 2 * m.pi) / (2 * m.pi)) * T_t

# Perpendicular velocity components at departure point
v_p_A = h_A / r_A
v_p_t_A = h_t / r_A

# Radial velocity components at departure point
v_r_A = mu / h_A * e_A * m.sin(theta_A)
v_r_t_A = mu / h_t * e_t * m.sin(theta_A)

# Total velocity magnitudes at departure point
v_A = m.sqrt(v_p_A**2 + v_r_A**2)
v_t_A = m.sqrt(v_p_t_A**2 + v_r_t_A**2)

# Perpendicular velocity components at arrival point
v_p_B = h_B / r_B
v_p_t_B = h_t / r_B

# Radial velocity components at arrival point
v_r_B = mu / h_B * e_B * m.sin(theta_B)
v_r_t_B = mu / h_t * e_t * m.sin(theta_B)

# Total velocity magnitudes at arrival point
v_B = m.sqrt(v_p_B**2 + v_r_B**2)
v_t_B = m.sqrt(v_p_t_B**2 + v_r_t_B**2)

# Flight path angles (degrees)
phi_A = m.degrees(m.atan2(v_r_A, v_p_A))
phi_t_A = m.degrees(m.atan2(v_r_t_A, v_p_t_A))
phi_B = m.degrees(m.atan2(v_r_B, v_p_B))
phi_t_B = m.degrees(m.atan2(v_r_t_B, v_p_t_B)) 

# Delta-v required for the transfer (magnitude)
Delta_v_A = m.sqrt((abs(v_p_t_A-v_p_A))**2 + (abs(v_r_t_A-v_r_A))**2)
Delta_v_B = m.sqrt((abs(v_p_t_B-v_p_B))**2 + (abs(v_r_t_B-v_r_B))**2)
Delta_v = Delta_v_A + Delta_v_B

# Angle of the required velocity change (degrees)
gamma_A = m.degrees(m.atan2(v_r_t_A - v_r_A, v_p_t_A - v_p_A))
gamma_B = m.degrees(m.atan2(v_r_t_B - v_r_B, v_p_t_B - v_p_B))


print(f"Transfer from {departure_planet} to {arrival_planet}:")
print(f"Delta-v required: {Delta_v:.3f} km/s")     
print(f"Eccentric anomaly at departure (E_A): {E_t_A:.6f} rad")
print(f"Eccentric anomaly at arrival   (E_B): {E_t_B:.6f} rad")  
print(f"Mean anomaly at departure (M_A): {M_t_A:.6f} rad")
print(f"Mean anomaly at arrival   (M_B): {M_t_B:.6f} rad")
print()

# --- Plotting ---
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(8,8))

# Plot the Sun
ax.plot(0, 0, 'yo', label='Sun', markersize=12)

# Function to get orbit points
def get_orbit(a, e, num=500):
    theta = np.linspace(0, 2*np.pi, num)
    r = a * (1 - e**2) / (1 + e * np.cos(theta))
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

# Plot departure planet orbit
xA, yA = get_orbit(a_A, e_A)
ax.plot(xA, yA, 'b-', label=f'{departure_planet} Orbit')

# Plot arrival planet orbit
xB, yB = get_orbit(a_B, e_B)
ax.plot(xB, yB, 'g-', label=f'{arrival_planet} Orbit')

# Plot transfer orbit if it exists
xT, yT = get_orbit(a_t, e_t)
ax.plot(xT, yT, 'r--', label='Transfer Orbit')

# Plot current positions of planets
xA_pos = r_A * np.cos(theta_A)
yA_pos = r_A * np.sin(theta_A)
ax.plot(xA_pos, yA_pos, 'bo', label=f'{departure_planet} (Location at departure)')

xB_pos = r_B * np.cos(theta_B)
yB_pos = r_B * np.sin(theta_B)
ax.plot(xB_pos, yB_pos, 'go', label=f'{arrival_planet} (Location at arrival)')

# Set aspect ratio and labels
ax.set_aspect('equal')
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
ax.set_title('Orbits and Transfer')
ax.legend()
ax.grid(True)

# Set limits to show all orbits
max_r = max(r_p_A, r_a_A, r_p_B, r_a_B, r_A, r_B)
ax.set_xlim(-1.1*max_r, 1.1*max_r)
ax.set_ylim(-1.1*max_r, 1.1*max_r)

plt.show()


