
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


import csv

def export_combined_table_csv(delta_v, time, theta_A_list, theta_B_list, filename, departure_planet, arrival_planet):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        writer.writerow(["theta_A\\theta_B"] + [str(theta_B) for theta_B in theta_B_list])
        # Write rows
        for i, theta_A in enumerate(theta_A_list):
            row = [str(theta_A)]
            for j in range(len(theta_B_list)):
                dv = delta_v[i][j].strip()
                t = time[i][j].strip()
                if dv == "N/A" or t == "N/A" or dv == "" or t == "":
                    cell = "-"
                else:
                    cell = f"{dv} ({t})"
                row.append(cell)
            writer.writerow(row)

def tabulate_orbital_transfer(departure_planet, arrival_planet, interval=45):
    delta_v_results = []
    time_results = []
    theta_range = list(range(0, 360, interval))
    # Set the true anomaly (in degrees) for departure and arrival planets
    for theta_A_deg in theta_range: # True anomaly of departure planet
        delta_v_row = []
        time_row = []
        for theta_B_deg in theta_range: # True anomaly of arrival planet
            try:
                # Orbit parameters for departure and arrival planets
                r_p_A = planets[departure_planet]['r_p']
                r_a_A = planets[departure_planet]['r_a']
                r_p_B = planets[arrival_planet]['r_p']
                r_a_B = planets[arrival_planet]['r_a']

                # Convert degrees to radians
                theta_A = m.radians(theta_A_deg)
                theta_B = m.radians(theta_B_deg)

                # Calculate orbital elements
                a_A = (r_p_A + r_a_A) / 2
                e_A = (r_a_A - r_p_A) / (r_a_A + r_p_A)
                p_A = a_A * (1 - e_A**2)
                h_A = m.sqrt(p_A * mu)
                a_B = (r_p_B + r_a_B) / 2
                e_B = (r_a_B - r_p_B) / (r_a_B + r_p_B)
                p_B = a_B * (1 - e_B**2)
                h_B = m.sqrt(p_B * mu)

                # Calculate position at theta_A and theta_B
                r_A = p_A / (1 + e_A * m.cos(theta_A))
                r_B = p_B / (1 + e_B * m.cos(theta_B))

                # Solve for transfer orbit parameters using Lambert's problem
                e_t = (r_B - r_A) / (r_A * m.cos(theta_A) - r_B * m.cos(theta_B))
                if abs(e_t) > 1:
                    raise ValueError()
                p_t = r_A * r_B * (m.cos(theta_A) - m.cos(theta_B)) / (r_A * m.cos(theta_A) - r_B * m.cos(theta_B))
                a_t = p_t / (1 - e_t**2)
                h_t = m.sqrt(p_t * mu)
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

                # conver true anomaly to eccentric anomaly
                E_t_A = 2 * m.atan2(m.tan(theta_t_A / 2) * m.sqrt((1 - e_t) / (1 + e_t)), 1)
                E_t_B = 2 * m.atan2(m.tan(theta_t_B / 2) * m.sqrt((1 - e_t) / (1 + e_t)), 1)

                # Convert eccentric anomaly to mean anomaly
                M_t_A = E_t_A - e_t * m.sin(E_t_A)
                M_t_B = E_t_B - e_t * m.sin(E_t_B)

                # Calculate transfer time
                T_t = 2 * m.pi * m.sqrt(a_t**3 / mu)
                if M_t_B >= M_t_A:
                    delta_t = (M_t_B - M_t_A) / (2 * m.pi) * T_t
                else:
                    delta_t = ((M_t_B - M_t_A + 2 * m.pi) / (2 * m.pi)) * T_t

                # Calculate perpendicular velocity components at departure
                v_p_A = h_A / r_A
                v_p_t = h_t / r_A

                # Calculate radial velocity components at departure
                v_r_A = mu / h_A * e_A * m.sin(theta_A)
                v_r_t = mu / h_t * e_t * m.sin(theta_A)

                # Calculate perpendicular velocity components at arrival
                v_p_B = h_B / r_B
                v_p_t_B = h_t / r_B

                # calculate radial velocity components at arrival
                v_r_B = mu / h_B * e_B * m.sin(theta_B)
                v_r_t_B = mu / h_t * e_t * m.sin(theta_B)

                # Calculate total delta-v for the transfer and save results
                Delta_v_A = m.sqrt((abs(v_p_t-v_p_A))**2 + (abs(v_r_t-v_r_A))**2)
                Delta_v_B = m.sqrt((abs(v_p_t_B-v_p_B))**2 + (abs(v_r_t_B-v_r_B))**2)
                Delta_v = Delta_v_A + Delta_v_B
                days = delta_t / (3600 * 24)
                if days > 10000:
                    delta_v_row.append("   N/A  ")
                    time_row.append("   N/A  ")
                else:
                    delta_v_row.append(f"{Delta_v:7.2f}")
                    time_row.append(f"{days:7.0f}")
            except Exception:
                delta_v_row.append("   N/A  ")
                time_row.append("   N/A  ")
        delta_v_results.append(delta_v_row)
        time_results.append(time_row)
    # Print delta-v table
    print(f"Results for transfer from {departure_planet} to {arrival_planet}")
    print(f"Rows correspond to theta_A for location of {departure_planet} at departure and columns correspond to theta_B for location of {arrival_planet} at arrival")
    header = "theta_A\\theta_B" + "".join([f"{theta_B_deg:>9}" for theta_B_deg in theta_range])
    print("Delta-v (km/s):")
    print(header)
    print("-" * len(header))
    for i, theta_A_deg in enumerate(theta_range):
        row_str = f"{theta_A_deg:>14}" + "".join([f"{val:>9}" for val in delta_v_results[i]])
        print(row_str)
    print("\nTransfer time (days):")
    print(header)
    print("-" * len(header))
    for i, theta_A_deg in enumerate(theta_range):
        row_str = f"{theta_A_deg:>14}" + "".join([f"{val:>9}" for val in time_results[i]])
        print(row_str)
    print()
    # Export combined table as CSV
    csv_filename = f"orbital_transfers_{departure_planet.replace(' ', '_')}_to_{arrival_planet.replace(' ', '_')}.csv"
    export_combined_table_csv(delta_v_results, time_results, theta_range, theta_range, csv_filename, departure_planet, arrival_planet)

planets_to_use = ['Mercury', 'Venus', 'Earth', 'Mars', 'The Asteroid Belt']
for dep in planets_to_use:
    for arr in planets_to_use:
        if dep != arr:
            tabulate_orbital_transfer(dep, arr)

