import os
import numpy as np
import math
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from warnings import catch_warnings

act_img_ratio = 0.1718 # Approximate ratio of length without bucket from specs in inches over length without bucket from image in mm

#act_img_ratio = 1     # Uncomment to use image link lengths instead of approximated actual lengths
# Known values (For the moment just measured in millimeters from image multiplied by ratio)
rba = 142.27*act_img_ratio
rdb = 115.2*act_img_ratio
rca = 49.46*act_img_ratio
thetaca = math.radians(293.04)
rgc = 148.5*act_img_ratio # Can't really tell from picture, but should not affect results as long as counted for accordingly
rde = 85.65*act_img_ratio
ref = 197.12*act_img_ratio
rfc = 267.12*act_img_ratio
thetafc = math.radians(5.44) # angle from positive x-axis to ground link, which is Rea
betad = math.radians(92.47) # Angle from DB to DA
rpb = 439.41*act_img_ratio # distance from point b to the bucket coupler point
betab = math.radians(7.41) # Angle from rpb to rdb


# Inputs
rdc = 173*act_img_ratio
rdotdg = 100 # rate of change of length of hydraulic piston
rddotdg = 0 # rate of change of rate of change of length of hydraulic piston

rdg = rdc - rgc # Given input length

# Finds most common element in a list
def most_common(lst: list):
    return max(set(lst), key=lst.count)


"""
Unknowns:
x[0]: thetaba
x[1]: thetadb
x[2]: thetagc
x[3]: thetadg
x[4]: thetade
x[5]: thetaef
"""
unknown_length_indices = []
used_inputs = []
num_equations = 6

# System of scalar equations to solve
def angle_eqs(x):
    return [
        rba*np.cos(x[0]) + rdb*np.cos(x[1]) - rca*np.cos(thetaca) - rgc*np.cos(x[2]) - rdg*np.cos(x[3]),
        rba*np.sin(x[0]) + rdb*np.sin(x[1]) - rca*np.sin(thetaca) - rgc*np.sin(x[2]) - rdg*np.sin(x[3]),
        rgc*np.cos(x[2]) + rdg*np.cos(x[3])  -rde*np.cos(x[4]) - ref*np.cos(x[5]) - rfc*np.cos(thetafc),
        rgc*np.sin(x[2]) + rdg*np.sin(x[3]) - rde*np.sin(x[4]) - ref*np.sin(x[5]) - rfc*np.sin(thetafc),
        x[2]-x[3],
        x[1]+betad-x[4]
    ]

def angular_velocity_eqs(thetaba, thetadb, thetagc, rdotdg, thetadg, thetade, thetaef, rdc):
    rdg = rdc - rgc
    return [
        [-rba*np.sin(thetaba), -rdb*np.sin(thetadb), rgc*np.sin(thetagc), rdg*np.sin(thetadg), 0, 0],
        [rba*np.cos(thetaba), rdb*np.cos(thetadb), -rgc*np.cos(thetagc), -rdg*np.cos(thetadg), 0, 0],
        [0, 0, -rgc*np.sin(thetagc), -rdg*np.sin(thetadg), rde*np.sin(thetade), ref*np.sin(thetaef)],
        [0, 0, rgc*np.cos(thetagc), rdg*np.cos(thetadg), -rde*np.cos(thetade), -ref*np.cos(thetaef)],
        [0, 0, 1, -1, 0, 0],
        [0, 1, 0, 0, -1, 0]
    ], [
        rdotdg*np.cos(thetadg),
        -rdotdg*np.sin(thetadg),
        -rdotdg*np.cos(thetadg),
        -rdotdg*np.sin(thetadg),
        0,
        0
    ]

def angular_acceleration_eqs(thetaba, thetadb, thetagc, rdotdg, thetadg, thetade, thetaef, rdc, rddotdg, omegaba, omegadb, omegagc, omegadg, omegade, omegaef):
    rdg = rdc - rgc
    return [
        [-rba*np.sin(thetaba), -rdb*np.sin(thetadb), rgc*np.sin(thetagc), rdg*np.sin(thetadg), 0, 0],
        [rba*np.cos(thetaba), -rdb*np.cos(thetadb), -rgc*np.cos(thetagc), -rdg*np.cos(thetadg), 0, 0],
        [0, 0, -rgc*np.sin(thetagc), -rdg*np.sin(thetadg), rde*np.sin(thetade), ref*np.sin(thetaef)],
        [0, 0, rgc*np.cos(thetagc), rdg*np.cos(thetadg), -rde*np.cos(thetade), -ref*np.cos(thetaef)],
        [0, 0, 1, -1, 0, 0],
        [0, 1, 0, 0, -1, 0]
    ], [
        rba*(omegaba**2)*np.cos(thetaba) + rdb*(omegadb**2)*np.cos(thetadb) - rgc*(omegagc**2)*np.cos(thetagc) + rddotdg*np.cos(thetadg) - 2*rdotdg*omegadg*np.sin(thetadg) - rdg*(omegadg**2)*np.cos(thetadg),
        rba*(omegaba**2)*np.sin(thetaba) + rdb*(omegadb**2)*np.sin(thetadb) - rgc*(omegagc**2)*np.sin(thetagc) + rddotdg*np.sin(thetadg) + 2*rdotdg*omegadg*np.cos(thetadg) - rdg*(omegadg**2)*np.sin(thetadg),
        rgc*(omegagc**2)*np.cos(thetagc) - rddotdg*np.cos(thetadg) + 2*rdotdg*omegadg*np.sin(thetadg) + rdg*(omegadg**2)*np.cos(thetadg) - rde*(omegade**2)*np.cos(thetade) - ref*(omegaef**2)*np.cos(thetaef),
        rgc*(omegagc**2)*np.sin(thetagc) - rddotdg*np.sin(thetadg) - 2*rdotdg*omegadg*np.cos(thetadg) + rdg*(omegadg**2)*np.sin(thetadg) - rde*(omegade**2)*np.sin(thetade) - ref*(omegaef**2)*np.sin(thetaef),
        0,
        0
    ]

# Rename pi for shorter call
pi = math.pi

# For finding right guesses
initial_guess_min = [0, 0, 0, 0, 0, 0]
initial_guess_max = [2*pi, 2*pi, 2*pi, 2*pi, 2*pi, 2*pi]
current_guess = [0, 0, 0, 0, 0, 0]

# Create empty set to store unique solutions
gen_solutions = set()

# Check if solution is in approximate range of motion
def result_check(result) -> bool:
    if result[0] > 3*pi/4 or result[0] < pi/4:

        return False
    if result[1] > pi/2 and result[1] < 3*pi/2:

        return False
    if result[2] > pi/2:

        return False
    if result[3] > pi/2:

        return False
    if result[4] > pi:

        return False
    if result[5] < pi/2 or result[5] > pi:

        return False
    return True

# Check error of solution
"""def calc_err(result, eq=eqs):
    err = eq(result)
    err_sum = 0
    for i in range(len(err)):
        err_sum += float(np.abs(err[i]))
    return err_sum"""

def calc_err(result, eq=angle_eqs, unkn_len_ind=unknown_length_indices):
    err = eq(result)
    err_sum = 0
    for i in range(len(err)):
        if i in unkn_len_ind:
            calc_err_amt = float(np.abs(err[i]))
        else:
            calc_err_amt = float(np.abs(err[i])) % 2*pi
            calc_err_amt = min(calc_err_amt, 2*pi - calc_err_amt)
        err_sum += calc_err_amt
    return err_sum
    

# Populate the solution set with valid solutions
def solution_gen(current_guess: list[float] | None = None, guess_per_var: int = 5, index: int = 0, initial_min: list[float] = initial_guess_min, initial_max: list[float] = initial_guess_max, output: set[tuple[float, ...]] = gen_solutions, eq=angle_eqs, num_var: int = num_equations, unkn_len_ind=unknown_length_indices):
    if current_guess is None:
        current_guess = []
    if len(initial_min) != num_var and len(initial_max) != num_var:
        exit("Initial and/or final guess are wrong length")
    
    if index == num_var:
        with catch_warnings(record=True) as recorded_warnings:
            sol = fsolve(eq, current_guess)
            curr_error = calc_err(sol, eq)
            if result_check(sol) and curr_error < 0.1:
                output.add(tuple(sol))
        return
    
    
    guesses = np.linspace(initial_min[index], initial_max[index], guess_per_var)

    # Finds solutions with all different combinations of guesses
    for guess in guesses:
        solution_gen(current_guess=current_guess+[guess], guess_per_var=guess_per_var, index=index+1, initial_min=initial_min, initial_max=initial_max, output=output, eq=eq, num_var=num_var, unkn_len_ind=unkn_len_ind)
    return


# Absolute position of bucket coupler point with point C as the origin, positive x-axis to the right
def link_point_pos(angles):

    point_c = (0, 0)

    point_a = (-rca*math.cos(thetaca), -rca*math.sin(thetaca))

    point_f = (rfc*math.cos(thetafc), rfc*math.sin(thetafc))

    point_b = (point_a[0] + rba*math.cos(angles[0]), point_a[1] + rba*math.sin(angles[0]))

    point_e = (point_f[0] + ref*math.cos(angles[5]), point_f[1] + ref*math.sin(angles[5]))

    point_g = (rgc*math.cos(angles[2]), rgc*math.sin(angles[2]))

    point_d = (point_e[0] + rde*math.cos(angles[4]), point_e[1] + rde*math.sin(angles[4]))

    bucket = (point_b[0] + rpb*math.cos(angles[1]-betab), point_b[1] + rpb*math.sin(angles[1]-betab))

    return point_a, point_b, point_c, point_d, point_e, point_f, point_g, bucket



# Generate solutions
solution_gen()

used_inputs.append(rdc)

# Turn set into a list so it is indexable
gen_solutions = list(gen_solutions)

min_error = 1
min_sol = None
for sol in gen_solutions:
    err = calc_err(sol)
    if err < min_error:
        min_error = err
        min_sol = sol


s = min_sol

if s is None:
    exit("No solutions found")

sol_domain = np.linspace(rdc, 1.9*rgc, 10000)
rdc_store = rdc
sol_set = [s]
for x in sol_domain[1:]:
    rdc = x
    rdg = rdc - rgc
    with catch_warnings(record=True) as record_warnings:
        try:    
            new_sol = fsolve(angle_eqs, sol_set[-1])
            sol_set.append(new_sol)
            used_inputs.append(rdc)
        except:
            continue
    




rdc = rdc_store
rdg = rdc - rgc

point_ax = []
point_ay = []
point_bx = []
point_by = []
point_cx = []
point_cy = []
point_dx = []
point_dy = []
point_ex = []
point_ey = []
point_fx = []
point_fy = []
point_gx = []
point_gy = []
bucket_x = []
bucket_y = []

omegaba = []
omegadb = []
omegagc = []
omegadg = []
omegade = []
omegaef = []

alphaba = []
alphadb = []
alphagc = []
alphadg = []
alphade = []
alphaef = []

ang_vel_sol_set = []
ang_acc_sol_set = []
for i in range(len(sol_set)):
    point_pos = link_point_pos(sol_set[i])
    point_ax.append(point_pos[0][0])
    point_ay.append(point_pos[0][1])
    point_bx.append(point_pos[1][0])
    point_by.append(point_pos[1][1])
    point_cx.append(point_pos[2][0])
    point_cy.append(point_pos[2][1])
    point_dx.append(point_pos[3][0])
    point_dy.append(point_pos[3][1])
    point_ex.append(point_pos[4][0])
    point_ey.append(point_pos[4][1])
    point_fx.append(point_pos[5][0])
    point_fy.append(point_pos[5][1])
    point_gx.append(point_pos[6][0])
    point_gy.append(point_pos[6][1])
    bucket_x.append(point_pos[7][0])
    bucket_y.append(point_pos[7][1])

    ang_vel_coeff, ang_vel_const = angular_velocity_eqs(sol_set[i][0], sol_set[i][1], sol_set[i][2], rdotdg, sol_set[i][3], sol_set[i][4], sol_set[i][5], used_inputs[i])
    ang_vel_sol_set.append(np.linalg.solve(ang_vel_coeff, ang_vel_const))

    ang_acc_coeff, ang_acc_const = angular_acceleration_eqs(sol_set[i][0], sol_set[i][1], sol_set[i][2], rdotdg, sol_set[i][3], sol_set[i][4], sol_set[i][5], used_inputs[i], rddotdg, ang_vel_sol_set[i][1], ang_vel_sol_set[i][1], ang_vel_sol_set[i][2], ang_vel_sol_set[i][3], ang_vel_sol_set[i][4], ang_vel_sol_set[i][5])
    ang_acc_sol_set.append(np.linalg.solve(ang_acc_coeff, ang_acc_const))


ba_ang_vel = [point[0] for point in ang_vel_sol_set]
ba_ang_acc = [point[0] for point in ang_acc_sol_set]

def r_degrs(num):
    return num/pi*180 % 360

s = list(s)
for i in range(len(s)):
    if not i in unknown_length_indices:
        s[i] = r_degrs(s[i])
print(f"Theta_BA = {round(s[0], 3)} degrees")
print(f"Theta_DB = {round(s[1], 3)} degrees")
print(f"Theta_GC = {round(s[2], 3)} degrees")
print(f"Theta_DG = {round(s[3], 3)} degrees")
print(f"Theta_DE = {round(s[4], 3)} degrees")
print(f"Theta_EF = {round(s[5], 3)} degrees")

plot_path = os.path.join(os.path.split(os.path.split(os.path.abspath(__file__))[0])[0], "imgs")

# Link Points Plot
plt.rcParams['font.size'] = 12
fig1, ax1 = plt.subplots(layout='constrained')
ax1.grid()
ax1.plot(bucket_x, bucket_y, 'black', label = 'Bucket Coupler Point')
ax1.plot(point_ax, point_ay, color='orange', marker='o', label = 'Point A')
ax1.plot(point_bx, point_by, 'red', label = 'Point B')
ax1.plot(point_cx, point_cy, color='green', marker='o', label = 'Point C')
ax1.plot(point_dx, point_dy, 'blue', label = 'Point D')
ax1.plot(point_ex, point_ey, 'purple', label = 'Point E')
ax1.plot(point_fx, point_fy, color='gray', marker='o', label = 'Point F')
ax1.plot(point_gx, point_gy, 'maroon', label = 'Point G')
ax1.set_xlim(-20, 80)
ax1.set_ylim(-20, 80)
ax1.set_aspect('equal', adjustable='box')
ax1.set_xlabel('Absolute x-coordinate (inches)')
ax1.set_ylabel('Absolute y-coordinate (inches)')
ax1.legend(bbox_to_anchor=(1.1, 1), loc='upper left')
plt.savefig(os.path.join(plot_path, "LinkPointCoordinates.png"), dpi=400)


plt.rcParams['font.size'] = 12
fig1, ax1 = plt.subplots(layout='constrained')
ax1.grid()
ax1.plot(used_inputs, ba_ang_vel, color='black', linestyle='solid', label = '$\\omega_{ba}$')
ax1.set_xlabel('Input length (inches)')
ax1.set_ylabel('Angular velocity (rad/s)')
ax1.legend(bbox_to_anchor=(1.1, 1), loc='upper left')
plt.savefig(os.path.join(plot_path, "LinkBAAngVel.png"), dpi=400)


plt.rcParams['font.size'] = 12
fig1, ax1 = plt.subplots(layout='constrained')
ax1.grid()
ax1.plot(used_inputs, ba_ang_acc, color='black', linestyle='solid', label = '$\\alpha_{ba}$')
ax1.set_xlabel('Input length (inches)')
ax1.set_ylabel('Angular acceleration (rad/s^2)')
ax1.legend(bbox_to_anchor=(1.1, 1), loc='upper left')
plt.savefig(os.path.join(plot_path, "LinkBAAngAcc.png"), dpi=400)