import numpy as np
import math
from scipy.optimize import fsolve
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


# Inputs
rdc = 173*act_img_ratio

#print(rba, rdb, rca, rgc, rde, ref, rfc, rdc)
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

num_equations = 6

# System of scalar equations to solve
def eqs(x):
    return [

        rba*np.cos(x[0]) + rdb*np.cos(x[1]) - rca*np.cos(thetaca) - rgc*np.cos(x[2]) - rdg*np.cos(x[3]),
        rba*np.sin(x[0]) + rdb*np.sin(x[1]) - rca*np.sin(thetaca) - rgc*np.sin(x[2]) - rdg*np.sin(x[3]),
        rgc*np.cos(x[2]) + rdg*np.cos(x[3])  -rde*np.cos(x[4]) - ref*np.cos(x[5]) - rfc*np.cos(thetafc),
        rgc*np.sin(x[2]) + rdg*np.sin(x[3]) - rde*np.sin(x[4]) - ref*np.sin(x[5]) - rfc*np.sin(thetafc),
        x[2]-x[3],
        x[1]+betad-x[4]
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
    if result[0] > 135 or result[0] < 45:

        return False
    if result[1] > 90 and result[1] < 270:

        return False
    if result[2] > 90:

        return False
    if result[3] > 90:

        return False
    if result[4] > 180:

        return False
    if result[5] < 90 or result[5] > 180:

        return False
    return True

# Check error of solution
def calc_err(result, eq=eqs):
    err = eq(result)
    err_sum = 0
    for i in range(len(err)):
        err_sum += float(np.abs(err[i]))
    return err_sum
    

# Populate the solution set with valid solutions
def solution_gen(current_guess: list[float] | None = None, guess_per_var: int = 5, index: int = 0, initial_min: list[float] = initial_guess_min, initial_max: list[float] = initial_guess_max, output: set[tuple[float, ...]] = gen_solutions, eq=eqs, num_var: int = num_equations, unkn_len_ind=unknown_length_indices):
    if current_guess is None:
        current_guess = []
    if len(initial_min) != num_var and len(initial_max) != num_var:
        exit("Initial and/or final guess are wrong length")
    
    if index == num_var:
        with catch_warnings(record=True):
            sol = fsolve(eq, current_guess)
            curr_error = calc_err(sol, eq)
            for numb in range(len(sol)):
                if not numb in unkn_len_ind:
                    sol[numb] = math.degrees(sol[numb]) % 360
                sol[numb] = round(sol[numb], 1)
            if result_check(sol) and curr_error < 1:
                output.add(tuple(sol))
        return
    
    guesses = np.linspace(initial_min[index], initial_max[index], guess_per_var)

    # Finds solutions with all different combinations of guesses
    for guess in guesses:
        solution_gen(current_guess=current_guess+[guess], guess_per_var=guess_per_var, index=index+1, initial_min=initial_min, initial_max=initial_max, output=output, eq=eq, num_var=num_var, unkn_len_ind=unkn_len_ind)
    return

# Generate solutions
solution_gen()

# Turn set into a list so it is indexable
gen_solutions = list(gen_solutions)


# Grab first solution (There will probably only be one)
s = gen_solutions[0]


# Prints results
print(f"Theta_BA = {s[0]} degrees")
print(f"Theta_DB = {s[1]} degrees")
print(f"Theta_GC = {s[2]} degrees")
print(f"Theta_DG = {s[3]} degrees")
print(f"Theta_DE = {s[4]} degrees")
print(f"Theta_EF = {s[5]} degrees")