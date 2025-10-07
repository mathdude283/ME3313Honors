import numpy as np
import math
from scipy.optimize import fsolve
from warnings import catch_warnings

#act_img_ratio = 0.1718 # Approximate ratio of length without bucket from specs in inches over length without bucket from image in mm

act_img_ratio = 1
# Known values (For the moment just measured in millimeters from image)
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
rdc = 173

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

"""# System of scalar equations to solve
def eqs(x):
    return [

        rba*np.cos(x[0]) + rdb*np.cos(x[1]) - rca*np.cos(thetaca) - rgc*np.cos(x[2]) - rdg*np.cos(x[3]),
        rba*np.sin(x[0]) + rdb*np.sin(x[1]) - rca*np.sin(thetaca) - rgc*np.sin(x[2]) - rdg*np.sin(x[3]),
        rgc*np.cos(x[2]) + rdg*np.cos(x[3])  -rde*np.cos(x[4]) - ref*np.cos(x[5]) - rfc*np.cos(thetafc),
        rgc*np.sin(x[2]) + rdg*np.sin(x[3]) - rde*np.sin(x[4]) - ref*np.sin(x[5]) - rfc*np.sin(thetafc),
        x[2]-x[3],
        x[1]+betad-x[4]
    ]"""

# test
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

var_choices = []


# Finds best guess for each variable
for var in range(len(current_guess)):
    # Generates list of guesses from min to max
    guesses = np.linspace(initial_guess_min[var], initial_guess_max[var], 1000)

    results = [] # Holds all results including Nones
    results_wo_none = [] # Holds all results except for Nones


    # Finds results from each guess, converting to degrees and constraining to 360 degrees for angles

    for guess in guesses:
        current_guess[var] = guess


        # If no 
        with catch_warnings(record=True) as w:
            sol = fsolve(eqs, current_guess)
            if w:
                results.append(None)
                continue
        if not var in unknown_length_indices:
            sol[var] = math.degrees(sol[var]) % 360
        results.append(float(np.round(sol[var], 3)))
        results_wo_none.append(float(np.round(sol[var], 3)))

    res_set = set(results_wo_none)
    var_choices.append(list(res_set))
    # Sets guess for current variable to guess resulting in most common solution
    current_guess[var] = guesses[results.index(most_common(results_wo_none))]


    



# Finds solution for guesses
s = fsolve(eqs, current_guess)

# Converts angles to degrees from 0 to 360 and rounds all results to 3 decimal places
for numb in range(len(s)):
    if not numb in unknown_length_indices:
        s[numb] = math.degrees(s[numb]) % 360
    s[numb] = round(s[numb], 3)

# Prints results
print(f"Theta_BA = {s[0]} degrees")
print(f"Theta_DB = {s[1]} degrees")
print(f"Theta_GC = {s[2]} degrees")
print(f"Theta_DG = {s[3]} degrees")
print(f"Theta_DE = {s[4]} degrees")
print(f"Theta_EF = {s[5]} degrees")