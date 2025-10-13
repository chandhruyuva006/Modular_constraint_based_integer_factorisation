import math
import time
import ast
import sympy
import random
from itertools import product
from tabulate import tabulate
from math import gcd
from functools import reduce

def modinv(a, m):
    """
    Compute modular inverse of a modulo m, using Extended Euclidean Algorithm.
    Returns x such that (a*x) % m == 1
    Raises ValueError if inverse does not exist.
    """
    g, x, y = extended_gcd(a, m)
    if g != 1:
        raise ValueError(f"No modular inverse for {a} mod {m}")
    return x % m

def extended_gcd(a, b):
    """Extended Euclidean Algorithm. Returns (gcd, x, y) such that a*x + b*y = gcd"""
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = extended_gcd(b % a, a)
        return (g, x - (b // a) * y, y)

def merge_residues(M1, R1, M2, R2):
    """
    Merge two residue sets using CRT (handles non-coprime moduli)
    Returns (M_eff, R_eff) where:
      - M_eff = lcm(M1, M2)
      - R_eff = list of residues modulo M_eff compatible with both sets
    """
    g = gcd(M1, M2)
    lcm = (M1 * M2) // g
    R_eff = []

    for r1, r2 in product(R1, R2):
        if (r1 - r2) % g != 0:
            # incompatible residues for non-coprime moduli
            continue
        # CRT formula for coprime or compatible residues
        # solve b ≡ r1 mod M1, b ≡ r2 mod M2
        if g == 1:
            # simple case: coprime
            inv = modinv(M1, M2)
            x = (r1 + M1 * ((r2 - r1) * inv % M2)) % lcm
            R_eff.append(x)
        else:
            # non-coprime but compatible: solve scaled system
            # divide everything by g
            M1g, M2g = M1//g, M2//g
            r1g, r2g = r1//g, r2//g
            try:
                inv = modinv(M1g, M2g)
                xg = (r1g + M1g * ((r2g - r1g) * inv % M2g)) % (M1g*M2g)
                x = xg * g + (r1 % g)
                R_eff.append(x)
            except ValueError:
                # modular inverse does not exist, skip
                continue
    R_eff = sorted(set(R_eff))  # remove duplicates
    return lcm, R_eff

# ===== Example usage =====

# Example residue sets
##M1, R1 = 210, [42, 48, 78, 132, 162, 168]
##M2, R2 = 17, [0, 2, 3, 7, 8, 9, 10, 14, 15]
##
##M_eff, R_eff = merge_residues(M1, R1, M2, R2)

#print(M_eff, R_eff, len(R_eff))


# ----------------- Load modular constraints -----------------
with open("signature_2310.txt", "r") as f:
    mod2310_data = ast.literal_eval(f.read())
mod2310_data = {int(k): v for k, v in mod2310_data.items()}

with open("signature_17.txt", "r") as f:
    mod17_data = ast.literal_eval(f.read())
mod17_data = {int(k): v for k, v in mod17_data.items()}

with open("signature_13.txt", "r") as f:
    mod13_data = ast.literal_eval(f.read())
mod13_data = {int(k): v for k, v in mod13_data.items()}

with open("signature_360.txt", "r") as f:
    mod360_data = ast.literal_eval(f.read())
mod360_data = {int(k): v for k, v in mod360_data.items()}

# ----------------- Two-point factoring -----------------


import math
import time

def two_point_factor(num):
    iterations = 0
    x = 3
    if num % 2 == 0:
        return None  # Early exit if num is even
    
    b_max = num // x
    sqrt_num = math.isqrt(num)
    b_min = 2 * sqrt_num
    

    # Get residues safely, fallback if missing
    n_mod_2310 = num % 2310
    residues_2310 = mod2310_data.get(n_mod_2310, list(range(2310)))

    n_mod_17 = num % 17
    residues_17 = mod17_data.get(n_mod_17, list(range(17)))

    M_eff, R_eff = merge_residues(2310, residues_2310, 17, residues_17)

    n_mod_13 = num % 13
    residues_13 = mod13_data.get(n_mod_13, list(range(13)))

    M_eff, R_eff = merge_residues(M_eff, R_eff, 13, residues_13)
    print(M_eff, len(R_eff))

    n_mod_360 = num % 360
    residues_360 = mod360_data.get(n_mod_360, list(range(360)))

    M_eff, R_eff = merge_residues(M_eff, R_eff, 360, residues_360)
    print(M_eff, len(R_eff))
    
    step_size = M_eff

    b_start = (b_min // step_size) - 1
    start_time = time.perf_counter()

    # ---------------- Optimized Main Loop ----------------
    while b_min <= b_max:
        for i in range(0, num // 3):
            p = 3
            if num % p == 0:
                elapsed = time.perf_counter() - start_time
                return ("trial", i, num // i, elapsed, iterations)
            
            # Precompute the initial value of b
            base_b = (b_start + i) * step_size

            iterations += 1
            
            for r in R_eff:
                iterations += 1
                b = base_b + r
                d = b * b - 4 * num  # b^2 - 4*num
                
                if d >= 0:
                    root = math.isqrt(d)  # Integer square root of d
                    if root * root == d:  # Check if d is a perfect square
                        factor1 = (b + root) // 2
                        factor2 = (b - root) // 2
                        elapsed = time.perf_counter() - start_time
                        return ("b_candidate", factor1, factor2, elapsed, iterations)
            
            p += 2

    return None



# ----------------- Generate semiprimes -----------------
primes = [sympy.randprime(10**9, 10**10) for _ in range(50)]
semiprimes = [random.choice(primes) * random.choice(primes) for _ in range(5)]




#----------------- Factorization and table -----------------
results = []
print("--- Starting Factorization ---")
for N in semiprimes:
    res = two_point_factor(N)
    if res:
        point, f1, f2, elapsed, iters = res
        results.append([N, point, f1, f2, f"{elapsed:.6f}", iters])
    else:
        results.append([N, "None", "-", "-", "-", "-"])

# ----------------- Print formatted table -----------------
headers = ["N", "point", "Factor1", "Factor2", "Time(s)", "Iterations"]
print(tabulate(results, headers=headers, tablefmt="grid"))
