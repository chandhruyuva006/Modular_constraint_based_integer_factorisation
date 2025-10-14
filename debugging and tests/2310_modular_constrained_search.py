import math
import time
import ast
import sympy
import random
from itertools import product
from tabulate import tabulate  

# ----------------- Load modular constraints -----------------
with open("signature_2310.txt", "r") as f:
    mod2310_data = ast.literal_eval(f.read())
mod2310_data = {int(k): v for k, v in mod2310_data.items()}

# ----------------- Two-point factoring -----------------
def two_point_factor(num, trial_limit=10000):
    
    iterations = 0

    sqrt_num = math.isqrt(num)
    b_min = 2 * sqrt_num
    step_size = 2310

    # Get residues safely, fallback if missing
    n_mod_2310 = num % 2310
    residues = mod2310_data.get(n_mod_2310, list(range(2310)))

    # Always include 0
    if 0 not in residues:
        residues.insert(0, 0)

    base = ((b_min // step_size) - 1) * step_size
    if base < b_min:
        base += step_size

    # ---------------- Generator for constrained b ----------------
    def b_generator():
        k = 0
        while True:
            for r in residues:
                b_candidate = base + k * step_size + r
                if b_candidate >= b_min:
                    # Debug: print the residue combination being used
                    #print(f"Debug: b_candidate = {b_candidate} (residue {r}, k={k})")
                    yield b_candidate
            k += 1

    b_gen = b_generator()
    trial_iter = iter(range(2, trial_limit))
    start_time = time.perf_counter()
    # ---------------- Synchronous loop ----------------
    while True:
        iterations += 1
        trial_done = b_done = False

        # point 2: constrained b
        try:
            b = next(b_gen)
            d = b * b - 4 * num
            if d >= 0:
                root = math.isqrt(d)
                if root * root == d:
                    factor1 = (b + root) // 2
                    factor2 = (b - root) // 2
                    elapsed = time.perf_counter() - start_time
                    return ("b_candidate", factor1, factor2, elapsed, iterations)
        
        except StopIteration:
            b_done = True
        # point 1: trial division
        try:
            p = next(trial_iter)
            if num % p == 0:
                elapsed = time.perf_counter() - start_time
                return ("trial", p, num // p, elapsed, iterations)
        except StopIteration:
            trial_done = True

        if trial_done and b_done:
            break

    return None

# ----------------- Generate semiprimes -----------------
primes = [sympy.randprime(10**7, 10**8) for _ in range(50)]
semiprimes = [random.choice(primes) * random.choice(primes) for _ in range(5)]

# ----------------- Factorization and table -----------------
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
