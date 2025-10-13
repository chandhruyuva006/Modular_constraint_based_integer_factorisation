import math
import time
import ast
import sympy
import random
from itertools import product
from tabulate import tabulate  # pip install tabulate if needed

# ----------------- Load modular constraints -----------------
def load_mod_file(filename):
    """Load a modular signature file into a dictionary with int keys."""
    with open(filename, "r") as f:
        data = ast.literal_eval(f.read())
    data = {int(k): v for k, v in data.items()}

    # Ensure 0 is always included in the residue list (not the key!)
    for key in data:
        if 0 not in data[key]:
            data[key].insert(0, 0)
    return data

# Load your files
mod2310_data = load_mod_file("2310 signature.txt")
mod13_data   = load_mod_file("13 signature.txt")
mod17_data   = load_mod_file("17 signature.txt")

# List of (modulus, residue_dict)
moduli_data = [
    (2310, mod2310_data),
    (13, mod13_data),
    (17, mod17_data)
]

# ----------------- Two-horse factoring -----------------
def two_horse_factor(num, trial_limit=10000):
    start_time = time.perf_counter()
    iterations = 0

    sqrt_num = math.isqrt(num)
    b_min = 2 * sqrt_num

    # Prepare residues for each modulus
    residues_list = []
    for M, data in moduli_data:
        n_mod = num % M
        residues = data.get(n_mod, list(range(M)))
        # Ensure 0 is always included
        if 0 not in residues:
            residues.insert(0, 0)
        residues_list.append((M, residues))

    # Step size: product of all moduli (naive; can optimize later with CRT)
    step_size = 1
    for M, _ in residues_list:
        step_size *= M

    # ---------------- Generator for constrained b ----------------
    def b_generator():
        k = 0
        while True:
            # Generate all combinations of residues for debugging
            for combo in product(*[res for _, res in residues_list]):
                b_candidate = b_min + k * step_size + sum(combo)
                if b_candidate >= b_min:
                    # Debug: print the combination used
                    #print(f"Debug: b_candidate={b_candidate}, combo={combo}, k={k}")
                    yield b_candidate
            k += 1

    b_gen = b_generator()
    trial_iter = iter(range(2, trial_limit))

    # ---------------- Synchronous loop ----------------
    while True:
        iterations += 1
        trial_done = b_done = False

        # Horse 1: trial division
        try:
            p = next(trial_iter)
            if num % p == 0:
                elapsed = time.perf_counter() - start_time
                return ("trial", p, num // p, elapsed, iterations)
        except StopIteration:
            trial_done = True

        # Horse 2: constrained b
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
    res = two_horse_factor(N)
    if res:
        horse, f1, f2, elapsed, iters = res
        results.append([N, horse, f1, f2, f"{elapsed:.6f}", iters])
    else:
        results.append([N, "None", "-", "-", "-", "-"])

# ----------------- Print formatted table -----------------
headers = ["N", "Horse", "Factor1", "Factor2", "Time(s)", "Iterations"]
print(tabulate(results, headers=headers, tablefmt="grid"))
