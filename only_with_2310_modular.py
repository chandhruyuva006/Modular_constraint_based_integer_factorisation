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
def two_point_factor(num):

    
    iterations = 0
    x= 3
    if num%2 == 0:
        return None
    b_max= num//x

    sqrt_num = math.isqrt(num)
    b_min = 2 * sqrt_num
    step_size = 2310
    

    # Get residues safely, fallback if missing
    n_mod_2310 = num % 2310
    residues = mod2310_data.get(n_mod_2310, list(range(2310)))

   
    b_start = (b_min//2310)-1
    
    
    start_time = time.perf_counter()
    # ---------------- Synchronous loop ----------------
    while b_min<= b_max:
        start_time = time.perf_counter()
        for i in range(0, num//3):
            p=3
            if num%3 == 0:
                elapsed = time.perf_counter() - start_time
                return ("trial", i, num // i, elapsed, iterations)
            
            
            for r in residues:
                iterations += 1
                b = ((b_start + i) * step_size) +r
                d = b * b - 4 * num
                if d>=0:
                    root = math.isqrt(d)
                    if root * root == d:
                        factor1 = (b + root) // 2
                        factor2 = (b - root) // 2
                        elapsed = time.perf_counter() - start_time
                        return ("b_candidate", factor1, factor2, elapsed, iterations)
            
            p+=2

    return None

# ----------------- Generate semiprimes -----------------
primes = [sympy.randprime(10**6, 10**7) for _ in range(50)]
semiprimes = [random.choice(primes) * random.choice(primes) for _ in range(50)]




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
