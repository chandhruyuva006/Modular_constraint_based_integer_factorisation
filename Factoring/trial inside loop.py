# This file is part of the Integer Factorisation project.
# 
# Integer Factorisation project is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
#  any later version.
#
# Integer Factorisation project is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#


import math
import time
import ast
import sympy
import random
from itertools import product
from tabulate import tabulate
from math import gcd
from functools import reduce


"""
###############################################################################
The current code is a crude product from the initial analysis and completely
limited by my understanding of math and python concepts. Acidically a Biology major.
But is should complete the run more than not. It is just an accelerated fermat like 
method with a better if (residue in allowed residue) loop loop loop. 
I welcome anyone who could understand the constraints and help me optimise the code. 
                                    Thank You. 
############################################################################### 

"""

"""
This code is based on the fact that. 
for any odd prime p and q they will give either 1 or 3 mod 4.
for example p = 47 q =61. p = 3 mod 4 and q = 1 mod 4. 
if the semiprime is made up of factors which are odd.
then N should be either 1 mod 4 or 3 mod 4.
so if the N is 3 mod 4 then the factor sum p+q should be divided by 4.
if the N is 1 mod 4 then the factor sum p+q should be divided by 2.

but here we extended to other composite and primorial moduli.
"""
"""
        any q 
    m  |  q mod m  |  q mod m
    ---------------------------
any  p |(p*q)mod m | (p*q)mod m
p    p |(p*q)mod m | (p*q)mod m

for the N (mod m) cant expect any value,  (p*q)mod m = N mod m  value is defined based on 
factors p mod m and q mod m.    

so is the p+q mod m. it is also restricted to a handful of values only.

for example when it is an odd semiprime so the factors should also be odd.
their sum (p+q) should be always even. 
we cant expect an odd factor sum. which  is just saying the factor sum is 
restricted to 0 mod 2.

we can extend to this any chosen and all moduli under N.

lets say an odd semiprime. N =r (mod 210). the r can never be a even number or a multiple of 5.
if (mod m)N = r is an even number or 5 then the N should be a multiple of 2 or 5.
so the we only have about  less than 105 residues for r when  N = r mod 210. 

so is the factor sum p+q values.  but they are restricted to even further down.

"""

"""
Let me explain the method using a trivial and well understood case.
The possibilities of odd semiprime N mod 4 is either 1 or 4.

modular multiplication table for moduli 4
this table gives us all the possible combinations of p and q and r = N mod 4 space.
           any q
       4 | 1 | 3
     -----------
any    1 | 1 | 3  -  if N (mod 4) = 3 then the factor sum should be 0.
p      3 | 1 | 1      -
                         -
                            -
modular addition table for    -
the moduli 4 in the factor    -
sum space.                   -  
                            - 
the factor sums           -
       4 | 1 |  3        -
      -----------       - 
       1 | 2 |  0 <---   
       3 | 0 |  2
       
so if the N mod 4 is 1 then the factor sum is either 0 mod 4 or 2 mod 4. 
But if the N mod 4 is 3 then the factor sum is strictly restricted to just 0 mod 4.
and we also have one more information about the factors. one factor 1 mod 4 and another 3 mod 4.
But when the moduli become large so are the combinations and
this combinatorial information is not so useful in the factoring process.

searching for through primes is hard because we dont know which is prime
and which is not. but searching for a composite (p+q) is much easy as we 
know most of the things about it from the modular arithmetic.

But we can use the factor sum information,so we dont to check other numbers 
while searching in the (2*(sqrt(N)) to N) space for 
the factor sum when using the x^2 - bx +N =0 equation to find factors.

"""
""" 
The key fact we can use any moduli to reduce the search space.
But the problem is when the moduli is a prime like 7,13,1009. 
on average they give about (m+1)/2 residues. and we can combine multiple 
the residues of different moduli using CRT as they are coprime all the 
residue combinations are allowed for the factor sum. But the reduced search space
vs the CRT combination explosion will nullify the effective factor process.

But when we combine a highly composite moduli like a primorial with a small 
prime moduli, we will get a bigger M_eff and in most cases some of the residus 
cancel out in the CRT process due to no possible number which can satisfy 
both constraints and the allowed factor sum residues will be a tiny fraction
of the combined modulus. Thus reducing the search space and time.
This is the heart of the semiprime modular restriction conjecture. 

The previous investigations on this matter may brushed this off because of the 
fact that we cant combine more and more moduli to get a bigger modulus to 
reduce the search space because as we add more moduli the residue possibilites
get out of hand quickly. also the reason previous investigations may have abandoned it 
because it cant compete with the heuristic algorithms interms of actual run time.
and also they may tried only prime moduli for reducing the search space which i also
did in my previous implementations. 

The key lies on the composites and primorials which in combinations can quickly
lead to the exact factor sum and dont exponentially blow up the combinations as 
some residues which are always incompatible during the CRT process. 

"""




# modinv never failed during the 1000s of runs. But may fail sometime due to the 
# combination of two incompatible composite moduli.
# we can stack any moduli but primorial moduli filters better.
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
            # incompatible residues for non-coprime moduli CRT with non coprime moduli theory cehck
            continue
        # CRT formula for coprime or compatible residues
        # solve b ≡ r1 mod M1, b ≡ r2 mod M2
        if g == 1:
            # simple case: coprime, just use regular CRT
            inv = modinv(M1, M2)
            x = (r1 + M1 * ((r2 - r1) * inv % M2)) % lcm
            R_eff.append(x)
        else:
            # non-coprime but compatible: solve scaled system
            # divide everything by g lcm method ,,,,, needed to check theory for the edge cases
            M1g, M2g = M1//g, M2//g
            r1g, r2g = r1//g, r2//g
            try:
                inv = modinv(M1g, M2g)
                xg = (r1g + M1g * ((r2g - r1g) * inv % M2g)) % (M1g*M2g)
                x = xg * g + (r1 % g)
                R_eff.append(x)
            except ValueError:
                # modular inverse does not exist, skip ,, use to resolve like 2,4 or 6,12 like this
                continue
    R_eff = sorted(set(R_eff))  # remove duplicates other wise next step computer boom
    return lcm, R_eff

# ===== Example usage =====

# Example residue sets to check if the code is working
##M1, R1 = 210, [42, 48, 78, 132, 162, 168]
##M2, R2 = 17, [0, 2, 3, 7, 8, 9, 10, 14, 15]
##
##M_eff, R_eff = merge_residues(M1, R1, M2, R2)

#print(M_eff, R_eff, len(R_eff))


# ----------------- Load modular constraints -----------------

with open("signature_3.txt", "r") as f:
    mod3_data = ast.literal_eval(f.read())
mod3_data = {int(k): v for k, v in mod3_data.items()}

with open("signature_4.txt", "r") as f:
    mod4_data = ast.literal_eval(f.read())
mod4_data = {int(k): v for k, v in mod4_data.items()}

with open("signature_30.txt", "r") as f:
    mod30_data = ast.literal_eval(f.read())
mod30_data = {int(k): v for k, v in mod30_data.items()}

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

with open("signature_29.txt", "r") as f:
    mod29_data = ast.literal_eval(f.read())
mod29_data = {int(k): v for k, v in mod29_data.items()}

with open("signature_36.txt", "r") as f:
    mod36_data = ast.literal_eval(f.read())
mod36_data = {int(k): v for k, v in mod36_data.items()}

with open("signature_16.txt", "r") as f:
    mod16_data = ast.literal_eval(f.read())
mod16_data = {int(k): v for k, v in mod16_data.items()}

with open("signature_18.txt", "r") as f:
    mod18_data = ast.literal_eval(f.read())
mod18_data = {int(k): v for k, v in mod18_data.items()}

with open("signature_20.txt", "r") as f:
    mod20_data = ast.literal_eval(f.read())
mod20_data = {int(k): v for k, v in mod20_data.items()}

with open("signature_22.txt", "r") as f:
    mod22_data = ast.literal_eval(f.read())
mod22_data = {int(k): v for k, v in mod22_data.items()}

with open("signature_24.txt", "r") as f:
    mod24_data = ast.literal_eval(f.read())
mod24_data = {int(k): v for k, v in mod24_data.items()}

with open("signature_26.txt", "r") as f:
    mod26_data = ast.literal_eval(f.read())
mod26_data = {int(k): v for k, v in mod26_data.items()}

with open("signature_28.txt", "r") as f:
    mod28_data = ast.literal_eval(f.read())
mod28_data = {int(k): v for k, v in mod28_data.items()}

with open("signature_30.txt", "r") as f:
    mod30_data = ast.literal_eval(f.read())
mod30_data = {int(k): v for k, v in mod30_data.items()}

with open("signature_34.txt", "r") as f:
    mod34_data = ast.literal_eval(f.read())
mod34_data = {int(k): v for k, v in mod34_data.items()}

with open("signature_36.txt", "r") as f:
    mod36_data = ast.literal_eval(f.read())
mod36_data = {int(k): v for k, v in mod36_data.items()}

with open("signature_37.txt", "r") as f:
    mod37_data = ast.literal_eval(f.read())
mod37_data = {int(k): v for k, v in mod37_data.items()}

with open("signature_38.txt", "r") as f:
    mod38_data = ast.literal_eval(f.read())
mod38_data = {int(k): v for k, v in mod38_data.items()}

with open("signature_40.txt", "r") as f:
    mod40_data = ast.literal_eval(f.read())
mod40_data = {int(k): v for k, v in mod40_data.items()}

with open("signature_15.txt", "r") as f:
    mod15_data = ast.literal_eval(f.read())
mod15_data = {int(k): v for k, v in mod15_data.items()}

with open("signature_210.txt", "r") as f:
    mod210_data = ast.literal_eval(f.read())
mod210_data = {int(k): v for k, v in mod210_data.items()}

with open("signature_6.txt", "r") as f:
    mod6_data = ast.literal_eval(f.read())
mod6_data = {int(k): v for k, v in mod6_data.items()}

with open("signature_8.txt", "r") as f:
    mod8_data = ast.literal_eval(f.read())
mod8_data = {int(k): v for k, v in mod8_data.items()}

with open("signature_10.txt", "r") as f:
    mod10_data = ast.literal_eval(f.read())
mod10_data = {int(k): v for k, v in mod10_data.items()}

with open("signature_12.txt", "r") as f:
    mod12_data = ast.literal_eval(f.read())
mod12_data = {int(k): v for k, v in mod12_data.items()}

with open("signature_6.txt", "r") as f:
    mod6_data = ast.literal_eval(f.read())
mod6_data = {int(k): v for k, v in mod6_data.items()}


with open("signature_11.txt", "r") as f:
    mod11_data = ast.literal_eval(f.read())
mod11_data = {int(k): v for k, v in mod11_data.items()}

with open("signature_14.txt", "r") as f:
    mod14_data = ast.literal_eval(f.read())
mod14_data = {int(k): v for k, v in mod14_data.items()}

with open("signature_32.txt", "r") as f:
    mod32_data = ast.literal_eval(f.read())
mod32_data = {int(k): v for k, v in mod32_data.items()}

with open("signature_34.txt", "r") as f:
    mod34_data = ast.literal_eval(f.read())
mod34_data = {int(k): v for k, v in mod34_data.items()}

with open("signature_42.txt", "r") as f:
    mod42_data = ast.literal_eval(f.read())
mod42_data = {int(k): v for k, v in mod42_data.items()}

with open("signature_9.txt", "r") as f:
    mod9_data = ast.literal_eval(f.read())
mod9_data = {int(k): v for k, v in mod9_data.items()}

with open("signature_48.txt", "r") as f:
    mod48_data = ast.literal_eval(f.read())
mod48_data = {int(k): v for k, v in mod48_data.items()}

with open("signature_19.txt", "r") as f:
    mod19_data = ast.literal_eval(f.read())
mod19_data = {int(k): v for k, v in mod19_data.items()}

with open("signature_26.txt", "r") as f:
    mod26_data = ast.literal_eval(f.read())
mod26_data = {int(k): v for k, v in mod26_data.items()}

with open("signature_46.txt", "r") as f:
    mod46_data = ast.literal_eval(f.read())
mod46_data = {int(k): v for k, v in mod46_data.items()}

with open("signature_60.txt", "r") as f:
    mod60_data = ast.literal_eval(f.read())
mod60_data = {int(k): v for k, v in mod60_data.items()}

with open("signature_100.txt", "r") as f:
    mod100_data = ast.literal_eval(f.read())
mod100_data = {int(k): v for k, v in mod100_data.items()}

with open("signature_5.txt", "r") as f:
    mod5_data = ast.literal_eval(f.read())
mod5_data = {int(k): v for k, v in mod5_data.items()}

with open("signature_70.txt", "r") as f:
    mod70_data = ast.literal_eval(f.read())
mod70_data = {int(k): v for k, v in mod70_data.items()}

# ----------------- Two-point factoring -----------------
# just to filter smaller prime factors. just a halfbaked failsafe

def b_candidates(num):
    # Get residues safely, fallback if missing
    # Stack the moduli one after another based on the size of N
    # dont stack too much prime moduli as it may lead to a CRT combination explosion
    # try to combine primorial and composite like 22,26,38
    # Then it is just a trial division with a huge computational overhead

    n_mod_3 = num % 3
    residues_3 = mod3_data.get(n_mod_3, list(range(3)))


    n_mod_4 = num % 4
    residues_4 = mod4_data.get(n_mod_4, list(range(4)))

    M_eff, R_eff = merge_residues(3, residues_3, 4, residues_4)

    n_mod_5 = num % 5
    residues_5 = mod5_data.get(n_mod_5, list(range(5)))

    M_eff, R_eff = merge_residues(M_eff, R_eff, 5, residues_5)

    n_mod_6 = num % 6
    residues_6 = mod6_data.get(n_mod_6, list(range(6)))

    M_eff, R_eff = merge_residues(M_eff, R_eff, 6, residues_6)

    n_mod_16 = num % 16
    residues_16 = mod16_data.get(n_mod_16, list(range(16)))

    M_eff, R_eff = merge_residues(M_eff, R_eff, 16, residues_16)

    n_mod_30 = num % 30
    residues_30 = mod30_data.get(n_mod_30, list(range(30)))

    M_eff, R_eff = merge_residues(M_eff, R_eff, 30, residues_30)


    n_mod_36 = num % 36
    residues_36 = mod36_data.get(n_mod_36, list(range(36)))

    M_eff, R_eff = merge_residues(M_eff, R_eff, 36, residues_36)


    n_mod_60 = num % 60
    residues_60 = mod60_data.get(n_mod_60, list(range(60)))

    M_eff, R_eff = merge_residues(M_eff, R_eff, 60, residues_60)



   

    n_mod_100 = num % 100
    residues_100 = mod100_data.get(n_mod_100, list(range(100)))

    M_eff, R_eff = merge_residues(M_eff, R_eff, 100, residues_100)

    n_mod_210 = num % 210
    residues_210 = mod210_data.get(n_mod_210, list(range(210)))

    M_eff, R_eff = merge_residues(M_eff, R_eff, 210, residues_210)

    n_mod_26 = num % 26
    residues_26 = mod26_data.get(n_mod_26, list(range(26)))

    M_eff, R_eff = merge_residues(M_eff, R_eff, 26, residues_26)

    n_mod_22 = num % 22
    residues_22 = mod22_data.get(n_mod_22, list(range(22)))

    M_eff, R_eff = merge_residues(M_eff, R_eff, 22, residues_22)

    n_mod_34 = num % 34
    residues_34 = mod34_data.get(n_mod_34, list(range(34)))

    M_eff, R_eff = merge_residues(M_eff, R_eff, 34, residues_34)

    n_mod_38 = num % 38
    residues_38 = mod38_data.get(n_mod_38, list(range(38)))

    M_eff, R_eff = merge_residues(M_eff, R_eff, 38, residues_38)

    return M_eff, R_eff




def two_point_factor(num):
    start_time = time.perf_counter()
    
    
    iterations = 0
    x = 3
    if num % 2 == 0:
        return None  # Early exit if num is even
    
    b_max = num // x
    sqrt_num = math.isqrt(num)
    b_min = 2 * sqrt_num
    
    

    SQUARE_RES_MOD64 = { (x*x) % 64 for x in range(64) }

    M_eff, R_eff = b_candidates(num)
    
    step_size = M_eff

    b_start = max(0, (b_min - min(R_eff)) // step_size)
    
    

    #return M_eff, len(R_eff), len(R_eff)/M_eff, end_time - start_time

    # ----------------Simple Main Loop ----------------
    i=0
    p=5
    step = 2
    while b_min <= b_max:
                
       
        base_b = (b_start + i) * step_size


        

        
        # accelerated fermat method using the modular constraint information
        for r in R_eff:
            

            #-------------trial division ----------------
            if num%p ==0:
                elapsed = time.perf_counter() - start_time
                return ("trial", p, num//p, elapsed, M_eff, len(R_eff), len(R_eff)/M_eff, iterations)
            p += step
            step = 6-step
            
            b = base_b + r
            d = b * b - 4 * num  # b^2 - 4*num

            if d<0: continue
            
            
            iterations += 1
            if (d & 63) not in SQUARE_RES_MOD64: continue
                
            root = math.isqrt(d)  # Integer square root of d
            if root * root == d:  # Check if d is a perfect square
                factor1 = (b - root) // 2
                factor2 = (b + root) // 2
                elapsed = time.perf_counter() - start_time
                return ("b_candidate", factor1, factor2, elapsed, M_eff, len(R_eff), len(R_eff)/M_eff, iterations)
        
            
        i+=1
        if (b_start +i) * step_size > num//3:
            break

    return None



# ----------------- Generate semiprimes -----------------
primes = [sympy.randprime(10**9, 10**11) for _ in range(50)]
semiprimes = [random.choice(primes) * random.choice(primes) for _ in range(50)]
semiprimes = []
for _ in range(5):
    while True:
        p = sympy.randprime(10**9, 10**10)
        q = sympy.randprime(10**9, 10**10)
        if p != q:
            semiprimes.append(p * q)
            break

#----------------- Factorization and table -----------------
results = []
print("--- Starting Factorization ---")
for N in semiprimes:
    res = two_point_factor(N)
    
    if res:
        point, f1, f2, elapsed, M, R , r_m, iters = res
        results.append([N, point, f1, f2, f"{elapsed:.6f}", M, R, f"{r_m: .6f}", iters])
    else:
        results.append([N, "None", "-", "-", "-", "-"])

# Just a table to see the results
# ----------------- Print formatted table -----------------
headers = ["N", "point", "Factor1", "Factor2", "Time(s)","Comb_M", "R","R/M", "Iterations"]
print(tabulate(results, headers=headers, tablefmt="grid"))

import csv
from datetime import datetime


# Save results to CSV with unique filename
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
filename = f"factorization_results_{timestamp}.csv"

with open(filename, mode="w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(headers)
    writer.writerows(results)

print(f"\nResults saved to {filename}")
