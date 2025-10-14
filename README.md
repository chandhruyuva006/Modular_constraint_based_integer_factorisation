# Modular Patterns in Semiprimes: Deterministic Congruence Framework for Factorization

##  Overview
This repository accompanies the paper  
**“Modular Patterns in Semiprimes: Empirical Observations and Conjectures on Factor Sums and Deterministic Congruence Patterns in Semiprimes Enabling Optimised Factorization.”**

It presents a **deterministic, modular-arithmetic-based approach** to integer factorization that identifies structured, non-random patterns in semiprime residue spaces.  
The central conjecture — the **Semiprime Modular Restriction Conjecture (SMRC)** — posits that the sum of a semiprime’s factors `S_N = p + q` obeys a restricted set of modular constraints that can be deterministically reconstructed using combinations of modular congruences and the **Chinese Remainder Theorem (CRT)**.

---

##  Implementation
The entire framework is implemented in **pure Python (standard library only)** — no external dependencies are required.
Once the theory (just a simple modular math with a matrix) is understood, then the implementation is pretty straightforward. 

### Files
#### use the factoring_using_stacking_modular_signature.py file and 
#### mod_signature_synthetic.ipynb 
i have described the theory as much i can possible in the comments.


---

##  Workflow Summary
### **Step 1: Precompute Modular Signatures**
This step constructs and stores residue–constraint mappings for each modulus `m`:
Using the mod_signature_synthetic.py if you need to tweak with other moduli.
Moduli upto 510510 was tested. but the file sizes are huge. that is the only concern
I am aware of the better strategies are there to reduce the file size and even 
on the fly signature generation for small moduli. but considering the file sizes 
I think we can safely use 2310 signature. the primorial 30030 file is 60 MB.
But stacking the various smaller m leads to M_eff upto 177657480 in several cases.


### **Step 2: Factorisation**

### Use the factoring_using_stacking_modular_signature.py file in the Factoring folder
##### can use other files in the repo but need to edit a few logics and the files dont communicate 
##### to each other except the txt
##### files and the factoring_using_stacking_modular_signature.py
##### in a Dell Inspiron 3558 i3 5005U upto 18 digit odd semiprimes can be factored under 
##### a second in this crude implementation and 
feel free to suggest some improvements. 
And all suggestions are welcomed.

# Disclaimer : It is a crude implementation and no consideration was given to
# optimal run time, this implementation is used only to see if there any flaw
# in the theory and counter example.
# due to the assumptions like the prime factors are not divisible by 5 and we have
# fairly close primes to create the semiprime as we expect in the fermat style algorithms.
# so a tight bound was used like 10^5 to 10^10 primes only to construct the semiprimes.
# the program should not run if we put a composite with more than 3 factors. 
# simply work with moderate length semiprimes N = p*q only. 
   
