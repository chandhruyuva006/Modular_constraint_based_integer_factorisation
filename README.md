# Modular Patterns in Semiprimes: Deterministic Congruence Framework for Factorization

## 📘 Overview
This repository accompanies the paper  
**“Modular Patterns in Semiprimes: Empirical Observations and Conjectures on Factor Sums and Deterministic Congruence Patterns in Semiprimes Enabling Optimised Factorization.”**

It presents a **deterministic, modular-arithmetic-based approach** to integer factorization that identifies structured, non-random patterns in semiprime residue spaces.  
The central conjecture — the **Semiprime Modular Restriction Conjecture (SMRC)** — posits that the sum of a semiprime’s factors `S_N = p + q` obeys a restricted set of modular constraints that can be deterministically reconstructed using combinations of modular congruences and the **Chinese Remainder Theorem (CRT)**.

---

## 🧠 Theoretical Foundation
For any odd semiprime \( N = p \times q \):

\[
p + q \equiv r_m \pmod{m}, \quad \text{where } r_m \in R_m
\]

- \( R_m \) is a **restricted residue set**, deterministically derivable from \( N \mod m \).  
- By combining multiple moduli \( m_1, m_2, \dots, m_k \), one constructs a **system of congruences** that uniquely encodes the true factor sum \( S_N \).  
- Once \( S_N \) is identified, the factors are recovered via the quadratic equation:

\[
x^2 - S_Nx + N = 0
\]

The system theoretically achieves completeness when the combined CRT modulus exceeds \( N \), guaranteeing unique factor recovery.

---

## ⚙️ Implementation
The entire framework is implemented in **pure Python (standard library only)** — no external dependencies are required.

### Core Modules
| File | Description |
|------|--------------|
| `signature_generator.py` | Generates modular signatures \( R_m \) for chosen moduli |
| `factorization_engine.py` | Performs factorization using precomputed modular signatures and CRT combination |
| `utils.py` | Contains helper functions for modular arithmetic, CRT, and result formatting |
| `data/signatures/` | Directory storing precomputed modular signature JSON files |
| `results/` | Contains empirical test results and runtime data |

---

## 🧩 Workflow Summary
### **Step 1: Precompute Modular Signatures**
This step constructs and stores residue–constraint mappings for each modulus `m`:

```bash
python signature_generator.py
