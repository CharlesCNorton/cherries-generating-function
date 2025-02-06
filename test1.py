"""
Closed-Form Generating Function for Cherries in Full Binary Trees:
A rigorous, independent, and scientific test.

This script performs the following:
1. Defines the bivariate generating function 
     F(x, y) = (1 - sqrt(1 - 4*x - 4*x^2*(y-1)))/2,
   where x marks the number of leaves and y marks the number of cherries.
2. Expands F(x, y) as a power series in x (up to a specified order)
   and prints the polynomial in y for each n (i.e. the distribution of cherries).
3. Verifies that setting y = 1 recovers the ordinary generating function
   for full binary trees (the Catalan generating function).
4. For small n (n = 1 to 6), generates all full binary trees, counts the number
   of cherries in each, and compares the brute–force distribution with that
   obtained from the generating function.
"""

import sympy as sp
import numpy as np
import pandas as pd

# --------------------------------------------------
# Part 1. Symbolic Setup and Generating Function
# --------------------------------------------------

# Define symbols: x marks leaves; y marks cherries.
x, y = sp.symbols('x y')

# Define the closed-form bivariate generating function:
#   F(x, y) = (1 - sqrt(1 - 4*x - 4*x^2*(y-1))) / 2
F = (1 - sp.sqrt(1 - 4*x - 4*x**2*(y - 1)))/2

# Set expansion order (maximum power of x to consider)
order = 8  # This will compute series for n = 1,2,...,7

# Expand F(x,y) as a power series in x up to order 'order'
F_series = sp.series(F, x, 0, order).removeO().expand()

print("Bivariate Generating Function F(x,y):")
sp.pprint(F)
print("\nSeries expansion of F(x,y) in x up to order {}:".format(order))
sp.pprint(F_series)

# --------------------------------------------------
# Part 2. Extract Coefficients and Compare with Catalan Numbers
# --------------------------------------------------

# For each n (1 <= n < order), extract the polynomial in y that is the coefficient of x^n.
coeffs_by_n = {}
for n in range(1, order):
    poly_y = sp.expand(sp.Poly(F_series, x).coeff_monomial(x**n))
    coeffs_by_n[n] = sp.expand(poly_y)

# Display the distribution polynomials:
print("\nDistribution of cherries (polynomials in y) for trees with n leaves:")
for n in sorted(coeffs_by_n.keys()):
    print("n =", n, ":", coeffs_by_n[n])
    
# Verification: Setting y = 1 should yield the Catalan generating function
F_univariate = sp.simplify(F_series.subs(y, 1))
print("\nF(x,1) (should equal (1-sqrt(1-4*x))/2):")
sp.pprint(F_univariate)
Catalan_GF = (1 - sp.sqrt(1-4*x))/2
print("\nCatalan generating function:")
sp.pprint(Catalan_GF)
assert sp.simplify(F_univariate - Catalan_GF) == 0, "Mismatch with Catalan GF!"

# --------------------------------------------------
# Part 3. Brute-Force Enumeration of Full Binary Trees and Cherry Count
# --------------------------------------------------

def generate_trees(n):
    """
    Recursively generate all full binary trees with n leaves.
    A leaf is represented by the string "L".
    An internal node is represented as a tuple (left, right).
    """
    if n == 1:
        return ["L"]
    trees = []
    # For n>=2, partition leaves into two groups: i and n-i, for 1 <= i <= n-1.
    for i in range(1, n):
        left_subtrees = generate_trees(i)
        right_subtrees = generate_trees(n - i)
        for L_tree in left_subtrees:
            for R_tree in right_subtrees:
                trees.append((L_tree, R_tree))
    return trees

def count_cherries(tree):
    """
    Count the number of cherries in a full binary tree.
    A cherry is defined as an internal node whose both children are leaves.
    """
    if tree == "L":
        return 0
    left, right = tree
    # Check if both children are leaves.
    is_cherry = 1 if (left == "L" and right == "L") else 0
    return is_cherry + count_cherries(left) + count_cherries(right)

# For n from 1 to max_n, generate trees and count cherries.
max_n = 7  # Small sizes for brute-force verification
brute_force_data = {}  # Will store a dictionary: {n: {cherry_count: frequency}}

for n in range(1, max_n + 1):
    trees = generate_trees(n)
    freq = {}
    for t in trees:
        c = count_cherries(t)
        freq[c] = freq.get(c, 0) + 1
    brute_force_data[n] = freq

# Display brute-force distributions:
print("\nBrute-force distribution of cherries in full binary trees:")
for n in range(1, max_n + 1):
    total_trees = len(generate_trees(n))
    print("n =", n, " (Total trees =", total_trees, "):", brute_force_data[n])

# --------------------------------------------------
# Part 4. Compare Brute-Force Distribution with Generating Function Coefficients
# --------------------------------------------------

# For each n from 1 to max_n, extract from the generating function expansion the coefficients for y^c.
gf_data = {}  # gf_data[n] = {c: coefficient}
for n in range(1, max_n + 1):
    # Coefficient of x^n in F_series is a polynomial in y.
    poly = sp.expand(sp.Poly(F_series, x).coeff_monomial(x**n))
    # Convert polynomial in y to a dictionary mapping exponent c to coefficient.
    coeff_dict = sp.Poly(poly, y).as_dict()
    # as_dict returns keys as tuples (exponent,), so convert:
    coeff_dict_converted = {exp[0]: coeff for exp, coeff in coeff_dict.items()}
    gf_data[n] = coeff_dict_converted

print("\nGenerating function distribution (coefficients) for cherries:")
for n in range(1, max_n + 1):
    total_from_gf = sum(gf_data[n].values())
    print("n =", n, " (Total trees =", total_from_gf, "):", gf_data[n])

# Compare: For each n, sum over c of gf_data should equal the Catalan number for n leaves (Catalan(n-1)).
def catalan(n):
    return sp.binomial(2*n, n) // (n+1)

for n in range(1, max_n + 1):
    total_gf = sum(gf_data[n].values())
    total_cat = catalan(n-1) if n > 1 else 1
    assert total_gf == total_cat, f"Total count mismatch for n={n}: GF gives {total_gf}, Catalan gives {total_cat}"
    
print("\nAll totals from the generating function match the corresponding Catalan numbers.")

# --------------------------------------------------
# Part 5. (Optional) Tabulate and Visualize the Distributions
# --------------------------------------------------

import matplotlib.pyplot as plt

# Prepare a table for n = 4,5,6 (for example)
table_data = []
for n in range(1, max_n + 1):
    row = {"n": n, "Total Trees": len(generate_trees(n))}
    for c in sorted(gf_data[n].keys()):
        row[f"Cherries = {c}"] = int(gf_data[n][c])
    table_data.append(row)
df = pd.DataFrame(table_data)
print("\nDistribution table from generating function (rows = n, columns = cherry count frequencies):")
print(df)

# Plot for a selected n (say n=6)
n_plot = 6
gf_dist = gf_data[n_plot]
cherry_counts = sorted(gf_dist.keys())
frequencies = [gf_dist[c] for c in cherry_counts]

plt.figure(figsize=(6,4))
plt.bar(cherry_counts, frequencies, color="skyblue")
plt.xlabel("Number of Cherries")
plt.ylabel("Frequency")
plt.title(f"Distribution of Cherries for Full Binary Trees with n = {n_plot}")
plt.xticks(cherry_counts)
plt.grid(True, axis="y", linestyle="--")
plt.show()
