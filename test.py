import numpy as np

# Set a random seed for reproducibility
np.random.seed(123)

# Classical optimal solution (single point)
classical_opt = 1.2e6  # e.g., a known optimal solution near 1e6

# Generate "All hybrid solutions" data:
# We'll create two clusters of solutions for illustrative purposes.
# One cluster around ~3e6 and another around ~1.2e7.
n_per_cluster = 5000
cluster1 = np.random.normal(loc=3e6, scale=5e5, size=n_per_cluster)   # Cluster centered at 3e6
cluster2 = np.random.normal(loc=1.2e7, scale=5e5, size=n_per_cluster) # Cluster centered at 12e6
all_hybrid_solutions = np.concatenate([cluster1, cluster2])

# "Feasible solutions":
# Let's define "feasible" arbitrarily as those solutions that have an objective 
# value less than, say, 5e6. In a real scenario, feasibility criteria should 
# reflect actual constraints or solution qualities.
feasible_solutions = all_hybrid_solutions[all_hybrid_solutions < 5e6]

# Now we have:
# classical_opt - a single float value
# all_hybrid_solutions - a large array of all generated solutions
# feasible_solutions - a subset of the above

# Example of how you might plot:

import matplotlib.pyplot as plt

plt.figure(figsize=(8,6))
plt.hist(all_hybrid_solutions, bins=50, color='blue', alpha=0.7, label='All hybrid solutions')
plt.hist(feasible_solutions, bins=50, color='orange', alpha=0.7, label='Feasible solutions')
plt.axvline(classical_opt, color='red', linestyle='--', label='Classical optimal solution')
plt.xscale('log')  # If you want a log scale on x; or adjust depending on your original plot
plt.yscale('log')  # If you want a log scale on y as in the provided figure
plt.xlabel('Objective function value')
plt.ylabel('Number of solutions')
plt.legend()
plt.tight_layout()
plt.show()

