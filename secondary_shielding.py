import pandas as pd
import numpy as np
from scipy.optimize import minimize

df = pd.read_csv("primary_shielding.csv")

densities = df['density'].values

densities_min = np.min(densities)
densities_max = np.max(densities)
densities_normalized = (densities - densities_min) / (densities_max - densities_min)

required_density = 1.5  # Example (in g/cm^3)
required_density_normalized = (required_density - densities_min) / (densities_max - densities_min)

distances = np.abs(densities_normalized - required_density_normalized)

closest_indices = np.argsort(distances)[:5]

closest_materials = df.iloc[closest_indices]
print("Closest materials:")
print(closest_materials)

def composite_density(weights, densities):
    return np.sum(weights * densities)

def objective_function(weights):
    return (composite_density(weights, closest_materials['density'].values) - required_density) ** 2

constraints = [{'type': 'eq', 'fun': lambda w: np.sum(w) - 1}]

bounds = [(0, 1) for _ in range(len(closest_materials))]

initial_weights = [1 / len(closest_materials)] * len(closest_materials)

result = minimize(objective_function, initial_weights, bounds=bounds, constraints=constraints)

optimized_weights = result.x
print("\nOptimized weights for composite:")
for i, material in enumerate(closest_materials['chemical formula'].values):
    print(f"{material}: {optimized_weights[i]:.4f}")

final_density = composite_density(optimized_weights, closest_materials['density'].values)
print(f"\nFinal composite density: {final_density:.4f} g/cm^3")
