import pandas as pd
import numpy as np
from scipy.optimize import minimize

df = pd.read_csv("heat_shield.csv")

properties = df[["energy above hull", "formation energy", "band gap"]].values

properties_mean = np.mean(properties, axis=0)
properties_std = np.std(properties, axis=0)
properties_normalized = (properties - properties_mean) / properties_std

required_properties = np.array([0.0, -2.5, 2.3])
required_properties_normalized = (required_properties - properties_mean) / properties_std

weights = np.array([1.0, 1.0, 1.0])  # Higher weight for band gap & anisotropy
weighted_diff = (properties_normalized - required_properties_normalized) * weights
distances = np.linalg.norm(weighted_diff, axis=1)  # Weighted Euclidean distance

closest_indices = np.argsort(distances)[:5]
closest_materials = df.iloc[closest_indices]

print("Closest materials:")
print(closest_materials)

def composite_properties(weights, properties):
    return np.sum(weights[:, np.newaxis] * properties, axis=0)

def objective_function(weights):
    composite = composite_properties(weights, closest_materials[["energy above hull", "formation energy", "band gap"]].values)
    return np.sum(((composite - required_properties_normalized) * weights) ** 2) 

constraints = [{'type': 'eq', 'fun': lambda w: np.sum(w) - 1}]
bounds = [(0, 1) for _ in range(len(closest_materials))]

initial_weights = np.ones(len(closest_materials)) / len(closest_materials)

result = minimize(objective_function, initial_weights, bounds=bounds, constraints=constraints)
optimized_weights = result.x

print("\nOptimized weights for composite:")
for i, material in enumerate(closest_materials['chemical formula'].values):
    print(f"{material}: {optimized_weights[i]:.4f}")

final_properties = composite_properties(optimized_weights, closest_materials[["energy above hull", "formation energy", "band gap"]].values)

final_properties_original = final_properties * properties_std + properties_mean  

print(f"\nFinal composite properties:")
print(f'Energy Above Hull: {final_properties_original[0]:.4f}')
print(f'Formation Energy: {final_properties_original[1]:.4f}')
print(f'Band Gap: {final_properties_original[2]:.4f}')
