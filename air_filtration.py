import pandas as pd
import numpy as np
from scipy.optimize import minimize

df = pd.read_csv("air_filtration.csv")

properties = df[["anisotropy", "density", "band gap", "bulk modulus", "shear modulus"]].values

properties_mean = np.mean(properties, axis=0)
properties_std = np.std(properties, axis=0)
properties_normalized = (properties - properties_mean) / properties_std

required_properties = np.array([0.0, 4.0, 5.0, 200, 70])
required_properties_normalized = (required_properties - properties_mean) / properties_std

weights = np.array([2.0, 1.0, 1.0, 1.0, 1.0])  
weighted_diff = (properties_normalized - required_properties_normalized) * weights
distances = np.linalg.norm(weighted_diff, axis=1) 

closest_indices = np.argsort(distances)[:5]
closest_materials = df.iloc[closest_indices]

print("Closest materials:")
print(closest_materials)

def composite_properties(weights, properties):
    return np.sum(weights[:, np.newaxis] * properties, axis=0)

def objective_function(weights):
    composite = composite_properties(weights, closest_materials[["anisotropy", "density", "band gap", "bulk modulus", "shear modulus"]].values)
    return np.sum(((composite - required_properties_normalized) * weights) ** 2)

constraints = [{'type': 'eq', 'fun': lambda w: np.sum(w) - 1}]
bounds = [(0, 1) for _ in range(len(closest_materials))]

initial_weights = np.ones(len(closest_materials)) / len(closest_materials)

result = minimize(objective_function, initial_weights, bounds=bounds, constraints=constraints)
optimized_weights = result.x

print("\nOptimized weights for composite:")
for i, material in enumerate(closest_materials['chemical formula'].values):
    print(f"{material}: {optimized_weights[i]:.4f}")

final_properties = composite_properties(optimized_weights, closest_materials[["anisotropy", "density", "band gap", "bulk modulus", "shear modulus"]].values)


final_properties_original = final_properties * properties_std + properties_mean  

print(f"\nFinal composite properties:")
print(f'Anisotropy: {final_properties_original[0]:.4f}')
print(f'Density: {final_properties_original[1]:.4f}')
print(f'Band Gap: {final_properties_original[2]:.4f}')
print(f'Bulk Modulus: {final_properties_original[3]:.4f}')
print(f'Shear Modulus: {final_properties_original[4]:.4f}')
