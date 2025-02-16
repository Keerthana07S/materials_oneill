Keerthana_API = "insert api here"

import pandas as pd 
import numpy as np 
import pymatgen.core as mat 
from pymatgen.ext.matproj import MPRester

mpr = MPRester(Keerthana_API)

def orbital_habitat_hulls():
    docs = mpr.materials.summary.search(
        is_metal=True,
        fields=["bulk_modulus", "shear_modulus", "homogeneous_poisson", "formula", "density"]
    )
    
    data = []
    for doc in docs:
        bulk_mod = doc.bulk_modulus
        shear_mod = doc.shear_modulus
        poisson = doc.homogeneous_poisson
        chemical_formula = doc.formula
        density = doc.density
        
        if (bulk_mod is None or shear_mod is None or 
            bulk_mod == 0 or shear_mod == 0 or 
            bulk_mod['vrh'] < 0 or shear_mod['vrh'] < 0):
            print(f"Invalid or zero bulk_modulus or shear_modulus for material {doc.material_id}. Skipping this material.")
            continue
        else:
            bulk_mod = bulk_mod['vrh']
            shear_mod = shear_mod['vrh']

        # Ensure Poisson’s ratio is not None
        if poisson is None:
            print(f"Missing Poisson’s ratio for material {doc.material_id}. Skipping this material.")
            continue

        # Calculate the material properties
        try:
            E = (3 * bulk_mod * shear_mod) / (bulk_mod + 2 * shear_mod)
            tensile = 0.25 * E
            yield_strength = (2 * shear_mod) / (1 + poisson)  # Poisson is safely added by 1
            fracture_toughness = tensile * (np.sqrt(E))
            thermal_expansion_coefficient = shear_mod / bulk_mod

            # Add material properties to the data list
            data.append({
                "chemical_formula": chemical_formula,
                "tensile_strength": tensile,
                "yield_strength": yield_strength,
                "fracture_toughness": fracture_toughness,
                "elastic_modulus": E,
                "thermal_expansion_coefficient": thermal_expansion_coefficient,
                "density": density
            })
        except ZeroDivisionError as e:
            print(f"Error in calculation for material {doc.material_id}: {e}. Skipping this material.")
            continue

    # Convert to DataFrame and save as CSV
    if data:
        df = pd.DataFrame(data)
        print(df)
        df.to_csv('orbital_habitat_hulls.csv', index=False)
    else:
        print("No valid data to write to CSV.")

def pressurized_modules():
    docs = mpr.materials.summary.search(
        band_gap=[3.0, 1000.0],
        fields=["bulk_modulus", "shear_modulus", "energy_above_hull", "formation_energy_per_atom", "formula_pretty", "band_gap", "density"]
    )
    
    data = []
    
    for doc in docs:
        bulk_mod = doc.bulk_modulus
        shear = doc.shear_modulus
        eabovehull = doc.energy_above_hull
        formation_energy = doc.formation_energy_per_atom
        formula = doc.formula_pretty
        band_gap = doc.band_gap
        density = doc.density
        
        if (bulk_mod is None or shear is None or 
            bulk_mod == 0 or shear == 0 or 
            bulk_mod['vrh'] < 0 or shear['vrh'] < 0):
            print(f"Invalid or zero bulk_modulus or shear_modulus for material {doc.material_id}. Skipping this material.")
            continue
        else:
            bulk_mod = bulk_mod['vrh']
            shear = shear['vrh']
        
        if (eabovehull is None or formation_energy is None or formula is None or band_gap is None or density is None):
            print("Invalid material. Skipping this material")
            continue

        data.append({
            "chemical formula": formula, 
            "bulk modulus":bulk_mod,
            "shear modulus":shear,
            "energy above hull": eabovehull, 
            "formation energy": formation_energy, 
            "band gap": band_gap, 
            "density": density
        })
        
    if data:
        df = pd.DataFrame(data)
        print(df)
        df.to_csv('pressurized_modules.csv', index=False)
    else:
        print("No valid data to write to CSV.")

def solar_panel_materials():
    docs = mpr.materials.summary.search(
        band_gap = [1.1, 1.5], 
        fields = ["formula_pretty", "energy_above_hull", "is_stable", "band_gap", "energy_per_atom", "density", "e_electronic", "weighted_surface_energy"]
    )
    
    data = []
    
    for doc in docs:
        formula = doc.formula_pretty
        eabovehull = doc.energy_above_hull
        is_stable = doc.is_stable
        bandgap = doc.band_gap
        energyperatom = doc.energy_per_atom
        density = doc.density
        e_electronic = doc.e_electronic
        weighted_surface_energy = doc.weighted_surface_energy

        if (is_stable is None or formula is None or eabovehull is None or bandgap is None or energyperatom is None or density is None):
            print("Invalid material. Skipping this material")
            continue
        else:
            data.append({
                "chemical formula": formula, 
                "energy above hull": eabovehull, 
                "is stable": is_stable, 
                "bandgap": bandgap, 
                "energyperatom": energyperatom, 
                "density": density 
            })
    if data:
        df = pd.DataFrame(data)
        print(df)
        df.to_csv('solar_panel_materials.csv', index=False)
    else:
        print("No valid data to write to CSV.")
    
def battery_storage():
    docs = mpr.materials.summary.search(fields = ["formula_pretty", "density", "energy_per_atom", "formation_energy_per_atom", "bulk_modulus", "shear_modulus", "is_stable", "energy_above_hull"])
    data = []
    
    for doc in docs:
        formula = doc.formula_pretty
        density = doc.density
        energyperatom = doc.energy_per_atom
        formatione = doc.formation_energy_per_atom 
        bulkmod = doc.bulk_modulus
        shearmod = doc.shear_modulus
        stability = doc.is_stable
        eabovehull = doc.energy_above_hull
        
        if (bulkmod is None or shearmod is None or 
            bulkmod == 0 or shearmod == 0 or 
            bulkmod['vrh'] < 0 or shearmod['vrh'] < 0):
            print(f"Invalid or zero bulk_modulus or shear_modulus for material {doc.material_id}. Skipping this material.")
            continue
        else:
            bulk_mod = bulkmod['vrh']
            shear = shearmod['vrh']
            
        if (formula is None or density is None or energyperatom is None or formatione is None or bulkmod is None or shearmod is None or stability is None or eabovehull is None):
            print("Invalid material. Skipping this material")
            continue
        else:
            data.append({
                "chemical formula": formula, 
                "energy above hull": eabovehull, 
                "density": density,
                "energy per atom": energyperatom,
                "formation energy per atom":formatione,
                "bulk modulus": bulk_mod,
                "shear modulus": shear,
                "stability": stability
            })
    if data:
        df = pd.DataFrame(data)
        print(df)
        df.to_csv('battery_storage.csv', index=False)
    else:
        print("No valid data to write to CSV.")
        
        
def primary_shielding():
    docs = mpr.materials.summary.search(elements = ["C", "H"], fields = ["formula_pretty", "density"])
    data = []
    for doc in docs:
        formula = doc.formula_pretty
        density = doc.density
        
        if (density is None or formula is None):
            print("Invalid material. Skipping this material")
            continue
        else:
            data.append({
                "chemical formula": formula, 
                "density": density
            })
    if data:
        df = pd.DataFrame(data)
        print(df)
        df.to_csv('primary_shielding.csv', index=False)
    else:
        print("No valid data to write to CSV.")    
        
def secondary_shielding():
    docs = mpr.materials.summary.search(is_metal=True, 
        elements = ["H", "C"], fields = ["density", "formula_pretty", "composition"])
    data = []
    for doc in docs:
        formula = doc.formula_pretty 
        density = doc.density
        if (density is None or formula is None):
            print("Invalid material. Skipping this material")
            continue
        else:
            data.append({
                "chemical formula": formula, 
                "density": density
            })
    if data:
        df = pd.DataFrame(data)
        print(df)
        df.to_csv('secondary_shielding.csv', index=False)
    else:
        print("No valid data to write to CSV.") 
    
    
def water_filtration_purification():
    docs = mpr.materials.summary.search(elements = ["H", "C"], fields = ["formula_pretty", "density", "universal_anisotropy", "is_stable"])
    data = []
    
    for doc in docs:
        formula = doc.formula_pretty
        density = doc.density 
        anisotropy = doc.universal_anisotropy 
        stability = doc.is_stable
        
        if (formula is None or density is None or anisotropy is None or stability is None):
            print("Invalid material. Skipping this material")
            continue 
        else:
            data.append({
                "chemical formula":formula, 
                "density": density,  
                "anisotropy": anisotropy, 
                "stability": stability
            })
    if data:
        df = pd.DataFrame(data)
        print(df)
        df.to_csv('water_filtration_purification.csv', index=False)
    else:
        print("No valid data to write to CSV.") 

def air_filtration():
    docs = mpr.materials.summary.search(elements = ["O"], fields = ["formula_pretty", "universal_anisotropy", "density", "band_gap", "bulk_modulus", "shear_modulus"])
    data =[]
    for doc in docs:
        formula = doc.formula_pretty, 
        anisotropy = doc.universal_anisotropy,  
        density = doc.density, 
        band_gap = doc.band_gap, 
        bulkmod = doc.bulk_modulus, 
        shearmod = doc.shear_modulus

        if (bulkmod is None or shearmod is None or 
                bulkmod == 0 or shearmod == 0):
                print(f"Invalid or zero bulk_modulus or shear_modulus for material {doc.material_id}. Skipping this material.")
                continue
        else:
            bulk_mod = bulkmod[0]['vrh']
            shear = shearmod['vrh']
        if (formula is None or anisotropy is None  or density is None or band_gap is None):
            print("Invalid material. Skipping this material")
            continue
        else: 
            data.append({
                "chemical formula":formula, 
                "anisotropy": anisotropy,  
                "density": density, 
                "band gap": band_gap, 
                "bulk modulus": bulk_mod, 
                "shear modulus": shear
            })
    if data:
        df = pd.DataFrame(data)
        print(df)
        df.to_csv('air_filtration.csv', index=False)
    else:
        print("No valid data to write to CSV.") 
    
def biological_growth():
    docs = mpr.materials.summary.search(band_gap = [3.1, 1000], fields = ["band_gap", "formula_pretty", "density", "universal_anisotropy"])
    data = []
    for doc in docs:
        formula = doc.formula_pretty 
        band_gap = doc.band_gap
        density = doc.density 
        anisotropy = doc.universal_anisotropy 
        
        if (formula is None or band_gap is None or density is None or anisotropy is None):
            print("Invalid material. Skipping this material")
            continue
        else:
            data.append({
                "chemical formula":formula, 
                "band gap": band_gap, 
                "density": density, 
                "anisotropy": anisotropy
            })
    if data:
        df = pd.DataFrame(data)
        print(df)
        df.to_csv('biological_growth.csv', index=False)
    else:
        print("No valid data to write to CSV.")    

def orbital_tug():
    docs = mpr.materials.summary.search(elements=["W"], fields=["formula_pretty", "is_stable", "shear_modulus", "bulk_modulus"])
    data = []
    
    for doc in docs:
        formula = doc.formula_pretty, 
        stability = doc.is_stable, 
        shearmod = doc.shear_modulus, 
        bulkmod = doc.bulk_modulus
        
        if (bulkmod is None or shearmod is None or 
                bulkmod == 0 or shearmod == 0):
                print(f"Invalid or zero bulk_modulus or shear_modulus for material {doc.material_id}. Skipping this material.")
                continue
        else:
            bulk_mod = bulkmod['vrh']
            shear_mod = shearmod[0]['vrh']
        
        if(formula is None or  stability is None):
            print("Invalid material. Skipping this material")
            continue
        else:
            data.append({
                "chemical formula":formula, 
                "stability": stability, 
                "shear modulus":shear_mod,
                "bulk modulus": bulk_mod
            })
    if data:
        df = pd.DataFrame(data)
        print(df)
        df.to_csv('orbital_tug.csv', index=False)
    else:
        print("No valid data to write to CSV.")   
            
def heat_shield():
    docs = mpr.materials.summary.search(band_gap = [0.5, 3.0], fields = ["formula_pretty", "energy_above_hull", "formation_energy_per_atom", "band_gap"])
    data=[]
    for doc in docs:
        formula = doc.formula_pretty, 
        energy_above_hull = doc.energy_above_hull
        formation_energy = doc.formation_energy_per_atom, 
        band_gap = doc.band_gap
        
        if (formula is None or energy_above_hull is None or formation_energy is None or band_gap is None):
            print("Invalid material. Skipping this material.")
        else:
            data.append({
                "chemical formula":formula, 
                "energy above hull":energy_above_hull, 
                "formation energy":formation_energy,
                "band gap": band_gap 
            })
    if data:
        df = pd.DataFrame(data)
        print(df)
        df.to_csv('heat_shield1.csv', index=False)
    else:
        print("No valid data to write to CSV.")   
    
def radiative_cooling():
    docs = mpr.materials.summary.search(elements = ["C", "H"], fields=["formula_pretty", "band_gap", "is_stable"])
    data=[]
    for doc in docs:
        formula = doc.formula_pretty, 
        band_gap = doc.band_gap, 
        stability = doc.is_stable
        
        if (formula is None or band_gap is None or stability is None):
            print("Invalid material. Skipping this material")
        else:
            data.append({
                "chemical formula":formula, 
                "band gap": band_gap, 
                "stability": stability
            })
    if data:
        df = pd.DataFrame(data)
        print(df)
        df.to_csv('radiative_cooling.csv', index=False)
    else:
        print("No valid data to write to CSV.")   
