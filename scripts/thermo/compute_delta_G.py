from equilibrator_api import ComponentContribution, Q_
from rdkit import Chem
import math
from data import COMPOUNDS_PATH, compounds
import numpy
import json 
import os

COMPOUNDS_STRUCTURE_PATH = COMPOUNDS_PATH+"/{compound}/structure/{compound}_{pubchem_id}.sdf"
COMPOUNDS_QUANTUM_PATH = COMPOUNDS_PATH+"/{compound}/structure/{compound}_{pubchem_id}.sdf"
COMPOUNDS_THERMO_PATH = COMPOUNDS_PATH+"/{compound}/thermo/{compound}_{pubchem_id}.json"

# cc = ComponentContribution()

# # set environment
# cc.p_h = Q_(7.4, "")          # unitless (pH)
# cc.p_mg = Q_(0.0, "")        # set if you want Mg corrections
# cc.ionic_strength = Q_(0.25, "M")
# cc.temperature = Q_(298.15, "K")

# # https://equilibrator.readthedocs.io/en/latest/equilibrator_examples.html

# def get_cpd(mol):
#     try:
#         inchi = Chem.MolToInchi(mol)
        
#         cpd = cc.get_compound_by_inchi(inchi)
        
#         print(cpd)
        
#         if cpd:
#             return cpd
#     except Exception:
#         print("Could not create equilibrator compound from InChI")

#     # raise ValueError("Could not create equilibrator compound from InChI")


# def dgf_prime_from_sdf(path):
#     """Return (mu_prime_kJmol, sigma_kJmol, compound_object) for the SDF file."""
#     mol = Chem.MolFromMolFile(path)
#     if mol is None:
#         raise ValueError(f"Could not parse SDF: {path}")

#     cpd = get_cpd(mol)

#     # standard formation mean (mu) and uncertainties returned by the API
#     mu, sigma_fin, sigma_inf = cc.standard_dg_formation(cpd)   # probably in kJ/mol
#     # apply pH / ionic strength / temp transformation; result as quantity
#     delta_q = cpd.transform(cc.p_h, cc.ionic_strength, cc.temperature, cc.p_mg)
#     delta_kJ = delta_q.m_as("kJ/mol")

#     mu_prime = mu + delta_kJ
#     delta_kJ_cov = sigma_fin @ sigma_fin.T + 1e6 * sigma_inf @ sigma_inf.T

#     # choose an uncertainty to carry forward:
#     # sigma_fin is fitted finite-sample uncertainty; sigma_inf is asymptotic error.
#     # Use sigma_fin (kJ/mol) if you want a finite estimate; both are in same units.
#     return mu_prime, sigma_fin, delta_kJ_cov, cpd

# def multi_dgf_prime_from_sdf(compounds):
#     """Return (mu_prime_kJmol, sigma_kJmol, compound_object) for the SDF file."""
#     mols = [Chem.MolFromMolFile(COMPOUNDS_STRUCTURE_PATH.format(compound=cname, pubchem_id=cdata["pubchem_id"])) for cname, cdata in compounds.items()]
#     if any(mol is None for mol in mols):
#         raise ValueError(f"Could not parse all SDFs")

#     inchis = [Chem.MolToInchi(mol) for mol in mols]
#     if any(not i for i in inchis):
#         raise ValueError(f"RDKit could not generate InChI for some compounds")
    
#     # get equilibrator compound object by InChI
#     cpds = zip(*map(cc.get_compound_by_inchi, inchis))
    
#     mu_list, sigmas_fin, sigmas_inf = zip(*map(cc.standard_dg_formation, cpds))
#     mu_list = numpy.array(mu_list)
#     sigmas_fin = numpy.array(sigmas_fin)
#     sigmas_inf = numpy.array(sigmas_inf)

#     # we now apply the Legendre transform to convert from the standard ΔGf to the standard ΔG'f
#     delta_q_list = numpy.array([
#         cpd.transform(cc.p_h, cc.ionic_strength, cc.temperature, cc.p_mg).m_as("kJ/mol")
#         for cpd in cpds
#     ])
    
#     delta_kJ_list = mu_list + delta_q_list
#     delta_kJ_cov = sigmas_fin @ sigmas_fin.T + 1e6 * sigmas_inf @ sigmas_inf.T

#     # choose an uncertainty to carry forward:
#     # sigma_fin is fitted finite-sample uncertainty; sigma_inf is asymptotic error.
#     # Use sigma_fin (kJ/mol) if you want a finite estimate; both are in same units.
#     return delta_kJ_list, sigmas_fin, delta_kJ_cov, cpds

# def reaction_dg_prime(compound_path_to_stoich):
#     """
#     compound_path_to_stoich: dict mapping SDF file path -> stoichiometric coefficient
#       (reactants negative, products positive), e.g. {'glc.sdf': -1, 'g6p.sdf': 1}
#     Returns (dGprime_kJmol, sigma_kJmol) where sigma is propagated.
#     """
#     total = 0.0
#     var_sum = 0.0
#     for path, stoich in compound_path_to_stoich.items():
#         mu_p, sigma, _ = dgf_prime_from_sdf(path)
#         total += stoich * mu_p
#         var_sum += (stoich * sigma) ** 2   # variance adds if errors independent

#     sigma_rxn = math.sqrt(var_sum)
#     return total, sigma_rxn

#     # Example usage:
#     # rxn = {
#     #     "../compounds/glucose/structure/glucose_5793.sdf": -1,   # reactant
#     #     "../compounds/glucose-6-phosphate/structure/g6p_… .sdf": 1  # product
#     # }
#     # dGprime_val, dGprime_sigma = reaction_dg_prime(rxn)
#     # print("ΔG' (kJ/mol) =", dGprime_val, "+/-", dGprime_sigma)

# # Example:
# if __name__ == "__main__":
    
#     # compute for all compounds in the data.compounds dictionary
    
#     for compound in compounds.keys():
        
#         print(f"Processing {compound}...")
        
#         if compound in ["strictosidine", "serotonin", "tryptamine", 
#                         "tetrahydroharman", "tetrahydroharmol", "tetrahydroharmine",
#                         "harmalol", "harmalan", "harmaline", 
#                         "h+", "nadp+", "nadph"]: 
#             continue
        
#         path = COMPOUNDS_STRUCTURE_PATH.format(compound=compound, pubchem_id=compounds[compound]["pubchem_id"])
#         mu_prime, sigma_fin, delta_kJ_cov, cpd = dgf_prime_from_sdf(path)
        
#         # save to JSON file
#         thermo_path = COMPOUNDS_THERMO_PATH.format(compound=compound, pubchem_id=compounds[compound]["pubchem_id"])
        
#         # Create folder if it doesn't exist
#         os.makedirs(os.path.dirname(thermo_path), exist_ok=True)
        
#         thermo_data = {
#             "mu_prime_kJmol": mu_prime,
#             "sigma_fin_kJmol": sigma_fin.tolist(),
#             "delta_kJ_cov": delta_kJ_cov.tolist(),
#             "cpd_id": cpd.id,
#             "cpd_inchi": cpd.inchi,
#             "cpd_smiles": cpd.smiles
#         }
        
#         with open(thermo_path, "w") as f:
#             json.dump(thermo_data, f, indent=4)
    
#     # print(dgf_prime_from_sdf("../compounds/co2/structure/co2_280.sdf"))
    
#     # from dgpredictor import DGReaction

#     # def get_reaction_dG(reactants, products):
#     #     rxn = DGReaction.from_smiles(reactants, products)
#     #     return rxn.get_dg()

#     # # Example: strictosidine synthase reaction
#     # reactants = ["C10H12N2O (tryptamine SMILES)", "secologanin SMILES"]
#     # products = ["strictosidine SMILES"]

#     # print(get_reaction_dG(reactants, products))
