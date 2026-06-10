#@title (Optional) Calculate Ligand Kd and dG

import json
import math

#@markdown Set the temperature in Celsius.
temperature = 20 #@param {type:"slider", min:0, max:100, step:1}

pred_value = 0.577207088470459
Kd = math.pow(10, pred_value)/1000000
dG = 8.3144626*(273.15+temperature)*math.log(Kd)/1000
IC50 = 0

print(f'Predicted dissociation eq constant  : {Kd:.3e} mol/L  (= {Kd*1000000:.1f} uM)')
print(f'Predicted binding Gibbs free energy : {dG:.3f} kJ/mol')
print(f'                                      {dG/4.184:.3f} kcal/mol')
print(f'Predicted IC50                      : {IC50:.3f} kJ/mol')


# There are two main predictions in the affinity output: affinity_pred_value and affinity_probability_binary. They are trained on largely different datasets, with different supervisions, and should be used in different contexts.

# The affinity_probability_binary field should be used to detect binders from decoys, for example in a hit-discovery stage. It's value ranges from 0 to 1 and represents the predicted probability that the ligand is a binder.

# The affinity_pred_value aims to measure the specific affinity of different binders and how this changes with small modifications of the molecule (note that this implies that it should only be used when comparing different active molecules, not inactives). This should be used in ligand optimization stages such as hit-to-lead and lead-optimization. It reports a binding affinity value as log10(IC50), derived from an IC50 measured in μM. Lower values indicate stronger predicted binding, for instance:

# IC50 of 10e-9 M ⟶our model outputs 3 (strong binder)
# IC50 of 10e-6 M ⟶ our model outputs 0 (moderate binder)
# IC50 of 10e-4 M ⟶ our model outputs 2 (weak binder / decoy)

# You can convert the model's output to pIC50 in kcal/mol by using y --> (6 - y) * 1.364 where y is the model's prediction.
