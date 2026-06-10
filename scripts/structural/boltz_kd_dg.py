#@title (Optional) Calculate Ligand Kd and dG

import json
import math

#@markdown Set the temperature in Celsius.
temperature = 20 #@param {type:"slider", min:0, max:100, step:1}

pred_value = 0.577207088470459
Kd = math.pow(10, pred_value)/1000000
dG = 8.3144626*(273.15+temperature)*math.log(Kd)/1000

print(f'Predicted dissociation eq constant  : {Kd:.3e} mol/L  (= {Kd*1000000:.1f} uM)')
print(f'Predicted binding Gibbs free energy : {dG:.3f} kJ/mol')
print(f'                                      {dG/4.184:.3f} kcal/mol')

