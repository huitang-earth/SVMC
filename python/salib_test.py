from SALib.sample import sobol  as sobol_sample 
from SALib.analyze import sobol as sobol_analyse
from SALib.test_functions import Ishigami
import xarray as xr
import numpy as np
import datetime as dt

operation=2
input_file="/projappl/project_2006422/svm_data/output/veg_sensitivty/SVM_Qvidja_hr_vegsen.all.C5.nc"


problem = {
    'names': ['conduct', 'psi', 'alpha', 'gamma'],
    'num_vars': 4,
    'bounds': [[0.5e-16, 10e-16], [-3.0, -1], [0.03, 0.12], [0.1, 5]],
    'dists': ['unif', 'unif', 'unif', 'unif']
}

if operation==1:
  param_values = sobol_sample.sample(problem, 128)
  np.savetxt("param_values.txt", param_values)


if operation==2:
  ## Run model (example)
  #Y = Ishigami.evaluate(param_values)
  #print(Y)

  ## Perform analysis
  svm=xr.open_dataset(input_file)
  
  ## Variables to be used for sensitivity analaysis:
  # Selected daily/monthly GPP or LH or Soil moisture
  # Overall performance (yearly, or the whole periods): correlation/biases with observed GPP, LH, soil moisture?

  sc_start = dt.datetime(2019, 5, 1, 0, 0) 
  sc_end   = dt.datetime(2019, 9, 30, 23, 59)

  Y=svm["GPP"].sel(time=slice(sc_start, sc_end)).resample(time='1M').mean()
  print(Y)

  for i in range(0,5):
      Si = sobol_analyse.analyze(problem, Y[:,i,0,0].data, print_to_console=True)
#      print(Si)

  # Print the first-order sensitivity indices

#Y = np.zeros([param_values.shape[0]])
#for i, X in enumerate(param_values):
#    Y[i] = Ishigami.evaluate(X)
#np.savetxt("output.txt", Y)
#Y = np.loadtxt("outputs.txt", float)
#Si = sobol.analyze(problem, Y)


