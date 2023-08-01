# modified from PyDDM samples
from pyddm import Sample
import pandas
import os
import pyddm as ddm
import numpy as np
from pyddm import InitialCondition
from pyddm import Model, Fittable
from pyddm.functions import fit_adjust_model, display_model
from pyddm.models import NoiseConstant, BoundConstant, OverlayChain, OverlayNonDecision, OverlayPoissonMixture, LossRobustLikelihood
from pyddm import Drift  

class DriftCoherenceRewBias(Drift):
    name = "Drift depends linearly on coherence, with a reward bias"
    required_parameters = ["driftcoh", "rewbias"] # <-- Parameters we want to include in the model
    required_conditions = ["coh", "highreward"] # <-- Task parameters ("conditions"). Should be the same name as in the sample.

# We must always define the get_drift function, which is used to compute the instantaneous value of drift.
    def get_drift(self, conditions, **kwargs):
        rew_bias = self.rewbias * (1 if conditions['highreward'] == 1 else -1)
        return self.driftcoh * conditions['coh'] + rew_bias



path = '/Users/horgalab/Documents/Python'
dir_list = os.listdir(path)
csv_list = [];
for x in dir_list:
    if x.endswith(".csv"):
        csv_list.append(x)
print(csv_list)        
for file in csv_list:
    #with open("mae_s_6_c_4.csv", "r") as f:
    with open(file, "r") as f:
        df_rt = pandas.read_csv(f)
    
    df_rt = df_rt[df_rt["participant"] == 1] # Only monkey 1
      
    # Remove short and long RTs, as in 10.1523/JNEUROSCI.4684-04.2005.
    # This is not strictly necessary, but is performed here for
    # compatibility with this study.
    df_rt = df_rt[df_rt["rt"] > .1] # Remove trials less than 100ms
    df_rt = df_rt[df_rt["rt"] < 3] # Remove trials greater than 5000ms
      
    # Create a sample object from our data.  This is the standard input
    # format for fitting procedures.  Since RT and correct/error are
    # both mandatory columns, their names are specified by command line
    # arguments.
    roitman_sample = Sample.from_pandas_dataframe(df_rt, rt_column_name="rt", correct_column_name="correct")
    
    
    model_DB = Model(name='drift varies with coherence',
                     drift=DriftCoherenceRewBias(driftcoh=Fittable(minval=0, maxval=20), rewbias=Fittable(minval=-1, maxval=1)),
                     noise=NoiseConstant(noise=1),
                     bound=BoundConstant(B=Fittable(minval=.1, maxval=3)),
                     overlay=OverlayChain(overlays=[OverlayNonDecision(nondectime=Fittable(minval=0, maxval=2)),
                                                    OverlayPoissonMixture(pmixturecoef=.02,
                                                                          rate=1)]),
                  dx=.01, dt=.0005, T_dur=2.999)
   
    
    fit_model_DB = fit_adjust_model(sample=roitman_sample, model=model_DB, lossfunction=LossRobustLikelihood, verbose=False)
    
    display_model(fit_model_DB)
    
    with open(file+"_model_DB_output.txt", "w") as f2:
        f2.write(repr(fit_model_DB))
    f2.close()
    f.close()
