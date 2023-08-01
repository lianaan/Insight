# Adapted from Simple demonstration of PyDDM.
import scipy.io
from pyddm import Sample
import pandas
import numpy as np
# Create a simple model with constant drift, noise, and bounds + add initial conditions (ICPointSourceCenter).
from pyddm import Model
from pyddm.models import DriftConstant, NoiseConstant, BoundConstant, OverlayNonDecision, OverlayChain, ICPointSourceCenter, OverlayPoissonMixture, LossRobustLikelihood
from pyddm.functions import fit_adjust_model, display_model
import pyddm as ddm
from pyddm import InitialCondition


from pyddm import Fittable, Fitted
from pyddm.models import LossRobustBIC
from pyddm.functions import fit_adjust_model

from pyddm import FitResult
from pyddm.functions import solve_partial_conditions

class DriftCoherence(ddm.models.Drift):
    name = "Drift depends linearly on coherence"
    required_parameters = ["driftcoh"] # <-- Parameters we want to include in the model
    required_conditions = ["coh"] # <-- Task parameters ("conditions"). Should be the same name as in the sample.
    
    # We must always define the get_drift function, which is used to compute the instantaneous value of drift.
    def get_drift(self, conditions, **kwargs):
        return self.driftcoh * conditions['coh']

class ICPointSideBiasRatio(InitialCondition):
     name = "A side-biased starting point expressed as a proportion of the distance between the bounds."
     required_parameters = ["x0"]
     required_conditions = ["left_is_correct"]
     def get_IC(self, x, dx, conditions):
         x0 = self.x0/2 + .5 #rescale to between 0 and 1
         # Bias > .5 for left side correct, bias < .5 for right side correct.
         # On original scale, positive bias for left, negative for right
         if not conditions['left_is_correct']:
             x0 = 1-x0
         shift_i = int((len(x)-1)*x0)
         assert shift_i >= 0 and shift_i < len(x), "Invalid initial conditions"
         pdf = np.zeros(len(x))
         pdf[shift_i] = 1. # Initial condition at x=x0*2*B.
         return pdf

model_base = Model(name='drift varies with coherence, also IC',
                 drift=DriftCoherence(driftcoh=Fittable(minval=0, maxval=20)),
                 noise=NoiseConstant(noise=1),
                 IC=ICPointSideBiasRatio(x0=Fittable(minval=-1, maxval=1)),
                 bound=BoundConstant(B=Fittable(minval=.1, maxval=3)),
                 overlay=OverlayChain(overlays=[OverlayNonDecision(nondectime=Fittable(minval=0, maxval=2)),
                                                OverlayPoissonMixture(pmixturecoef=.02,
                                                                      rate=1)]),
              dx=.01, dt=.0005, T_dur=2.999)

cond_all={"coh": [0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050,0.055, 0.060,0.065,0.070, \
     0.075,0.080, 0.085, 0.090, 0.095, 0.100, 0.105, 0.110, 0.115, 0.120, 0.125, 0.130, 0.135, 0.140, \
     0.145, 0.150, 0.155, 0.160, 0.165, 0.170, 0.175,0.180, 0.185, 0.190, 0.195,0.200, 0.205, 0.210, \
     0.215, 0.220, 0.225, 0.230, 0.235, 0.240, 0.245, 0.250, 0.255, 0.260, 0.265, 0.270, 0.275, 0.280, \
     0.285, 0.290, 0.295, 0.300],"left_is_correct": [0,1], "x0":[0,1]}
    

output_list_choices = []
output_list_RT = []
    

for s in range(1,23):
    print(s)
    for c in range(1,5): 
        print(c)

        with open("mae_s_" + str(s)  + "_c_" + str(c)  +".csv", "r") as f:
            df_rt = pandas.read_csv(f)
        
            df_rt = df_rt[df_rt["participant"] == 1] # Only monkey 1
              
            # Remove short and long RTs, as in 10.1523/JNEUROSCI.4684-04.2005.
            # This is not strictly necessary, but is performed here for
            # compatibility with this study.
            df_rt = df_rt[df_rt["rt"] > .1] # Remove trials less than 100ms
            df_rt = df_rt[df_rt["rt"] < 3] # Remove trials greater than 5000ms
              
            
            roitman_sample = Sample.from_pandas_dataframe(df_rt, rt_column_name="rt", correct_column_name="correct")
         
            # Load the model
            with open("mae_s_" + str(s) + "_c_" + str(c) + ".csv_model_IC_output.txt", "r") as g:
                model_loaded = eval(g.read())
                g.close()
        
            
            RT_pred = [];
            prob_corr_pred = [];
            for j in range(0, 121):
                if j in df_rt.index:
                    coh_val = df_rt.coh[j];
                    left_is_corr_val = df_rt.left_is_correct[j];
                    tr_corr = df_rt.correct[j];  
                    trial_pred = solve_partial_conditions(model_loaded,conditions={"coh": [coh_val], "left_is_correct": [left_is_corr_val]});
                    prob_corr_pred.append(trial_pred.prob_correct()) 
                    if tr_corr == 1: 
                      RT_pred.append(trial_pred.mean_decision_time())
                    else: 
                        tp_pdf_err = trial_pred.pdf_err()
                        tp_t_domain = model_loaded.t_domain()
                        val_time = sum(tp_pdf_err *tp_t_domain)/ sum(tp_pdf_err)
                        RT_pred.append(val_time) 
                else:
                    prob_corr_pred.append(np.nan)
                    RT_pred.append(np.nan)
                    
            output_list_choices.append(prob_corr_pred)
            output_list_RT.append(RT_pred)
            
            f.close()
           
scipy.io.savemat("output_list_RT_allEE.mat", {"output_list_RT": output_list_RT})  
scipy.io.savemat("output_list_choices_allEE.mat", {"output_list_choices": output_list_choices})
     



