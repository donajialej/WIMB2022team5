# WIMB2022team5
The repository for team 5 from the Collaboraitve Workshop on Women Health 2022, at Minnesota, USA. The team leaders: Ruiyan Luo (Department of Population Health Sciences, School of Public Health) and Alexandra Smirnov (Department of Mathematics \& Statistics) from Georgia State University, Atlanta, USA. The rest of the team: , ... , ....,  Alejandra D. Herrera-Reyes (School of Mathematical Sciences, University of Nottingham, Nottingham, UK), and Diana White (Department of Mathematics & Statistics, Clarkson University, Potsdam, NY, USA).

Our work is titled: __Parameter Estimation for COVID-19 SVIRD   Model Using Predictor-Corrector Algorithm__. 

To run this repository, you will need the [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html) and the [Optimization Toolbox](https://www.mathworks.com/products/optimization.html) from Matlab.

In `/fits`, you can find our estimates and bootstraps for all 4 data sets. They were geenerated with the following codes:

- Synthetic scenario 1: data set is `SynthDI2R.txt` and the  fitting code is `Param_Est_SVIRD_PCM_DI2_new.m`.
- Synthetic scenario 2: data set is `SynthDI1R.txt` and the fitting code is `Param_Est_SVIRD_PCM_DI1_new.m`.
- Georgia: data set is `GA_inc_cases_deaths_7_9_21_11_25_21_CDC_7day_average.txt` and the fitting code is `Param_Est_SVIRD_PCM_GA_normalized3.m`.
- California: data set is `California_7_9_21_11_25_21.txt` and the fitting code is `Param_Est_SVIRD_PCM_CAL_new.m`.


After running the previous codes, the figures can be generated using  `Figures_SVIRDpaper.m`, except for Figure 8 and Figure 9 in the paper. In `/figures`, you can find a .eps and a Matlab figure version of all the figures in the paper. 

To generate figure 8 and 9 run:
- `Param_Est_SVIRD_PCM_GA_mult_rep_rates` for GA
- `Param_Est_SVIRD_PCM_CAL_mult_rep_rates` for CA.
