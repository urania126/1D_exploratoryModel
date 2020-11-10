# 1D_exploratoryModel
(10/11/2020)  
We built this repository to make our source code available for everybody. 
You will find 4 core files: 
- `binormal.m`: a function to perform normal approximation to binomial
distribution when number of elements is over 1000. 
- `variables.m`: list of all parameters used in the model, including
models' scaffold, cellular skills, replicates tracking and data storing structures. 
- `dataProcessing.m`, which we use to calculate mean, 
standard deviation and mean standard error for all variables of interest. We also use this file to 
manipulate files and store results in local directories. 
- `plots_replicates.m`: code for generating plots, including population growth plot, 
a colormap for phenotypes frequencies, 
average rho evolution over whole simulation time, 
barplots for history of phenotypic frequencies, 
evolution of average rho in physical space
and tumor volume determination (when feasible).

We planned three different modelling approaches and 
each of them is store in an independent file: 
- Mass-action model considering unbounded growth: `model1Drep_nomigration_nosaturation.m`. 
- Mass-action model considering logistic growth: `model1Drep_nomigration_saturation.m`. 
- 1D Spatial model considering logistic growth: `model1Drep_migration_saturation.m`. 

Julia code files will be available soon. 
