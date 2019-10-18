close all
clear all
%---------------simulation parameters---------
nb_simul=1; %number of realizations to simulate
simul_type=1; %model to use: 1=parametric; 2=non-parametric

%-------------------load data------------------
%In this example we use a 17-year long rain type dataset for calibration and simulation.
%The model is calibrated using the 10 first years of the dataset.
%The simulation is performed for the years 11 to 15.
%The rain type data corresponding to years 11 to 15 are kept out of the simulation but can be used to assess the quality of the simulation. These data are stored in V_raintypes_valid
cd('Raintype_data\')
M_Cov=csvread('M_Cov_ERA5.csv');
M_Cov_calib=M_Cov(1:52560*10,:);
Cov_simul=M_Cov(52560*10+1:52560*15,:);
M=load('M_datevec.mat');
M_datevec_calib=M.M_datevec(1:52560*10,:);
M=load('V_raintypes_6clusters.mat');
V_raintypes_calib=M.V_raintypes(1:52560*10,:);
V_raintypes_valid=M.V_raintypes(52560*10+1:52560*15,:);
cd('..')
%----------------run simulation--------------------
%!!!! g2s (i.e. the interface of the MPS simulation software) must be installed before running simul_type=2. See https://github.com/GAIA-UNIL/G2S for details about the installation of g2s

[M_raintypes_simul]=Stochastic_raintype_simulation(M_Cov_calib, V_raintypes_calib, Cov_simul, nb_simul, simul_type);
%M_raintypes_simul contains the results of the simulation. lines are different realizations, columns are rain type time series 

