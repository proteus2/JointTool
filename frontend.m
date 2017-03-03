clear classes
clear all
close all
clc

%% add required modules
addpath(genpath('ext'))
addpath(genpath('staticAeroModule'))
addpath(genpath('strModule'))
addpath(genpath('aeroElasticModule'))

%% think about a way of putting states into an/another input file
%                          U  aoa beta    Ma      rho
state = class_aero_state(100, 10,   0,   0.0,   1.225);
 
%% Setup Aero Model
% crmWing4dAEDalus
% rectangularWing4dAEDalus
myAeroModel=aeroModel('rectangularWing4dAEDalus.xml',state);

%% Setup Strutural Model
myStrModel = structuralModel('RecBeam');

%% Setup Aerostructural Model
% myAE = aeroElasticModeldAED(myAeroModel,myStrModel);
% myAE = aeroElasticModelMatrix(myAeroModel,myStrModel);

%% Solve aerostruct Model
% myAE_sol = myAE.solve();
% myAE_sol = myAE.solve();

%% Linear solution (For comparison)
% myAE_linsol = myAE;
% myAE_linsol.strModel.settings.lin = 1;
% myAE_linsol = myAE_linsol.solve;
% myAE_linsol.plotResults

%% Dynamic simulation
t  = 0:100;
dt = (t(end) - t(1))/length(t);
F  = rand(size(myStrModel.Fs,1),length(t)); % Random tip force (SIZE --> Fs x Ntimesteps)
myStrModel = myStrModel.initializeDynSimulation(myStrModel.disp); % Use myAE_sol.strModel.disp to initialize at static trim point.
myStrModel_tsim = myStrModel.solveTimeStep(F,t,dt);

