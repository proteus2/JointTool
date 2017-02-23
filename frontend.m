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
myAeroModel=aeroModel('rectangularWing4dAEDalus.xml',state);

%% Setup Strutural Model
myStrModel = structuralModel('RecBeam');

%% Setup Aerostructural Model
myAE = aeroElasticModeldAED(myAeroModel,myStrModel);
% myAE = aeroElasticModelMatrix(myAeroModel,myStrModel);

%% Solve aerostruct Model
myAE_nlinsol = myAE.solve();
% myAE_nlinsol = myAE.solve();

%% Linear solution (For comparison)
% myAE_linsol = myAE;
% myAE_linsol.strModel.settings.lin = 1;
% myAE_linsol = myAE_linsol.solve;
% myAE_linsol.plotResults
