clear classes
clear all
close all
clc

%% add required modules
addpath(genpath('ext'))
addpath(genpath('linStaticAeroModule'))
addpath(genpath('linStaticStrModule'))
addpath(genpath('aeroElasticModule'))

%% think about a way of putting states into an/another input file
%                       U   aoa beta    Ma      rho
 state = class_aero_state(100,0.17453*180/pi, 0,      0.0,    1.225);
%% Setup Aero Model

myAeroModel=aeroModel('rectangularWing4dAEDalus.xml',state);
myAeroModel=myAeroModel.calculateForces();
%% Setup Strutural Model
% Class init.
myStrModel = structuralModel('recBeam');

%% Setup Aerostructural Model
myAE=aeroElasticModeldAED(myAeroModel,myStrModel);
 %%
% Solve aerostruct Model
myAE=myAE.solve();
