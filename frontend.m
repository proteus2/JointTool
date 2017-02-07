clear classes
clear all
close all
clc

%% add required modules
addpath(genpath('ext'))
addpath(genpath('linStaticAeroModule'))
addpath(genpath('linStaticStrModule'))
addpath(genpath('aeroStrModule'))

%% think about a way of putting states into an/another input file
%                       U   aoa beta    Ma      rho
 state = class_aero_state(100, 10*pi/180, 0,      0.0,    1.225);
 
%% Setup Aero Model
myAeroModel=aeroModel('crmWing4dAEDalus.xml',state);
myAeroModel=myAeroModel.calculateForces();
myAeroModel.plotCp();

%% Setup Strutural Model
% Class init.
myStrModel = structuralModel('CRM');

%% Setup Aerostructural Model
 myAeroStructModel=aeroStrModel(myAeroModel,myStrModel);