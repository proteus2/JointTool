clear classes
clear all
close all
clc

%% add required modules
addpath(genpath('ext'))
addpath(genpath('linStaticAeroModule'))
addpath(genpath('linStaticStrModule'))

%% think about a way of putting states into an/another input file
%                       U   aoa beta    Ma      rho

%% Setup Aero Model


%% Setup Strutural Model
% Class init.
myStrModel = structuralModel('CRM');