clear classes
clear all
close all
clc

rehash path

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

return

%% Dynamic simulation
% t  = 0:1E-3:10*1E-3;
% dt = mean(diff(t));
% % F  = 10000*rand(size(myStrModel.Fs,1),length(t)); % Random force (SIZE --> Fs x Ntimesteps)
% F = zeros(size(myStrModel.Fs,1),length(t));
% F(3,:)     = 1E6*sin(10*pi*t + 1); % Sinusoidal tip force
% F(end-3,:) = 1E6*sin(10*pi*t + 1); % Sinusoidal tip force
% myStrModel_tsim = myStrModel.initializeDynSimulation(myStrModel.disp,t); % Use myAE_sol.strModel.disp to initialize at static trim point.
% myStrModel_tsim = myStrModel_tsim.solveDynamics(F,t,dt);

%% dyn aero sim (start from rest)
% dt=0.0011;
% t=0:dt:0.2;
% initstate= class_aero_state(100, 0,   0,   0.0,   1.225);
% myAeroModel=myAeroModel.initializeDynSimulation(initstate,dt);
% 
% for i=1:length(t)
%     myAeroModel=myAeroModel.solveTimeStep(state,i);  
% end

%% aeroelastic dyn sim
myAE = aeroElasticModeldAED(myAeroModel,myStrModel);
state = class_aero_state(100, 10,   0,   0.0,   1.225);
myAE_tsim  = myAE.solveDynamic(state,zeros(366,1),0:1E-3:1);
