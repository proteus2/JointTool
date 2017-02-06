function [constant] = AddFuel(constant,fuel_level,GRAV,grav_cst,wingID)

% =====                                                              ==== %
%                          Fuel Input [14/04/2015]                        %
%                                                                         %
%  Define the location and mass of fuel (Kg)
% =====                                                              ==== %

global GRAPHFLAG xlsFileName

fext   = constant.fext;
lumped = constant.lumped;

% Load Fuel Data
% fueldata.gri --> x   y     z
% fueldata.mas --> m   I11   I22   I33   I21   I23   I13

fueldata       = xlsread(xlsFileName,'FuelData');

for i=1:length(constant.lumped.type)
    if strcmp(constant.lumped.type{i},'Ribs')
        RibLumNum = i;
        break
    end 
end

fuelind = ~isnan(fueldata(:,1));
fuel_ID = fueldata(fuelind,1);
if strcmp(wingID,'right')
    fuel_xyz(:,2) = fueldata(fuelind,3);
elseif strcmp(wingID,'left')
    fuel_xyz(:,2) = -fueldata(fuelind,3); 
end
fuel_xyz(:,3) = fueldata(fuelind,4);
fuel_xyz(:,1) = fueldata(fuelind,2);
fuel_range = fueldata(1,5:end);
fuel_mass = fueldata(fuelind,5:end);

if ~exist('RibLumNum') && sum(sum(isnan(fuel_xyz(:,[1,3]))))>0
    error('Please specify ribs in order to position the fuel masses')
end

% ----------------- Take into account as lumped mass ------------------- %

EntryNum = length(constant.lumped.type);
EntryNum_fext = length(constant.fext.type);

Ntanks = max(fuel_ID);

lumped.mass{EntryNum+Ntanks}=[];
lumped.location{EntryNum+Ntanks}=[];
lumped.IG{EntryNum+Ntanks}=[];
if GRAV
    fext.magnitude{EntryNum_fext+Ntanks} = [];
    fext.location{EntryNum_fext+Ntanks}  = [];
    fext.follower{EntryNum_fext+Ntanks}  = [];
    fext.alphaflag{EntryNum_fext+Ntanks}  = [];
end
for idx = 1:length(fuel_ID)
    
    % Store Mass Input (Lumped Formulation)
    lumped.type{EntryNum+fuel_ID(idx)}	      = ['Fuel_Tank_',num2str(fuel_ID(idx))];
    lumped.mass{EntryNum+fuel_ID(idx)}(end+1,1)   = interp1(fuel_range,fuel_mass(idx,:),fuel_level(fuel_ID(idx))); 
    massloc = fuel_xyz(idx,:);
    
    % Shift massloc with the xyzshift
    massloc =  massloc-constant.inp.xyzshift;
    
    % Check whether massloc is completely defined, otherwise interpolate
    % based on rib centroids
    if isnan(massloc(2))
       error('Please specify spanwise fuel location') 
    end
    if isnan(massloc(1))
       massloc(1) = interp1(constant.lumped.location{RibLumNum}(:,2),constant.lumped.location{RibLumNum}(:,1),massloc(2)); 
    end
    if isnan(massloc(3))
       massloc(3) = interp1(constant.lumped.location{RibLumNum}(:,2),constant.lumped.location{RibLumNum}(:,3),massloc(2)); 
    end    
        
    lumped.location{EntryNum+fuel_ID(idx)}(end+1,:)    = massloc;
    
	% ATTENTION: disable inertia. Unreliable input!
    IG = zeros(3,3);
    lumped.IG{EntryNum+fuel_ID(idx)}(:,end+(1:3))=IG;
    
    if GRAV
        fext.type{EntryNum_fext+fuel_ID(idx)} = ['Fuel Tank # ',num2str(fuel_ID(idx))];
        fext.magnitude{EntryNum_fext+fuel_ID(idx)}(end+1,:) = [0 0 -constant.general.nz*lumped.mass{EntryNum+fuel_ID(idx)}(end,1)*grav_cst 0 0 0];
        fext.location{EntryNum_fext+fuel_ID(idx)}(end+1,:)  = massloc;
        fext.follower{EntryNum_fext+fuel_ID(idx)}(1,end+1)  = 0;
        fext.alphaflag{EntryNum_fext+fuel_ID(idx)}(1,end+1)  = 1;
    end
end

lumped = mext_inp (constant.str,lumped);

% Output
constant.lumped = lumped;
constant.fext = fext;

FuelMass = eps;
for i=EntryNum+1:length(constant.lumped.type)
    FuelMass = FuelMass+sum(constant.lumped.mass{i});
end

% ---
if 1 && GRAPHFLAG == 1    % plot lumped mass all fuell as cube
    figure(1)
    colours = {'yellow','green','white','magenta','blue'};
    for i = EntryNum+1:length(constant.lumped.type)
        for j=1:length(constant.lumped.mass{i})
            PlotCube(constant.lumped.location{i}(j,1:3),constant.lumped.mass{i}(j)/FuelMass*constant.inp.xyz(end-1)/2.5,colours{mod(i-EntryNum-1,length(colours))+1})
        end
    end   
end
% ---

function PlotCube ( origin, size ,c)
x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*size+origin(1);
y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*size+origin(2);
z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*size+origin(3);
for i=1:6
    h=patch(x(:,i),y(:,i),z(:,i),c);
    set(h,'edgecolor','k')
end