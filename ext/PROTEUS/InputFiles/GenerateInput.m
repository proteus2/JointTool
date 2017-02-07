% =====                                                              ==== %
%                        CRM Input File [22/10/2015]                      %
%                                                                         %
%   Generates .txt files defining aerodynamic profiles, wing box nodes
%   location, elements' connectivity and laminate props.
%   Plots wing profiles, indicates spanwise location, shows wing box
% =====                                                              ==== %

function [constant] = GenerateInput(model)

global GRAPHFLAG FLAGSPAR TOTNSPARS xlsFileName
constant = [];

AircraftData              = xlsread(xlsFileName,'AircraftData');
[wing_data]               = xlsread(xlsFileName,'WingData','A1:T99');
wing_1g_data              = xlsread(xlsFileName,'Wing1g','A1:T99');
ModellingInput            = xlsread(xlsFileName,'ModellingInput');
fueldata                  = xlsread(xlsFileName,'FuelData');
[Num,~]                   = xlsread(xlsFileName,'NonStructuralMasses');

BUCKL   = ModellingInput(33,1);
TWIST1G = ModellingInput(34,1);

%=========================================================================%
% Define aircraft data
%=========================================================================%
% Half aircraft weight without wing
ac_weight = 0.5*sum(AircraftData(~isnan(AircraftData(:,1)),1))*9.81;% half aircraft weigth withouth wing (remove MLG, Engine/Pylon)
MTOW = 2*sum(Num(:,1))+2*sum(fueldata(3:end,end))+2*ac_weight/9.81;% Excluding wing mass
MZFW = 2*sum(Num(:,1))+2*ac_weight/9.81;% Excluding wing mass
MLW = 2*sum(Num(:,1))+AircraftData(4,4)*2*sum(fueldata(3:end,end))+2*ac_weight/9.81;% Excluding wing mass
LW = 2*sum(Num(:,1))+AircraftData(4,4)*2*sum(fueldata(3:end,end))+2*ac_weight/9.81;% Excluding wing mass
TSFC = AircraftData(1,4); % Thrust specific fuel consumption
ZMO = AircraftData(2,4); % Maximum operating altitude
IWW = AircraftData(3,4); % Initial guess for the wing weight to determine the gust amplitude

%=========================================================================%
% Wing geometry
%=========================================================================%
% Maximum AOA allowed
if ModellingInput(30,1)==1
    alpha_max   = deg2rad(ModellingInput(30,2));     % maximum AOA allowed
else
    alpha_max   = deg2rad(15);     % maximum AOA allowed
end

% Process wing data
xyz_LE_right_0   = wing_data(:,[3,2,4]);

% Apply dihedral (put as input)
Lambda = 6; % degrees
Rdih   = [1 0 0; 0 cosd(Lambda) -sind(Lambda); 0 sind(Lambda) cosd(Lambda)];

xyz_LE_right = (Rdih*(xyz_LE_right_0 - repmat(xyz_LE_right_0(1,:),size(xyz_LE_right_0,1),1))')' + repmat(xyz_LE_right_0(1,:),size(xyz_LE_right_0,1),1);
clear xyz_LE_right_0

xyz_LE_left      = xyz_LE_right.*repmat([1 -1 1],size(xyz_LE_right,1),1);
xyz_LE_left(1,:) = []; 
xyz_LE           = [flipud(xyz_LE_left); xyz_LE_right];

theta_right = deg2rad(wing_data(:,8));
theta_left  = theta_right(2:end);
theta       = [flipud(theta_left); theta_right];

c_right = wing_data(:,5);
c_left  = c_right(2:end);
c       = [flipud(c_left); c_right];

xyz_TE = xyz_LE+[cos(theta).*c,zeros(size(xyz_LE,1),1),-sin(theta).*c];

xyzqc  = xyz_LE+[0.25*cos(theta).*c,zeros(size(xyz_LE,1),1),-0.25*sin(theta).*c];

% Reference axis location w.r.t. LE in percentage of chord
xref_right = wing_data(:,6);
xref_left  = xref_right(2:end);
xref       = [flipud(xref_left); xref_right];

% Structural cross-section origin wrt LE in percentage chord (0 means on xref itself)
xsref_right = wing_data(:,7); 
xsref_left  = xsref_right(2:end);
xsref       = [flipud(xsref_left); xsref_right];

% Beam axis calculation
xyz = xref(:,ones(3,1)).*xyz_TE+(1-xref(:,ones(3,1))).*xyz_LE;       

% Boolean vector defining fixed nodes that are preserved in the structure
fixnodes_right = wing_data(:,11);
fixnodes_left  = fixnodes_right(2:end);
fixnodes       = [flipud(fixnodes_left); fixnodes_right];

% xyzshift = [xyz(1,1);0;xyz(1,3)];
xyzshift = xyz(xyz(:,2)==0,:);
    
xyz_LE(:,1) = xyz_LE(:,1) - xyzshift(1);
xyz_LE(:,3) = xyz_LE(:,3) - xyzshift(3);
% xyz_TE(:,1) = xyz_TE(:,1) - xyzshift(1);
% xyz_TE(:,3) = xyz_TE(:,3) - xyzshift(3);
xyzqc(:,1)  =  xyzqc(:,1) - xyzshift(1);                              
xyzqc(:,3)  =  xyzqc(:,3) - xyzshift(3);                                
xyz(:,1)    =    xyz(:,1) - xyzshift(1);                                 
xyz(:,3)    =    xyz(:,3) - xyzshift(3); 

% Manage flags
if size(wing_data,2)>11
    FLAGSPAR = 1;
    TOTNSPARS = size(wing_data,2)-11;
else
    FLAGSPAR = 0;
end

% Process 1g data
if TWIST1G == 1
    theta1g = wing_1g_data(:,3);
    y1g     = wing_1g_data(:,2);
end

% Plot wing planform (IF ACTIVE)
if 1 && GRAPHFLAG == 1  
    figure (1)
    hold all
    xlabel('x')
    ylabel('y')
    zlabel('z')
    plot3(xyz(:,1),xyz(:,2),xyz(:,3),'green --')
    view([45,30])
    axis equal
end

%=========================================================================%
% Load airfoil data
%=========================================================================%
constant = LoadAerofoilProfiles(constant,xyz);

% Plot the input wing surface
constant = PlotAeroSurf(constant,wing_data,xyz,theta,c,xref,xsref,FLAGSPAR);

%=========================================================================%
% Laminate distribution
%=========================================================================%
NlamSpan_right = ModellingInput(1,~isnan(ModellingInput(1,:)));
NlamSpan_left  = NlamSpan_right;
Nlam.Span      = [fliplr(NlamSpan_left) NlamSpan_right];

NlamChord_right = ModellingInput(2,~isnan(ModellingInput(2,:)));
NlamChord_left  = NlamChord_right;
Nlam.Chord      = [fliplr(NlamChord_left) NlamChord_right];

%=========================================================================%
% Check input
%=========================================================================%
% Nlam = CheckInput(wing_data,Nlam);

%=========================================================================%
% Include ribs
%=========================================================================%
[lumped] = AddRibs(wing_data,constant,xyz,theta,xref,xsref,c);

%=========================================================================%
% Wing discretisation
%=========================================================================%
Ns = 2*ModellingInput(3,1)-1; % # of structural  elements
Na = 2*ModellingInput(4,1)-1; % # of aerodynamic elements

[str,aero,general,Nlam] = prepr_stat(xyz,xyzqc,theta,c,xref,Ns,Na,fixnodes,Nlam,lumped,model);      % Format code inputs

if TWIST1G == 1
	str.theta1g = interp1(y1g,theta1g,str.xyz(2:3:end),'pchip');
end

aero.nbs        = Na;                    % Number of spanwise panels
aero.nbc        = ModellingInput(5,1);   % Number of chordwise panels
aero.lw         = ModellingInput(6,1);   % Length wake in number of chords
general.symasym = 1;                     % Symmetry
aero.Udtc       = 1/ModellingInput(7,1); % 1 divided by the number of timesteps per chord travelled, the wake panel size is determined by the minimum chord

%=========================================================================%
% Process lumped data ribs
%=========================================================================%
lumped = mext_inp(str,lumped);

cross.numel = ModellingInput(8,1); % Specify minimum number of cross-sectional elements
if cross.numel ~= 0
%=========================================================================%
% Laminate locations and numbering
%=========================================================================%
    [constant,Nlam] = GenerateLaminateDistribution(wing_data,constant,xyz,theta,c,Nlam,FLAGSPAR);
    
%=========================================================================%
% Add stringers
%=========================================================================%
    [stringer,FLAGSTRING] = AddStringer();
    
%=========================================================================%
% Generate wingbox and corresponding laminates
%=========================================================================%
    [cross,constant,panels] = GenerateWingbox(wing_data,constant,str,cross,xyz,xyz_LE,theta,xref,xsref,c,Nlam,stringer,FLAGSPAR,FLAGSTRING);
    
%=========================================================================%
% Process the buckling panels from panels (CHECK for double wing)
%=========================================================================%
    if BUCKL == 1
        [buckl] = prepr_buckl(constant,str,lumped,panels);
    end
    
%=========================================================================%
% Define the material properties
%=========================================================================%
    cd(['input_' model])
    constant = DefineMaterialProp(constant);                                   
    cd ..

%=========================================================================%
% Prepare material data for strain computation
%=========================================================================%
    constant = prepr_strain(constant);
    
else
    
%=========================================================================%
% Define the material properties
%=========================================================================%
    cd(['input_' model])
    [cross] = DefineBeamProp(cross);
    constant.mat.rho = cross.mat.rho;
    cd ..
    
end

%=========================================================================%
% Define non-structural masses
%=========================================================================%
[lumped] = AddNonStructuralMasses(xyzshift,constant,lumped,'right');
[lumped] = AddNonStructuralMasses(xyzshift,constant,lumped,'left');
[lumped] = shiftNSM(lumped,xyzshift,str);

fext = [];
[fext] = AddExternalForces(fext,xyzshift,constant,'right');
[fext] = AddExternalForces(fext,xyzshift,constant,'left');
[fext] = shiftFEXT(fext,xyzshift,str);

%=========================================================================%
% Define initial laminate distribution
%=========================================================================%
if cross.numel ~= 0
    cd(['input_' model])
    tmin = ModellingInput(9,1);           % min thickness allowed - 10 plies
    tmax = ModellingInput(10,1);          % max thickness allowed - 150 plies
    constant.lam.tmin = tmin;                % saved for thickness normalisation
    constant.lam.tmax = tmax;                % saved for thickness normalisation
    [constant] = GenerateLaminateIniGuess(constant,Nlam,FLAGSPAR);    % Generate initial stacking sequence (data store in constant.lam.ini)
    cd ..
end

%=========================================================================%
% Mission data (GET BACK TO THIS - CHANGE REQ.'D TO GUST INPUT DATA)
%=========================================================================%
aero.rho    = 1.225;                    % air density
aero.mu = 1.7894e-5; % Viscosity of air at standard sea-level atmosphere
aero.Ksr = AircraftData(7,4); % Scale factor on the parasite drag for surface roughness and imperfections

GustData = xlsread(xlsFileName,'Gust');

gust.type     = GustData(1,1); % 1: regular 1-cos, 2: CS25 gust
gust.Uref     = GustData(2,1);
gust.Npergust = GustData(3,1);
gust.truncc   = GustData(4,1);
gust.H        = GustData(5,:);

%=========================================================================%
% Optimisation Settings
%=========================================================================%
% ---
if 1   % Optimiser settings
    opt.ObjType      = ModellingInput(26,1);
    
    % Number of constraints
    NEigval          = ModellingInput(29,1);   % # of Eigenvalues to consider
    opt.TrimConst    = ModellingInput(30,1);
    opt.BlendConst   = ModellingInput(31,1);
    opt.AileffConst  = ModellingInput(32,1);
    opt.BucklConst   = BUCKL;
    opt.oneGConst    = TWIST1G;
    
    if cross.numel ~= 0
        opt.constraint.LamFeasibility     = 6*length(constant.lam.ID);                                             % 1: 5 feasibility constraints per laminate -> 5*(#laminates)
    end
    
    opt.constraint.Eigval             = NEigval;
    
    if cross.numel ~= 0
        NStrainConstr = 0;
        for i=1:length(cross.lam)
            for j=1:length(cross.lam{i})
                NStrainConstr = NStrainConstr + 4*length(unique(cross.lam{i}{j}(cross.type{i}{j}==1)));
            end
        end
        opt.constraint.Strain             = 4*NStrainConstr;                                           % 2: tensile and shear strain at both end of each element -> 4*(#elements)
        
        NStrainCritConstr = 0;
        for i=1:length(cross.lam)
            for j=1:length(cross.lam{i})
                NStrainCritConstr = NStrainCritConstr + 8*length(unique(cross.lam{i}{j}(cross.type{i}{j}==1)));
            end
        end
        opt.constraint.StrainCrit             = NStrainCritConstr;
    
        if opt.BucklConst
            opt.constraint.BucklPerLam    = 8;              % Number of buckling panels to consider per laminate section
            opt.constraint.Buckl          = opt.constraint.BucklPerLam*8*length(constant.lam.ID); % BucklPerLam panels per laminate, for each laminate, two modes x two critical loads x two ends of the laminate
        else
            opt.constraint.Buckl          = 0;
        end
    end

    if opt.AileffConst
        opt.constraint.Aileff             = 1;
    else
        opt.constraint.Aileff             = 0;
    end
    
    if opt.TrimConst
        opt.constraint.Trim           = 2*(aero.nbs+1);      % 2 trim constraints (positive and negative) per aerodynamic cross-section (local angle of attack)
    else
        opt.constraint.Trim           = 0;
    end
    
    if opt.oneGConst
        opt.constraint.oneG            = 2*length(str.xyz(2:3:end));
        general.oneGmax                = ModellingInput(34,2);
    else
        opt.constraint.oneG            = 0;
    end
    
    if opt.BlendConst
        
        opt.BlendOptions.TopSkin       = true;
        opt.BlendOptions.BotSkin       = true;
        opt.BlendOptions.Spar          = true;
        opt.BlendOptions.SphereCoefVkA = 0.2;
        opt.BlendOptions.SphereCoefVkD = 0.2;
        
      
        NLamTop = numel(cell2mat(constant.lam.TopID(:)));
        NLamBot = numel(cell2mat(constant.lam.BotID(:)));
        opt.constraint.Blending = 2*(NLamTop*(NLamTop-1)/2 + NLamBot*(NLamBot-1)/2); % 5: Laminate Blending Constraints
        
        if FLAGSPAR
            SparNum     = cell2mat(constant.lam.SparNum(:)')';                  % correspond with LamSpar ID
            SparNumbers = unique(SparNum);
            Nspars      = length(SparNumbers);                                  % total # of spars
            Nconst      = 0;
            for iSpar = 1:Nspars
                NLamSpar = numel(find(SparNumbers(iSpar)==SparNum));
                Nconst = Nconst + 2*(NLamSpar*(NLamSpar-1)/2);
            end
            opt.constraint.Blending = opt.constraint.Blending + Nconst;
        end
        
    else
        opt.BlendOptions.TopSkin       = false;
        opt.BlendOptions.BotSkin       = false;
        opt.BlendOptions.Spar          = false;
        opt.constraint.Blending = 0;
    end
    
    if cross.numel ~= 0
        % Set init. guess
        opt.dvFull    = [];
        for ii=1:length(constant.lam.ID)
            curdir = cd;
            cd('../lam_param_laminate')
            LPs = 0*lam2lampar(constant.lam.ini.layup{ii});
            cd(curdir)
            opt.dvFull = [opt.dvFull; [LPs; constant.lam.ini.t(ii)]];
        end
    end
end
% ---

%=========================================================================%
% Output
%=========================================================================%
if 1    % Update Constant
    
    constant.general      = general;
    constant.str          = str;
    constant.aero         = aero;
    constant.gust         = gust;
    constant.cross        = cross;
    if opt.BucklConst
        constant.buckl    = buckl;
    end
    constant.fext         = fext;
    constant.lumped       = lumped;
    constant.opt          = opt;
    
    constant.inp.fixnodes = fixnodes;
    constant.inp.xyz      = reshape(xyz',[],1);
    constant.inp.theta    = theta;
    constant.inp.c        = c;
    constant.inp.xref     = xref;
    constant.inp.xyzshift = xyzshift;
    constant.inp.xyzqc    = reshape(xyzqc',[],1);
    if opt.oneGConst
    	constant.inp.y1g      = y1g;
    	constant.inp.theta1g  = theta1g;
    end
    
    constant.model = model;
   
    constant.general.weight    = ac_weight;
    constant.general.alpha_max = alpha_max;
    constant.general.MTOW = MTOW;
    constant.general.LW = LW;
    constant.general.TSFC = TSFC;
    constant.general.ZMO = ZMO;
    constant.general.MZFW = MZFW;
    constant.general.MLW = MLW;
    constant.general.IWW = IWW;
    constant.general.romflag = 0;
end

%=========================================================================%
% END
%=========================================================================%
