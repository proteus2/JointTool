function varargout = genStrModel(model)

curdir = cd;

% Specify number of computational threads, note should not be greater than
% specified in the pbs script ()
if isempty(gcp('nocreate'))
    parpool(2,'IdleTimeout', 300);
end

global GRAPHFLAG xlsFileName

addpath(genpath([cd,'/ext/PROTEUS/InBuildCommands']));

% --- Inputs
xlsFileName = ['input_',model,'/',model,'_Input.xlsx'];

cd ext/PROTEUS/InputFiles
[ModellingInput] = xlsread(xlsFileName,'ModellingInput');
cd(curdir)

LIN          = ModellingInput(12,1);                    % Flag for linear or non-linear, 1: Linear, 0: Non-linear
TRIMDEF      = ModellingInput(13,1);                    % trimdef: 1: weight = general.weight+str_mass+non_str_mass; 2: weight = general.weight
DERS         = ModellingInput(14,1);                    % set to 1 to calculate the derivative
GRAV         = ModellingInput(15,1);                    % Gravity Y/N

GRAPHFLAG    = 0;

OUTFLAG      = ModellingInput(18,1);                    % 1: for classic analysis, 0 for obj and constraint

ANALYSISTYPE = 1;

% Generate input
cd ext/PROTEUS/InputFiles
loadcase = GenerateLoadcases(xlsFileName);
loadcase.trim = 0;
constant = GenerateInput(model);
cd(curdir)

constant.opt.Optimiser = 'None';
constant.model   = model;
constant.curdir  = curdir;

% Prepare design variable vector
dvinp.tail = constant.opt.dvFull;

% Analysis or Optimisation
if OUTFLAG == 0
    if DERS == 1
        [obj,g,dobj,dg] = analysis_loadcase(constant,dvinp,ANALYSISTYPE,LIN,TRIMDEF,GRAV,DERS,OUTFLAG,loadcase);
        [f0val,fval,df0dx,dfdx,ConstraintLabel] = PostProc(constant,loadcase,obj,g,dobj,dg);
    else
        [obj,g] = analysis_loadcase(constant,dvinp,ANALYSISTYPE,LIN,TRIMDEF,GRAV,DERS,OUTFLAG,loadcase);
        [f0val,fval,~,~,ConstraintLabel] = PostProc(constant,loadcase,obj,g);
    end
elseif OUTFLAG == 1
    if constant.opt.BucklConst
        [constantOut,lampar,stringer,crossmod,statics,strain,buckl] = analysis_loadcase(constant,dvinp,ANALYSISTYPE,LIN,TRIMDEF,GRAV,DERS,OUTFLAG,loadcase);
    else
        [constantOut,lampar,stringer,crossmod,statics,strain] = analysis_loadcase(constant,dvinp,ANALYSISTYPE,LIN,TRIMDEF,GRAV,DERS,OUTFLAG,loadcase);
    end
end

% Output
if OUTFLAG==0
    if DERS==1
        varargout{1} = f0val;
        varargout{2} = fval;
        varargout{3} = df0dx;
        varargout{4} = dfdx;
        varargout{5} = ConstraintLabel;
    else
        varargout{1} = f0val;
        varargout{2} = fval;
        varargout{3} = ConstraintLabel;
    end
elseif OUTFLAG==1
    varargout{1} = constantOut;
    varargout{2} = lampar; 
    varargout{3} = stringer;
    varargout{4} = crossmod;
    varargout{5} = statics;
    varargout{6} = strain;
    if constant.opt.BucklConst
        varargout{7} = buckl;
    end 
end
end
