function [varargout] = analysis_loadcase(constant,dvFull,ANALYSISTYPE,TRIMDEF,GRAV,DERS,outflag,loadcase)

if constant.cross.numel ~= 0 
    dvLoad.tail = dvFull.tail;
else
    dvLoad = dvFull;
end

mainloopsize  = length(loadcase.nz);
% ---
if 1    % Pre-allocation && pre-prints
    if outflag == 1
        constantOut = cell(mainloopsize,1);
        statics     = cell(mainloopsize,1);
        strain      = cell(mainloopsize,1);
        if constant.opt.BucklConst
            buckl   = cell(mainloopsize,1);
        end
    else
        displc  = cell(mainloopsize,1);
        alphalc = cell(mainloopsize,1);  
        glc     = cell(mainloopsize,1);
        objlc   = cell(mainloopsize,1);
        if DERS == 1
            dglc          = cell(mainloopsize,1);
            dobjlc        = cell(mainloopsize,1);
        end
    end
    
%     cprintf([0 0 1],strcat(' outflag =',num2str(outflag),'\n' ))
    displayTable = [{'LoadCase #'}                            {'EAS'}        {'Mach'}                 {'Altitude'}             {'Load Factor'};
        num2cell([1:length(loadcase.EAS)]') num2cell(loadcase.EAS) num2cell(loadcase.M)   num2cell(loadcase.H)     num2cell(loadcase.nz)];
    for i=1:length(loadcase.fuel_level{1})
        dummat = [{['Fuel Tank ',num2str(i)]}];
        for j = 1:length(loadcase.fuel_level)
            dummat = [dummat;{loadcase.fuel_level{j}(i)}];
        end
        displayTable(:,end+1) = dummat;
    end
    clear loadcaseTemp
%     display(displayTable)
end
% ---

% Precompute cross-sectional analysis
if constant.cross.numel ~= 0
    [constant,lampar,stringer,~,crossmod,FEAS] = analysis_crossmod_preproc(constant,dvLoad.tail,ANALYSISTYPE,DERS);
else
    [constant,crossmod] = analysis_crossmod_preproc_beam(constant);
    FEAS = 1;
    lampar = [];
    stringer = [];
end

% Define lamination parameter constraints
if (ANALYSISTYPE == 1 || ANALYSISTYPE == 3) && constant.cross.numel ~= 0
    g.lam.g1   = -(lampar.g1);
    g.lam.g2   = -(lampar.g2);
    g.lam.g3   = -(lampar.g3);
    g.lam.g4   = -(lampar.g4);
    g.lam.g5   = -(lampar.g5);
    g.lam.g6   = -(lampar.g6);
end

if DERS == 1 && constant.cross.numel ~= 0 && (ANALYSISTYPE == 1 || ANALYSISTYPE == 3)
    dg.lam.g1.dtail   = -lampar.dg1ddv;
    dg.lam.g2.dtail   = -lampar.dg2ddv;
    dg.lam.g3.dtail   = -lampar.dg3ddv;
    dg.lam.g4.dtail   = -lampar.dg4ddv;
    dg.lam.g5.dtail   = -lampar.dg5ddv;
    dg.lam.g6.dtail   = -lampar.dg6ddv;
end

%% ---
for i=1:mainloopsize
    
    constant2 = constant;

    %% Select Load Case
    constant2.general.nz = loadcase.nz(i);
    constant2.aero.V     = loadcase.EAS(i);
    constant2.aero.M     = loadcase.M(i);
    constant2.aero.alpha0 = deg2rad(loadcase.alpha0(i));              % initial angle of attack, overwritten by trim
    constant2.general.TAS = eas2tas(loadcase.EAS(i),loadcase.H(i));
    
    if constant2.gust.type == 2
       constant2.gust.Uref = tas2eas(constant.gust.Uref,loadcase.H(i));
       constant2.gust.alt  = loadcase.H(i);
    end
    
    %% --- Add gravitational forces
    fext = constant2.fext;
    grav_cst = 9.81;
    Index = length(constant2.fext.type) + 1;
    fext.type{Index} = 'Non-structural masses';
    if GRAV  && isfield(constant2,'lumped')
        EntryNum = 1;
        for ilumped = 1 : length(constant2.lumped.type)
            for j = 1:length(constant2.lumped.mass{ilumped})
                fext.magnitude{Index}(EntryNum,:) = [0 0 -constant2.general.nz*constant2.lumped.mass{ilumped}(j)*grav_cst 0 0 0];
                fext.location{Index} (EntryNum,:)  = constant2.lumped.location{ilumped}(j,:);
                fext.follower{Index} (1,EntryNum)  = 0;
                fext.alphaflag{Index}(1,EntryNum)  = 1;
                EntryNum = EntryNum+1;
            end
        end
    end
    constant2.fext = fext;
    
    % Load Fuel Data (If Req'd)
    cd([constant.curdir,'/ext/PROTEUS/InputFiles'])

    if sum(loadcase.fuel_level{i})~=0
        [constant2] = AddFuel(constant2,loadcase.fuel_level{i},GRAV,grav_cst,'right');
        [constant2] = AddFuel(constant2,loadcase.fuel_level{i},GRAV,grav_cst,'left');
    end
    
    %% Store Inputs
    constant2.fext   = fext_inp(constant2.str,constant2.fext);
    if ~strcmp(constant.model,'EMB') 
        constant2.lumped = mext_inp(constant2.str,constant2.lumped);
    else
        constant2.lumped = mext_inp(constant2.str,constant2.lumped,ENGINE);
    end
    cd(constant.curdir)
    
    %% load previous results
    dispi = zeros(constant.str.Ndof,1);
    if loadcase.trim(i)
        alphai = 2*pi/180*constant2.general.nz;
    else
        alphai = constant2.aero.alpha0;
    end
    
    %% run analysis
    if outflag == 1
        if constant.general.romflag ~= 2
            if constant.opt.BucklConst
%                 [constantOut{i},statics{i},strain{i},buckl{i}] = analysis_crossmod(constant2,crossmod,dispi,alphai,ANALYSISTYPE,loadcase.trim(i),TRIMDEF,GRAV,DERS,outflag);
                [constantOut{i},statics{i}] = analysis_crossmod(constant2,crossmod,dispi,alphai,ANALYSISTYPE,loadcase.trim(i),TRIMDEF,GRAV,DERS,outflag);
            else
%                 [constantOut{i},statics{i},strain{i}] = analysis_crossmod(constant2,crossmod,dispi,alphai,ANALYSISTYPE,loadcase.trim(i),TRIMDEF,GRAV,DERS,outflag);
                [constantOut{i},statics{i}] = analysis_crossmod(constant2,crossmod,dispi,alphai,ANALYSISTYPE,loadcase.trim(i),TRIMDEF,GRAV,DERS,outflag);
            end
        else
            [constantOut{i},statics{i}] = analysis_crossmod(constant2,crossmod,dispi,alphai,ANALYSISTYPE,loadcase.trim(i),TRIMDEF,GRAV,DERS,outflag);
        end
    elseif DERS==1
        [displc{i},alphalc{i},objlc{i},glc{i},dobjlc{i},dglc{i}] = analysis_crossmod(constant2,crossmod,dispi,alphai,ANALYSISTYPE,loadcase.trim(i),TRIMDEF,GRAV,DERS,outflag);
        
    else
        [displc{i},alphalc{i},objlc{i},glc{i}] = analysis_crossmod(constant2,crossmod,dispi,alphai,ANALYSISTYPE,loadcase.trim(i),TRIMDEF,GRAV,DERS,outflag);
    end
end

%% Arrange constraints
if outflag == 0
    for i=1:mainloopsize
        obj.mass{i} = objlc{i}.mass;
        obj.lift{i} = objlc{i}.lift;
        obj.drag{i} = objlc{i}.drag;
        if constant.opt.ObjType == 2
           obj.range{i} = objlc{i}.range;
        end
        
        g.exmax{i} = glc{i}.exmax;
        g.exmin{i} = glc{i}.exmin;
        g.gammamax{i} = glc{i}.gammamax;
        g.strainrcrit{i} = glc{i}.strainrcrit;
        if constant.opt.constraint.Eigval>0
            g.aestab{i} = glc{i}.aestab;
        end
        if constant.opt.BucklConst
            g.buckl{i} = glc{i}.buckl;
        end
        if constant.opt.AileffConst
            g.aileff{i} = glc{i}.aileff;
        end
        if constant.opt.TrimConst
            g.alphal{i} = glc{i}.alphal;
        end
        if constant.opt.oneGConst
            g.twist{i} = glc{i}.twist;
        end
        if isfield(glc{i},'moen');
            g.moen{i} = glc{i}.moen;
        end
        if DERS == 1
            dobj.mass.dtail{i} = dobjlc{i}.mass.dtail;
            dobj.lift.dtail{i} = dobjlc{i}.lift.dtail;
            dobj.drag.dtail{i} = dobjlc{i}.drag.dtail;
            if constant.opt.ObjType == 2
                dobj.range.dtail{i} = dobjlc{i}.range.dtail;
            end

            dg.exmax.dtail{i} = dglc{i}.exmax.dtail;
            dg.exmin.dtail{i} = dglc{i}.exmin.dtail;
            dg.gammamax.dtail{i} = dglc{i}.gammamax.dtail;
            dg.strainrcrit.dtail{i} = dglc{i}.strainrcrit.dtail;
            if constant.opt.constraint.Eigval>0
                dg.aestab.dtail{i} = dglc{i}.aestab.dtail;
            end
            if constant.opt.BucklConst
                dg.buckl.dtail{i} = dglc{i}.buckl.dtail;
            end
            if constant.opt.AileffConst
                dg.aileff.dtail{i} = dglc{i}.aileff.dtail;
            end
            if constant.opt.TrimConst
                dg.alphal.dtail{i} = dglc{i}.alphal.dtail;
            end
            if constant.opt.oneGConst
                dg.twist.dtail{i} = dglc{i}.twist.dtail;
            end
            if isfield(glc{i},'moen');
                if isfield(dglc{i}.moen,'fold')
                    dg.moen.fold.tot.dtail{i} = dglc{i}.moen.fold.tot.dtail;
                    dg.moen.fold.pos.dtail{i} = dglc{i}.moen.fold.pos.dtail;
                    dg.moen.fold.neg.dtail{i} = dglc{i}.moen.fold.neg.dtail;
                end
                if isfield(dglc{i}.moen,'twist')
                    dg.moen.twist.tot.dtail{i} = dglc{i}.moen.twist.tot.dtail;
                    dg.moen.twist.pos.dtail{i} = dglc{i}.moen.twist.pos.dtail;
                    dg.moen.twist.neg.dtail{i} = dglc{i}.moen.twist.neg.dtail;
                end
                if isfield(dglc{i}.moen,'camber')
                    dg.moen.camber.tot.dtail{i} = dglc{i}.moen.camber.tot.dtail;
                    dg.moen.camber.pos.dtail{i} = dglc{i}.moen.camber.pos.dtail;
                    dg.moen.camber.neg.dtail{i} = dglc{i}.moen.camber.neg.dtail;
                end
                if isfield(dglc{i}.moen,'shear')
                    dg.moen.shear.tot.dtail{i} = dglc{i}.moen.shear.tot.dtail;
                    dg.moen.shear.pos.dtail{i} = dglc{i}.moen.shear.pos.dtail;
                    dg.moen.shear.neg.dtail{i} = dglc{i}.moen.shear.neg.dtail;
                end
            end
        end
    end
end

%% Output Results
if outflag == 1
    varargout{1} = constantOut;
    varargout{2} = lampar;
    varargout{3} = stringer;
    varargout{4} = crossmod;
    varargout{5} = statics;
    if constant.general.romflag ~= 2
        varargout{6} = strain;
        if constant.opt.BucklConst
            varargout{7} = buckl;
        end
    end
else
    varargout{1} = obj;
    varargout{2} = g;
    if DERS==1
        varargout{3} = dobj;
        varargout{4} = dg;
    end
end
