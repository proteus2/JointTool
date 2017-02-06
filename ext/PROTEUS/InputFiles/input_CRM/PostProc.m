function [varargout] = PostProc(constant,loadcase,varargin)

%%% Output:
%%% varargout{1}: f0val, objective
%%% varargout{2}: fval, constraint vector
%%% varargout{3}: df0dx, derivate of objective wrt xval (in case gradients)
%%% varargout{4}: dfdx, derivative of constraint wrt xval (in case gradients)
%%% varargout{5}: ConstraintLabels, constraint name

%%% Display functionality:
%%% Display per loadcase, critical constraint value for each constraint

%%% NOTE: DO NOT FORGET TO NORMALIZE OBJECTIVE AND CONSTRAINTS

if length(varargin) == 5
    obj = varargin{1};
    g = varargin{2};
    dobj = varargin{3};
    dg = varargin{4};
    normal = varargin{5};
    ders = 1;
    normflag = 1;
elseif length(varargin) == 4
    obj = varargin{1};
    g = varargin{2};
    dobj = varargin{3};
    dg = varargin{4};
    ders = 1;
    normflag = 0;
elseif length(varargin) == 3
    obj = varargin{1};
    g = varargin{2};
    normal = varargin{3};
    ders = 0;
    normflag = 1;
elseif length(varargin) == 2
    obj = varargin{1};
    g = varargin{2};
    ders = 0;
    normflag = 0;
end

%% Example only tailoring

Nlam = length(constant.lam.ID);
numlc = length(loadcase.EAS);

f0val = -obj;
if ders == 1
   df0dx = -dobj.dtail; 
end

if normflag == 1
   f0val = f0val/abs(normal.obj);
   if ders == 1
       df0dx = df0dx/abs(normal.obj);
   end
end

fval = [g.lam.g1;g.lam.g2;g.lam.g3;g.lam.g4;g.lam.g5;g.lam.g6];
if ders == 1
    dfvaldx = [dg.lam.g1.dtail;dg.lam.g2.dtail;dg.lam.g3.dtail;dg.lam.g4.dtail;dg.lam.g5.dtail;dg.lam.g6.dtail];
end


constraintvec = mat2cell((1:Nlam)',ones(Nlam,1),1);
constraintvec = cellfun(@(x)num2str(x),constraintvec,'UniformOutput',0);
ConstraintLabel = [strcat('G1: lam. ',constraintvec);strcat('G2: lam. ',constraintvec);strcat('G3: lam. ',constraintvec);strcat('G4: lam. ',constraintvec);strcat('G5: lam. ',constraintvec);strcat('G6: lam. ',constraintvec)];

for i = 1:numlc
    if loadcase.gustflag(i)==1
        numgust = length(constant.gust.H);
    else
        numgust = 1;
    end
    
    % Exmax
    fval = [fval;(g.exmax{i}-constant.mat.exmax)/constant.mat.exmax];
    if ders == 1
        dfvaldx = [dfvaldx;dg.exmax.dtail{i}/constant.mat.exmax];
    end
    ConstraintEx = repmat({['Strain, exmax, LC ',num2str(i)]},size(g.exmax{i},1),1);
    if loadcase.gustflag(i)==1
        nex = size(g.exmax{i},1)/numgust;
        for j=1:numgust
            ConstraintEx((j-1)*nex+(1:nex)) = strcat(ConstraintEx((j-1)*nex+(1:nex)),', Gust ',num2str(j));
        end
    end
    ConstraintLabel = [ConstraintLabel;ConstraintEx];
    
    % Exmin
    fval = [fval;(g.exmin{i}-constant.mat.exmin)/constant.mat.exmin];
    if ders == 1
        dfvaldx = [dfvaldx;dg.exmin.dtail{i}/constant.mat.exmin];
    end
    ConstraintEx = repmat({['Strain, exmin, LC ',num2str(i)]},size(g.exmin{i},1),1);
    if loadcase.gustflag(i)==1
        nex = size(g.exmin{i},1)/numgust;
        for j=1:numgust
            ConstraintEx((j-1)*nex+(1:nex)) = strcat(ConstraintEx((j-1)*nex+(1:nex)),', Gust ',num2str(j));
        end
    end
    ConstraintLabel = [ConstraintLabel;ConstraintEx];
    
    % Gammamax
    fval = [fval;(g.gammamax{i}-constant.mat.gmax)/constant.mat.gmax];
    if ders == 1
        dfvaldx = [dfvaldx;dg.gammamax.dtail{i}/constant.mat.gmax];
    end
    ConstraintGamma = repmat({['Strain, gammamax, LC ',num2str(i)]},size(g.gammamax{i},1),1);
    if loadcase.gustflag(i)==1
        ngamma = size(g.gammamax{i},1)/numgust;
        for j=1:numgust
            ConstraintGamma((j-1)*nex+(1:ngamma)) = strcat(ConstraintGamma((j-1)*ngamma+(1:ngamma)),', Gust ',num2str(j));
        end
    end
    ConstraintLabel = [ConstraintLabel;ConstraintGamma];
    
    % Aeroelastic stability
    if constant.opt.constraint.Eigval>0
        fval = [fval;g.aestab{i}];
        if ders == 1
            dfvaldx = [dfvaldx;dg.aestab.dtail{i}];
        end
        ConstraintAestab = repmat({['Aeroelastic stability, LC ',num2str(i)]},size(g.aestab{i},1),1);
        ConstraintLabel = [ConstraintLabel;ConstraintAestab];
    end
    
    % Buckling
    if constant.opt.BucklConst
        fval = [fval;g.buckl{i}-1];
        if ders == 1
            dfvaldx = [dfvaldx;dg.buckl.dtail{i}];
        end
        if loadcase.gustflag(i)==1
            ConstraintBuckl = reshape((repmat(constraintvec,1,6*numgust*constant.opt.constraint.BucklPerLam*8))',[],1);
            ConstraintBuckl = strcat('Buckling, lam.',ConstraintBuckl,', LC',num2str(i));
            nbuckl = 6*constant.opt.constraint.BucklPerLam*8*Nlam;
            for j=1:numgust
                ConstraintBuckl((j-1)*nbuckl+(1:nbuckl)) = strcat(ConstraintBuckl((j-1)*nbuckl+(1:nbuckl)),', Gust ',num2str(j));
            end
        else
            ConstraintBuckl = reshape((repmat(constraintvec,1,numgust*constant.opt.constraint.BucklPerLam*8))',[],1);
            ConstraintBuckl = strcat('Buckling, lam.',ConstraintBuckl,', LC',num2str(i));
        end
        ConstraintLabel = [ConstraintLabel;ConstraintBuckl];
    end
    
    
    % Aileron effectiveness
    if constant.opt.AileffConst
        fval = [fval;(constant.aileron.mineff-g.aileff{i})/constant.aileron.mineff];
        if ders == 1
            dfvaldx = [dfvaldx;-dg.aileff.dtail{i}/constant.aileron.mineff];
        end
        ConstraintLabel = [ConstraintLabel;['Aileron effectiveness, LC ',num2str(i)]];    
    end
    
    
    % Local angle of attack
    if constant.opt.TrimConst
        fval = [fval;(g.alphal{i}-constant.general.alpha_max)/constant.general.alpha_max]; % Maximum angle of attack
        fval = [fval;-(g.alphal{i}+constant.general.alpha_max)/constant.general.alpha_max]; % Minimum angle of attack
        if ders == 1
            dfvaldx = [dfvaldx;dg.alphal.dtail{i}/constant.general.alpha_max];
        end
        if ders == 1
            dfvaldx = [dfvaldx;-dg.alphal.dtail{i}/constant.general.alpha_max];
        end
        ConstraintAlphal = repmat({['Max. local angle of attack, LC ',num2str(i)]},size(g.alphal{i},1),1);
        ConstraintAlphal = [ConstraintAlphal;repmat({['Min. local angle of attack, LC ',num2str(i)]},size(g.alphal{i},1),1)];
        ConstraintLabel = [ConstraintLabel;ConstraintAlphal];
    end
    
    % 1g Twist
    if i == 1
        if constant.opt.oneGConst
            fval = [fval;(g.twist{i}-(constant.str.theta1g+constant.general.oneGmax))/10]; % Maximum overshoot in twist
            fval = [fval;((constant.str.theta1g-constant.general.oneGmax)-g.twist{i})/10]; % Minimum undershoot in twist
            if ders == 1
                dfvaldx = [dfvaldx;dg.twist.dtail{i}/10];
            end
            if ders == 1
                dfvaldx = [dfvaldx;-dg.twist.dtail{i}/10];
            end
            ConstraintTwist = repmat({'1g Twist upperbound'},size(g.twist{i},1),1);
            ConstraintTwist = [ConstraintTwist;repmat({'1g Twist lowerbound'},size(g.twist{i},1),1)];
            ConstraintLabel = [ConstraintLabel;ConstraintTwist];
        end
    end
end

% Normalize sensitivities
if ders == 1
    df0dx = (constant.opt.NormMatrix\df0dx');
    dfvaldx   = (constant.opt.NormMatrix\dfvaldx')';
end

varargout{1} = f0val;
varargout{2} = fval;
if ders == 1
    varargout{3} = df0dx;
    varargout{4} = dfvaldx;
end
varargout{5} = ConstraintLabel;
%}

%% Example tailoring and morphing
%{
Nlam = length(constant.lam.ID);
numlc = length(loadcase.EAS);

ntail = 9*Nlam;
ntwist = size(constant.morph.inp.twist.anglefin,1); % Number of twist design variables per loadcase
nshear = size(constant.morph.inp.shear.anglefin,1); % Number of shear design variables per loadcase
nfold = size(constant.morph.inp.fold.anglefin,1); % Number of fold design variables per loadcase
ncamber = size(constant.morph.inp.camber.paramfin,1); % Number of camber design variables per loadcase
nspan = 1;

nmorph = (ntwist+nshear+nfold+ncamber+nspan);

ndv = ntail+numlc*nmorph;

f0val = obj;
if ders == 1
   df0dx = sparse(1,ndv);
   df0dx(1,1:ntail) = dobj.dtail;
   % In this case range is given for a single loadcase, so only dependent
   % upon that loadcase
%    for i=1:numlc
   df0dx(1,ntail+(1:ntwist)) = dobj.dphifin*constant.morph.twist.dangleddv;
   df0dx(1,ntail+ntwist+(1:nshear)) = dobj.dpsifin*constant.morph.shear.dangleddv;
   df0dx(1,ntail+ntwist+nshear+(1:nfold)) = dobj.dthetafin*constant.morph.fold.dangleddv;
   df0dx(1,ntail+ntwist+nshear+nfold+(1:ncamber)) = dobj.dparamfin*constant.morph.camber.dparamddv;
   df0dx(1,ntail+ntwist+nshear+nfold+ncamber+(1:nspan)) = dobj.dextfin*constant.morph.span.dextddv;
%    end
end

if normflag == 1
   f0val = f0val/normal.obj;
   if ders == 1
       df0dx = df0dx/normal.obj;
   end
end

fval = [g.lam.g1;g.lam.g2;g.lam.g3;g.lam.g4;g.lam.g5;g.lam.g6];
if ders == 1
    dfvaldx = [dg.lam.g1.dtail;dg.lam.g2.dtail;dg.lam.g3.dtail;dg.lam.g4.dtail;dg.lam.g5.dtail;dg.lam.g6.dtail];
    dfvaldx(length(fval),ntail+(1:numlc*nmorph)) = 0;
end


constraintvec = mat2cell((1:Nlam)',ones(Nlam,1),1);
constraintvec = cellfun(@(x)num2str(x),constraintvec,'UniformOutput',0);
ConstraintLabel = [strcat('G1: lam. ',constraintvec);strcat('G2: lam. ',constraintvec);strcat('G3: lam. ',constraintvec);strcat('G4: lam. ',constraintvec);strcat('G5: lam. ',constraintvec);strcat('G6: lam. ',constraintvec)];

for i = 1:numlc
    if loadcase.gustflag(i)==1
        numgust = length(constant.gust.H);
    else
        numgust = 1;
    end
    
    % Exmax
    fval = [fval;(g.exmax{i}-constant.mat.exmax)/constant.mat.exmax];
    if ders == 1
        dfvaldx = [dfvaldx;[dg.exmax.dtail{i},zeros(length(g.exmax{i}),(i-1)*nmorph),...
            dg.exmax.dphifin{i}*constant.morph.twist.dangleddv,...
            dg.exmax.dpsifin{i}*constant.morph.shear.dangleddv,...
            dg.exmax.dthetafin{i}*constant.morph.fold.dangleddv,...
            dg.exmax.dparamfin{i}*constant.morph.camber.dparamddv,...
            dg.exmax.dextfin{i}*constant.morph.span.dextddv,...
            zeros(length(g.exmax{i}),(numlc-i)*nmorph)]/constant.mat.exmax];
    end
    ConstraintEx = repmat({['Strain, exmax, LC ',num2str(i)]},size(g.exmax{i},1),1);
    if loadcase.gustflag(i)==1
        nex = size(g.exmax{i},1)/numgust;
        for j=1:numgust
            ConstraintEx((j-1)*nex+(1:nex)) = strcat(ConstraintEx((j-1)*nex+(1:nex)),', Gust ',num2str(j));
        end
    end
    ConstraintLabel = [ConstraintLabel;ConstraintEx];
    
    % Exmin
    fval = [fval;(g.exmin{i}-constant.mat.exmin)/constant.mat.exmin];
    if ders == 1
        dfvaldx = [dfvaldx;[dg.exmin.dtail{i},zeros(length(g.exmin{i}),(i-1)*nmorph),...
            dg.exmin.dphifin{i}*constant.morph.twist.dangleddv,...
            dg.exmin.dpsifin{i}*constant.morph.shear.dangleddv,...
            dg.exmin.dthetafin{i}*constant.morph.fold.dangleddv,...
            dg.exmin.dparamfin{i}*constant.morph.camber.dparamddv,...
            dg.exmin.dextfin{i}*constant.morph.span.dextddv,...
            zeros(length(g.exmin{i}),(numlc-i)*nmorph)]/constant.mat.exmin];
    end
    ConstraintEx = repmat({['Strain, exmin, LC ',num2str(i)]},size(g.exmin{i},1),1);
    if loadcase.gustflag(i)==1
        nex = size(g.exmin{i},1)/numgust;
        for j=1:numgust
            ConstraintEx((j-1)*nex+(1:nex)) = strcat(ConstraintEx((j-1)*nex+(1:nex)),', Gust ',num2str(j));
        end
    end
    ConstraintLabel = [ConstraintLabel;ConstraintEx];
    
    % Gammamax
    fval = [fval;(g.gammamax{i}-constant.mat.gmax)/constant.mat.gmax];
    if ders == 1
        dfvaldx = [dfvaldx;[dg.gammamax.dtail{i},zeros(length(g.gammamax{i}),(i-1)*nmorph),...
            dg.gammamax.dphifin{i}*constant.morph.twist.dangleddv,...
            dg.gammamax.dpsifin{i}*constant.morph.shear.dangleddv,...
            dg.gammamax.dthetafin{i}*constant.morph.fold.dangleddv,...
            dg.gammamax.dparamfin{i}*constant.morph.camber.dparamddv,...
            dg.gammamax.dextfin{i}*constant.morph.span.dextddv,...
            zeros(length(g.gammamax{i}),(numlc-i)*nmorph)]/constant.mat.gmax];
    end
    ConstraintGamma = repmat({['Strain, gammamax, LC ',num2str(i)]},size(g.gammamax{i},1),1);
    if loadcase.gustflag(i)==1
        ngamma = size(g.gammamax{i},1)/numgust;
        for j=1:numgust
            ConstraintGamma((j-1)*ngamma+(1:ngamma)) = strcat(ConstraintGamma((j-1)*ngamma+(1:ngamma)),', Gust ',num2str(j));
        end
    end
    ConstraintLabel = [ConstraintLabel;ConstraintGamma];
    
    % Rcrit
    fval = [fval;(g.strainrcrit{i}-1)];
    if ders == 1
        dfvaldx = [dfvaldx;[dg.strainrcrit.dtail{i},zeros(length(g.strainrcrit{i}),(i-1)*nmorph),...
            dg.strainrcrit.dphifin{i}*constant.morph.twist.dangleddv,...
            dg.strainrcrit.dpsifin{i}*constant.morph.shear.dangleddv,...
            dg.strainrcrit.dthetafin{i}*constant.morph.fold.dangleddv,...
            dg.strainrcrit.dparamfin{i}*constant.morph.camber.dparamddv,...
            dg.strainrcrit.dextfin{i}*constant.morph.span.dextddv,...
            zeros(length(g.strainrcrit{i}),(numlc-i)*nmorph)]];
    end
    ConstraintStrainR = repmat({['Strain, rcrit, LC ',num2str(i)]},size(g.strainrcrit{i},1),1);
    if loadcase.gustflag(i)==1
        ngamma = size(g.strainrcrit{i},1)/numgust;
        for j=1:numgust
            ConstraintStrainR((j-1)*ngamma+(1:ngamma)) = strcat(ConstraintStrainR((j-1)*ngamma+(1:ngamma)),', Gust ',num2str(j));
        end
    end
    ConstraintLabel = [ConstraintLabel;ConstraintStrainR];
    
    % Aeroelastic stability
    if constant.opt.constraint.Eigval>0
        fval = [fval;g.aestab{i}];
        if ders == 1
            dfvaldx = [dfvaldx;[dg.aestab.dtail{i},zeros(length(g.aestab{i}),(i-1)*nmorph),...
                dg.aestab.dphifin{i}*constant.morph.twist.dangleddv,...
                dg.aestab.dpsifin{i}*constant.morph.shear.dangleddv,...
                dg.aestab.dthetafin{i}*constant.morph.fold.dangleddv,...
                dg.aestab.dparamfin{i}*constant.morph.camber.dparamddv,...
                dg.aestab.dextfin{i}*constant.morph.span.dextddv,...
                zeros(length(g.aestab{i}),(numlc-i)*nmorph)]];
        end
        ConstraintAestab = repmat({['Aeroelastic stability, LC ',num2str(i)]},size(g.aestab{i},1),1);
        ConstraintLabel = [ConstraintLabel;ConstraintAestab];
    end
    
    % Buckling
    if constant.opt.BucklConst
        fval = [fval;g.buckl{i}-1];
        if ders == 1
            dfvaldx = [dfvaldx;[dg.buckl.dtail{i},zeros(length(g.buckl{i}),(i-1)*nmorph),...
                dg.buckl.dphifin{i}*constant.morph.twist.dangleddv,...
                dg.buckl.dpsifin{i}*constant.morph.shear.dangleddv,...
                dg.buckl.dthetafin{i}*constant.morph.fold.dangleddv,...
                dg.buckl.dparamfin{i}*constant.morph.camber.dparamddv,...
                dg.buckl.dextfin{i}*constant.morph.span.dextddv,...
                zeros(length(g.buckl{i}),(numlc-i)*nmorph)]];
        end
        if loadcase.gustflag(i)==1
            ConstraintBuckl = reshape((repmat(constraintvec,1,6*numgust*constant.opt.constraint.BucklPerLam*8))',[],1);
            ConstraintBuckl = strcat('Buckling, lam.',ConstraintBuckl,', LC',num2str(i));
            nbuckl = 6*constant.opt.constraint.BucklPerLam*8*Nlam;
            for j=1:numgust
                ConstraintBuckl((j-1)*nbuckl+(1:nbuckl)) = strcat(ConstraintBuckl((j-1)*nbuckl+(1:nbuckl)),', Gust ',num2str(j));
            end
        else
            ConstraintBuckl = reshape((repmat(constraintvec,1,numgust*constant.opt.constraint.BucklPerLam*8))',[],1);
            ConstraintBuckl = strcat('Buckling, lam.',ConstraintBuckl,', LC',num2str(i));
        end
        ConstraintLabel = [ConstraintLabel;ConstraintBuckl];
    end

    % Aileron effectiveness
    if constant.opt.AileffConst
        fval = [fval;(constant.aileron.mineff-g.aileff{i})/constant.aileron.mineff];
        if ders == 1
            dfvaldx = [dfvaldx;-[dg.aileff.dtail{i},zeros(length(g.aileff{i}),(i-1)*nmorph),...
                dg.aileff.dphifin{i}*constant.morph.twist.dangleddv,...
                dg.aileff.dpsifin{i}*constant.morph.shear.dangleddv,...
                dg.aileff.dthetafin{i}*constant.morph.fold.dangleddv,...
                dg.aileff.dparamfin{i}*constant.morph.camber.dparamddv,...
                dg.aileff.dextfin{i}*constant.morph.span.dextddv,...
                zeros(length(g.aileff{i}),(numlc-i)*nmorph)]/constant.aileron.mineff];
        end
    	ConstraintLabel = [ConstraintLabel;['Aileron effectiveness, LC ',num2str(i)]];    
    end
    
    % Local angle of attack
    if constant.opt.TrimConst
        fval = [fval;(g.alphal{i}-constant.general.alpha_max)/constant.general.alpha_max]; % Maximum angle of attack
        fval = [fval;-(g.alphal{i}+constant.general.alpha_max)/constant.general.alpha_max]; % Minimum angle of attack
        if ders == 1
            dfvaldx = [dfvaldx;[dg.alphal.dtail{i},zeros(length(g.alphal{i}),(i-1)*nmorph),...
                dg.alphal.dphifin{i}*constant.morph.twist.dangleddv,...
                dg.alphal.dpsifin{i}*constant.morph.shear.dangleddv,...
                dg.alphal.dthetafin{i}*constant.morph.fold.dangleddv,...
                dg.alphal.dparamfin{i}*constant.morph.camber.dparamddv,...
                dg.alphal.dextfin{i}*constant.morph.span.dextddv,...
                zeros(length(g.alphal{i}),(numlc-i)*nmorph)]/constant.general.alpha_max];
        end
        if ders == 1
            dfvaldx = [dfvaldx;-[dg.alphal.dtail{i},zeros(length(g.alphal{i}),(i-1)*nmorph),...
                dg.alphal.dphifin{i}*constant.morph.twist.dangleddv,...
                dg.alphal.dpsifin{i}*constant.morph.shear.dangleddv,...
                dg.alphal.dthetafin{i}*constant.morph.fold.dangleddv,...
                dg.alphal.dparamfin{i}*constant.morph.camber.dparamddv,...
                dg.alphal.dextfin{i}*constant.morph.span.dextddv,...
                zeros(length(g.alphal{i}),(numlc-i)*nmorph)]/constant.general.alpha_max];
        end
        ConstraintAlphal = repmat({['Max. local angle of attack, LC ',num2str(i)]},size(g.alphal{i},1),1);
        ConstraintAlphal = [ConstraintAlphal;repmat({['Min. local angle of attack, LC ',num2str(i)]},size(g.alphal{i},1),1)];
        ConstraintLabel = [ConstraintLabel;ConstraintAlphal];
    end
    
    % 1g Twist
    if i == 1
        if constant.opt.oneGConst
            fval = [fval;(g.twist{i}-(constant.str.theta1g+constant.general.oneGmax))/10]; % Maximum overshoot in twist
            fval = [fval;((constant.str.theta1g-constant.general.oneGmax)-g.twist{i})/10]; % Minimum undershoot in twist
            if ders == 1
                dfvaldx = [dfvaldx;[dg.twist.dtail{i},zeros(length(g.twist{i}),(i-1)*nmorph),...
                dg.twist.dphifin{i}*constant.morph.twist.dangleddv,...
                dg.twist.dpsifin{i}*constant.morph.shear.dangleddv,...
                dg.twist.dthetafin{i}*constant.morph.fold.dangleddv,...
                dg.twist.dparamfin{i}*constant.morph.camber.dparamddv,...
                dg.twist.dextfin{i}*constant.morph.span.dextddv,...
                zeros(length(g.twist{i}),(numlc-i)*nmorph)]/10];
            end
            if ders == 1
                dfvaldx = [dfvaldx;-[dg.twist.dtail{i},zeros(length(g.twist{i}),(i-1)*nmorph),...
                dg.twist.dphifin{i}*constant.morph.twist.dangleddv,...
                dg.twist.dpsifin{i}*constant.morph.shear.dangleddv,...
                dg.twist.dthetafin{i}*constant.morph.fold.dangleddv,...
                dg.twist.dparamfin{i}*constant.morph.camber.dparamddv,...
                dg.twist.dextfin{i}*constant.morph.span.dextddv,...
                zeros(length(g.twist{i}),(numlc-i)*nmorph)]/10];
            end
            ConstraintTwist = repmat({'1g Twist upperbound'},size(g.twist{i},1),1);
            ConstraintTwist = [ConstraintTwist;repmat({'1g Twist lowerbound'},size(g.twist{i},1),1)];
            ConstraintLabel = [ConstraintLabel;ConstraintTwist];
        end
    end
end

% Normalize sensitivities
if ders == 1
    df0dx = (constant.opt.NormMatrix\df0dx');
    dfvaldx   = (constant.opt.NormMatrix\dfvaldx')';
end

varargout{1} = f0val;
varargout{2} = fval;
if ders == 1
    varargout{3} = df0dx;
    varargout{4} = dfvaldx;
end
varargout{5} = ConstraintLabel;
%}