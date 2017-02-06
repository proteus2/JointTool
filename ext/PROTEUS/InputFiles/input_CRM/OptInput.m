function [ndv,xval,xmin,xmax,dvLabels,nconstr,constant] = OptInput(constant,loadcase)

global FLAGSPAR

%%% Write function such that:
%%% ndv: total number of design variables
%%% xval: Initial design variable vector
%%% xmin: lower bounds of design variables
%%% xmax: upper bounds of design variables
%%% dvLabels: cell vector storing design variable labels
%%% nconstr: total number of constraints
%%% constant.opt.NormMatrix: normalisation matrix for design variables
%%% constant.opt.NormFct: normalisation function for design variables
%%% constant.opt.DeNormFct: de-normalisation function for design variables

%% Example only tailoring:

Nlam = length(constant.lam.ID);
numlc = length(loadcase.EAS);
numgust = sum(loadcase.gustflag);

% Number of design variables
ndv = 9*Nlam; % In case of tailoring and thickness

% Initial guess
xval = constant.opt.dvFull;

% Design variable limits
xmin = repmat([-1*ones(8,1);constant.lam.tmin],Nlam,1);
xmax = repmat([1*ones(8,1);constant.lam.tmax],Nlam,1);

% Normalising Matrix
normvec = 2./(xmax-xmin);
NormMatrix = diag(normvec);
dv0 = (xmin+xmax)/2; % Centering Vector

NormFct   = @(dv)   NormMatrix*(dv-dv0);
DeNormFct = @(xval) NormMatrix\xval + dv0;

constant.opt.NormMatrix = NormMatrix;
constant.opt.NormFct    = NormFct;       % save normalisation fct
constant.opt.DeNormFct  = DeNormFct;     % save de-normalisation fct

% Normalise bounds and initial guess
xmin = NormFct(xmin);
xmax = NormFct(xmax);
xval = NormFct(xval);

% Design variable labels
dvLabels = cell(9*Nlam,1);
for ii = 1 : length(constant.lam.TopID)
    % Top skin
    for j = 1:length(constant.lam.TopID{ii})
        lamnum = constant.lam.TopID{ii}(j);
        genstr = ['Top skin: Lam ',num2str(lamnum),', '];
        dvLabels(9*(lamnum-1)+(1:8)) = mat2cell([repmat(genstr,8,1),['V1A';'V2A';'V3A';'V4A';'V1D';'V2D';'V3D';'V4D']],ones(8,1),(length(genstr)+3));
        dvLabels{9*(lamnum-1)+9} = [genstr,'t'];
    end
    
    % Bottom skin
    for j = 1:length(constant.lam.BotID{ii})
        lamnum = constant.lam.BotID{ii}(j);
        genstr = ['Bottom skin: Lam ',num2str(lamnum),', '];
        dvLabels(9*(lamnum-1)+(1:8),1) = mat2cell([repmat(genstr,8,1),['V1A';'V2A';'V3A';'V4A';'V1D';'V2D';'V3D';'V4D']],ones(8,1),(length(genstr)+3));
        dvLabels{9*(lamnum-1)+9} = [genstr,'t'];
    end
    
    % Spars
    if FLAGSPAR
        for j = 1:length(constant.lam.SparID{ii})
            lamnum = constant.lam.SparID{ii}(j);
            genstr = ['Spars: Lam ',num2str(lamnum),', '];
            dvLabels(9*(lamnum-1)+(1:8),1) = mat2cell([repmat(genstr,8,1),['V1A';'V2A';'V3A';'V4A';'V1D';'V2D';'V3D';'V4D']],ones(8,1),(length(genstr)+3));
            dvLabels{9*(lamnum-1)+9} = [genstr,'t'];
        end
    end
end

% Number of constraints
nconstr = constant.opt.constraint.LamFeasibility+...                        % Lamination parameters feasibility constraints
    +constant.opt.constraint.Blending+...                                   % Blending constraint, please verify!!!!
    +numlc*constant.opt.constraint.Eigval+...                               % Aeroelastic stability per loadcase
    +constant.opt.constraint.oneG+...                                       % 1g trim constraint
    +numlc*constant.opt.constraint.Trim+...                                 % Maximum and minimum angle of attack per loadcase
    +numlc*constant.opt.constraint.Aileff+...                               % Aileron effectiveness per loadcase
    +(numlc-numgust)*constant.opt.constraint.Strain+...                     % Strain constraints for non-gust loadcases
    +2*numgust*length(constant.gust.H)*constant.opt.constraint.Strain+...   % Strain constraints for gust loadcases (per gust length: positive, negative, min, max deflection)
    +(numlc-numgust)*constant.opt.constraint.Buckl+...                      % Buckling constraints for non-gust loadcases
    +2*numgust*length(constant.gust.H)*constant.opt.constraint.Buckl;       % Buckling constraints for gust loadcases (per gust length: positive, negative, min, max deflection)

%}
%% Example tailoring and morphing, where only the final morphing shape is a
% design variable, the initial morphing shape depends on other loadcases
%{
Nlam = length(constant.lam.ID);
numlc = length(loadcase.EAS);
numgust = sum(loadcase.gustflag);

ndv = 9*Nlam+... % Lamination parameters and thickness
    +numel(constant.morph.inp.twist.anglefin)+... % Twist morphing
    +numel(constant.morph.inp.shear.anglefin)+... % Shear morphing
    +numel(constant.morph.inp.fold.anglefin)+...  % Fold morphing
    +numel(constant.morph.inp.camber.paramfin)+...% Camber morphing
    +numel(constant.morph.inp.span.extfin);       % Span morphing

% Initial guess
xval = constant.opt.dvFull;

for i=1:numlc
    xval = [xval;
        deg2rad(constant.morph.inp.twist.anglefin(:,i));
        deg2rad(constant.morph.inp.shear.anglefin(:,i));
        deg2rad(constant.morph.inp.fold.anglefin(:,i));
        constant.morph.inp.camber.paramfin(:,i);
        constant.morph.inp.span.extfin(:,i)];
end

% Design variable limits
xmin = repmat([-1*ones(8,1);constant.lam.tmin],Nlam,1);     % Lamination parameters and thickness
xmax = repmat([1*ones(8,1);constant.lam.tmax],Nlam,1);      % Lamination parameters and thickness
xminmorph = repmat([deg2rad(constant.morph.inp.twist.low)';...
    deg2rad(constant.morph.inp.shear.low)';...
    deg2rad(constant.morph.inp.fold.low)';...
    constant.morph.inp.camber.low';...
    constant.morph.inp.span.low'],numlc,1); 
xmaxmorph = repmat([deg2rad(constant.morph.inp.twist.high)';...
    deg2rad(constant.morph.inp.shear.high)';...
    deg2rad(constant.morph.inp.fold.high)';...
    constant.morph.inp.camber.high';...
    constant.morph.inp.span.high'],numlc,1); 
xmin = [xmin;xminmorph];
xmax = [xmax;xmaxmorph];

% Normalising Matrix
normvec = 2./(xmax-xmin);
NormMatrix = diag(normvec);
dv0 = (xmin+xmax)/2; % Centering Vector

NormFct   = @(dv)   NormMatrix*(dv-dv0);
DeNormFct = @(xval) NormMatrix\xval + dv0;

constant.opt.NormMatrix = NormMatrix;
constant.opt.NormFct    = NormFct;       % save normalisation fct
constant.opt.DeNormFct  = DeNormFct;     % save de-normalisation fct

% Normalise bounds and initial guess
xmin = NormFct(xmin);
xmax = NormFct(xmax);
xval = NormFct(xval);

% Design variable labels
dvLabels = cell(ndv,1);
for ii = 1 : length(constant.lam.TopID)
    % Top skin
    for j = 1:length(constant.lam.TopID{ii})
        lamnum = constant.lam.TopID{ii}(j);
        genstr = ['Top skin: Lam ',num2str(lamnum),', '];
        dvLabels(9*(lamnum-1)+(1:8)) = mat2cell([repmat(genstr,8,1),['V1A';'V2A';'V3A';'V4A';'V1D';'V2D';'V3D';'V4D']],ones(8,1),(length(genstr)+3));
        dvLabels{9*(lamnum-1)+9} = [genstr,'t'];
    end
    
    % Bottom skin
    for j = 1:length(constant.lam.BotID{ii})
        lamnum = constant.lam.BotID{ii}(j);
        genstr = ['Bottom skin: Lam ',num2str(lamnum),', '];
        dvLabels(9*(lamnum-1)+(1:8),1) = mat2cell([repmat(genstr,8,1),['V1A';'V2A';'V3A';'V4A';'V1D';'V2D';'V3D';'V4D']],ones(8,1),(length(genstr)+3));
        dvLabels{9*(lamnum-1)+9} = [genstr,'t'];
    end
    
    % Spars
    if FLAGSPAR
        for j = 1:length(constant.lam.SparID{ii})
            lamnum = constant.lam.SparID{ii}(j);
            genstr = ['Spars: Lam ',num2str(lamnum),', '];
            dvLabels(9*(lamnum-1)+(1:8),1) = mat2cell([repmat(genstr,8,1),['V1A';'V2A';'V3A';'V4A';'V1D';'V2D';'V3D';'V4D']],ones(8,1),(length(genstr)+3));
            dvLabels{9*(lamnum-1)+9} = [genstr,'t'];
        end
    end
end

twistsec = find(constant.morph.inp.twist.sec);
for i=1:size(constant.morph.inp.twist.anglefin,1)
    dvLabelsTwist{i,1} = ['Twist morphing, Section ',num2str(twistsec(i))];
end

shearsec = find(constant.morph.inp.shear.sec);
for i=1:size(constant.morph.inp.shear.anglefin,1)
    dvLabelsShear{i,1} = ['Shear morphing, Section ',num2str(shearsec(i))];
end

foldsec = find(constant.morph.inp.fold.sec);
for i=1:size(constant.morph.inp.fold.anglefin,1)
    dvLabelsFold{i,1} = ['Fold morphing, Section ',num2str(foldsec(i))];
end

camberloc = find(constant.morph.inp.camber.loc);
for i=1:size(constant.morph.inp.camber.paramfin,1)
    dvLabelsCamber{i,1} = ['Camber morphing, Location ',num2str(camberloc(i))];
end

dvLabelsSpan= 'Span morphing';

dvLabelsMorph = [dvLabelsTwist;dvLabelsShear;dvLabelsFold;dvLabelsCamber;dvLabelsSpan];

for i=1:numlc
   dvLabels(9*Nlam+(i-1)*length(dvLabelsMorph)+(1:length(dvLabelsMorph))) = strcat(dvLabelsMorph,[', LC: ',num2str(i)]);
end

% Number of constraints
nconstr = constant.opt.constraint.LamFeasibility+...                        % Lamination parameters feasibility constraints
    +constant.opt.constraint.Blending+...                                   % Blending constraint, please verify!!!!
    +numlc*constant.opt.constraint.Eigval+...                               % Aeroelastic stability per loadcase
    +numlc*constant.opt.constraint.Trim+...                                 % Maximum and minimum angle of attack per loadcase
    +numlc*constant.opt.constraint.Aileff+...                               % Aileron effectiveness per loadcase
    +constant.opt.constraint.oneG+...                                       % 1g trim constraint
    +(numlc-numgust)*constant.opt.constraint.Strain+...                     % Strain constraints for non-gust loadcases
    +2*numgust*length(constant.gust.H)*constant.opt.constraint.Strain+...   % Strain constraints for gust loadcases (per gust length: positive, negative, min, max deflection)
    +(numlc-numgust)*constant.opt.constraint.StrainCrit+...                     % Strain constraints for non-gust loadcases
    +2*numgust*length(constant.gust.H)*constant.opt.constraint.StrainCrit+...   % Strain constraints for gust loadcases (per gust length: positive, negative, min, max deflection)
    +(numlc-numgust)*constant.opt.constraint.Buckl+...                      % Buckling constraints for non-gust loadcases
    +2*numgust*length(constant.gust.H)*constant.opt.constraint.Buckl;%+...    % Buckling constraints for gust loadcases (per gust length: positive, negative, min, max deflection)
%     +1;                                                                     % Morphing energy (number can vary depending on whether you split energy per loadcase, per morphing parameter)

%}
