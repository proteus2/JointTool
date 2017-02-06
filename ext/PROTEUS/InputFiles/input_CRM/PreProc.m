function [dvopt] = PreProc(constant,loadcase,dvinp,xval)

%%% Please update as required. Function should take a vector of design
%%% variables, xval, as input and provide the full input structure dvinp
%%% containing:

%%% dvopt.tail: 9x the number of laminates (4 lamination parameters, A, 4 lamination parameter, D, and thickness)
%%% In case of morphing, specify for each loadcase when required:
%%% dvopt.XX.XXstart: level of morphing the loadcase starts speeding up with
%%% dvopt.XX.XXini: level of morphing the loadcase reaches in speeding up
%%% dvopt.XX.XXfin: final level of morphing at end of morphing manoeuvre, in case morphing energy is not important, make ini and fin equal

xdenorm = constant.opt.DeNormFct(xval);
dvopt = dvinp; % Copy initial input to optimiser


%% Example only tailoring

% Update dvopt with optimiser outputs
dvopt.tail = xdenorm;
%}
%% Example tailoring and morphing
%{
Nlam = length(constant.lam.ID);
numlc = length(loadcase.EAS);

dvopt.tail = xdenorm(1:9*Nlam);

ntwist = size(constant.morph.inp.twist.anglefin,1); % Number of twist design variables per loadcase
nshear = size(constant.morph.inp.shear.anglefin,1); % Number of shear design variables per loadcase
nfold = size(constant.morph.inp.fold.anglefin,1); % Number of fold design variables per loadcase
ncamber = size(constant.morph.inp.camber.paramfin,1); % Number of camber design variables per loadcase
nspan = 1;

dvmorph = reshape(xdenorm(9*Nlam+1:end),[],numlc);
dvtwist = dvmorph(1:ntwist,:);
dvshear = dvmorph(ntwist+(1:nshear),:);
dvfold = dvmorph(ntwist+nshear+(1:nfold),:);
dvcamber = dvmorph(ntwist+nshear+nfold+(1:ncamber),:);
dvspan = dvmorph(ntwist+nshear+nfold+ncamber+(1:nspan),:);

for k=1:numlc
    % Twist:
    dvopt.twist.phifin(:,k) = constant.morph.twist.dangleddv*dvtwist(:,k); % Note design variables need to be converted to structural sections
    angleini = deg2rad(constant.morph.inp.twist.angleini(:,k));
    anglestart = zeros(sum(constant.morph.inp.twist.sec),1);
    for i=1:size(angleini,1)
        lcinp = constant.morph.inp.twist.anglelcinp(i,k); % Load previous loadcase if necessary
        if lcinp~=0
            angleini(i,1) = dvtwist(i,lcinp);
            anglestart(i,1) = dvtwist(i,lcinp);
        end
    end
    dvopt.twist.phistart(:,k) = constant.morph.twist.dangleddv*anglestart;
    dvopt.twist.phiini(:,k) = constant.morph.twist.dangleddv*angleini;
    
    % Shear:
    dvopt.shear.psifin(:,k) = constant.morph.shear.dangleddv*dvshear(:,k); % Note design variables need to be converted to structural sections
    angleini = deg2rad(constant.morph.inp.shear.angleini(:,k));
    anglestart = zeros(sum(constant.morph.inp.shear.sec),1);
    for i=1:size(angleini,1)
        lcinp = constant.morph.inp.shear.anglelcinp(i,k); % Load previous loadcase if necessary
        if lcinp~=0
            angleini(i,1) = dvshear(i,lcinp);
            anglestart(i,1) = dvshear(i,lcinp);
        end
    end
    dvopt.shear.psistart(:,k) = constant.morph.shear.dangleddv*anglestart;
    dvopt.shear.psiini(:,k) = constant.morph.shear.dangleddv*angleini;
    
    % Fold:
    dvopt.fold.thetafin(:,k) = constant.morph.fold.dangleddv*dvfold(:,k); % Note design variables need to be converted to structural sections
    angleini = deg2rad(constant.morph.inp.fold.angleini(:,k));
    anglestart = zeros(sum(constant.morph.inp.fold.sec),1);
    for i=1:size(angleini,1)
        lcinp = constant.morph.inp.fold.anglelcinp(i,k); % Load previous loadcase if necessary
        if lcinp~=0
            angleini(i,1) = dvfold(i,lcinp);
            anglestart(i,1) = dvfold(i,lcinp);
        end
    end
    dvopt.fold.thetastart(:,k) = constant.morph.fold.dangleddv*anglestart;
    dvopt.fold.thetaini(:,k) = constant.morph.fold.dangleddv*angleini;
    
    % Camber:
    dvopt.camber.paramfin(:,k) = constant.morph.camber.dparamddv*dvcamber(:,k); % Note design variables need to be converted to structural sections
    paramini = constant.morph.inp.camber.paramini(:,k);
    paramstart = zeros(length(constant.morph.camber.loc),1);
    paramind = zeros(length(constant.morph.camber.loc),1);
    for i=1:size(paramini,1)
        lcinp = constant.morph.inp.camber.paramlcinp(i,k); % Load previous loadcase if necessary
        if lcinp~=0
            paramini(i,1) = dvcamber(i,lcinp);
            paramupd = constant.morph.camber.dparamddv(:,i)*dvcamber(i,lcinp);
            paramstart(paramupd~=0,1) = paramstart(paramupd~=0,1)+paramupd(paramupd~=0);
            paramind(constant.morph.camber.dparamddv(:,i)~=0) = 1;
        end
    end
    paramstart(paramind==0) = constant.morph.camber.ini(paramind==0);
    dvopt.camber.paramstart(:,k) = paramstart;
    dvopt.camber.paramini(:,k) = constant.morph.camber.dparamddv*paramini;

    % Span:
    dvopt.span.extfin(:,k) = constant.morph.span.dextddv*dvspan(:,k); % Note design variables need to be converted to structural sections
    lcinp = constant.morph.inp.span.extlcinp(1,k); % Load previous loadcase if necessary
    if lcinp>0
        extini = dvspan(1,lcinp);
        extstart = dvspan(1,lcinp);
        dvopt.span.extstart(:,k) = constant.morph.span.dextddv*extstart;
        dvopt.span.extini(:,k) = constant.morph.span.dextddv*extini;
    elseif lcinp<0
        dvopt.span.extstart(:,k) = constant.morph.span.ext0;
        dvopt.span.extini(:,k) = constant.morph.span.ext0;
    else
        dvopt.span.extstart(:,k) = constant.morph.span.ext0;
        dvopt.span.extini(:,k) = constant.morph.span.ext0;
    end
end
%}
