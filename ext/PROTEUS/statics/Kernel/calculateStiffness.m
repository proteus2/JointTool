function [statics] = calculateStiffness(constant,statics,ders,tailflag)

% Extract parameters from structured arrays
frdof = constant.str.frdof;

fprintf('\n')
fprintf('--- Static Structure ---\n')
fprintf('Calculate stiffness.')

%=================================================================%
% STIFFNESS CALCULATION
%=================================================================%
statics.str.p = zeros(size(statics.str.p));
statics = structure(constant,statics,ders,tailflag,0,0);
statics = pgen_str(constant,statics,ders,0);
statics = fext (constant,statics,1,ders,0,tailflag,0); % Routine to create stiffness matrix from external forces
R       = statics.str.Fext; % Add aerodynamic forces
Rr      = R(frdof);
J       = statics.str.Ks - statics.str.Kfext;
Jr      = J(frdof,frdof);
plin    = Jr\Rr;

statics.str.p(frdof,1) = plin;

%%% Update of the structural, aerodynamic and external forces based on the
%%% structural deformations
statics.str.Fs   = statics.str.Ks*statics.str.p;
statics.str.Fext = statics.str.Fext + statics.str.Kfext*statics.str.p;

% Compute the corresponding internal forces
statics = structure_lin(constant,statics,ders);    

if ders==1
    for i = 1:numel(statics.str.C)
        dKsdC = reshape(statics.sens.pdC_Ks(:,i),size(statics.str.Ks,2),size(statics.str.Ks,1))';
        dpdC(:,i) = -(statics.str.Ks(frdof,frdof)-statics.str.Kfext(frdof,frdof))\dKsdC(frdof,frdof)*plin;
        dpfulldC = zeros(length(statics.str.p),1);
        dpfulldC(frdof) = dpdC(:,i);
        dFsdC(:,i) = statics.str.Ks*dpfulldC+dKsdC*statics.str.p;
    end

    statics.sens.dC_p  = zeros(length(statics.str.p),numel(statics.str.C));
    statics.sens.dC_p(frdof,:) = dpdC;
    statics.sens.dC_Ks = statics.sens.pdC_Ks;
    statics.sens.dC_Fs = dFsdC;
    statics.sens.dC_Fext = statics.str.Kfext*statics.sens.dC_p;
    statics.sens.dC_re = statics.sens.pdC_re+statics.sens.dp_re*statics.sens.dC_p;

    for i = 1:size(statics.sens.dKfextdt,2)
        dKextdt = reshape(statics.sens.dKfextdt(:,i),size(statics.str.Kfext,2),size(statics.str.Kfext,1))';
        dpdt(:,i) = (statics.str.Ks(frdof,frdof)-statics.str.Kfext(frdof,frdof))\(statics.sens.dFextdt(frdof,i)+dKextdt(frdof,frdof)*plin);

        dpfulldt = zeros(length(statics.str.p),1);
        dpfulldt(frdof) = dpdt(:,i);
        dFextdt(:,i) = statics.sens.dFextdt(:,i)+statics.str.Kfext*dpfulldt+dKextdt*statics.str.p;
    end
    statics.sens.dpdt  = zeros(length(statics.str.p),size(statics.sens.dKfextdt,2));
    statics.sens.dpdt(frdof,:) = dpdt;
    statics.sens.dFsdt = statics.str.Ks*statics.sens.dpdt;
    statics.sens.dredt = statics.sens.dp_re*statics.sens.dpdt;

    statics.sens.dKfextdt = statics.sens.dKfextdt;
    statics.sens.dFextdt = dFextdt;
end
%=========================================================================%
% END
%=========================================================================%
