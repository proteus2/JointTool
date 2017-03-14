function [dynamics] = genDynStrModel(constant,lampar,crossmod,statics,lin,varargin)

if ~isempty(varargin)
   echoflag = varargin{1};
else
   echoflag = 1; 
end

% Set flags
feas = 1;
sens = 0;
tailflag = 1;

% Store current directory
curdir = cd;

% Prepare sensitivities
if sens == 1
    if tailflag == 1
        if lin == 1
            % Linear
            cd('ext/PROTEUS/dynamics/structure')
            [presens.dKredddv,presens.dKrsddv] = selK((statics.sens.dC_Ks)*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv)-statics.sens.dKfextdt*lampar.dtddv,constant);
            cd(curdir)
            presens.dKlddv = statics3.sens.dC_Kl*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv);
        else
            % Non-linear
            cd('ext/PROTEUS/dynamics/structure')
            [presens.dKredddv,presens.dKrsddv] = selK(((statics.sens.dC_Ks-statics.sens.dC_Kfext)*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv)+(statics.sens.dKsdt-statics.sens.dKfextdt)*lampar.dtddv),constant);
            cd(curdir)
            presens.dKlddv = statics.sens.dC_Kl*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv);
        end
        presens.dKredddv = sparse(presens.dKredddv);
        presens.dKlddv = sparse(presens.dKlddv);
        presens.dKrsddv = sparse(presens.dKrsddv);
    end
else
    presens = [];
end

str.cam = interp1(constant.inp.xyz(2:3:end),constant.inp.Aerofoil.NodeCamber',constant.str.xyz(2:3:end))';

% Required inputs
if feas == 0 || lin == 1
    str.xyz = constant.str.xyz;
    str.xyz0 = constant.str.xyz;
    str.statp = zeros(length(statics.str.p),1);
else
    str.xyz = statics.str.xdef;
    str.statp = statics.str.p;
    str.xyz0 = constant.str.xyz;
end

xref_ini = interp1(constant.inp.xyz(2:3:end),constant.inp.xref,constant.str.xyz(2:3:end));

c_ini = interp1(constant.inp.xyz(2:3:end),constant.inp.c,constant.str.xyz(2:3:end));
str.c = c_ini;
str.xref = xref_ini.*str.c;
str.dxrefdc = diag(xref_ini);

str.thetaini = interp1(constant.inp.xyz(2:3:end),constant.inp.theta,constant.str.xyz(2:3:end));

str.K    = statics.str.Ks-statics.str.Kfext(constant.str.dof.wing.all,constant.str.dof.wing.all);
str.Kl   = statics.str.Kl;

if sens == 1
    if tailflag == 1
        str.dKddv = presens.dKredddv;
        str.dKrsddv = presens.dKrsddv;
        str.dKlddv = presens.dKlddv;
        if lin == 1
            str.dstatpddv = 0*(statics.sens.dC_p*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv)+statics.sens.dpdt*lampar.dtddv);
        else
            str.dstatpddv = statics.sens.dC_p*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv)+statics.sens.dpdt*lampar.dtddv;
        end
    end
end

str.elm.C = crossmod.C;
str.elm.A = crossmod.mA;
str.elm.I = crossmod.mI;
str.elm.Q = crossmod.mQ;

if feas == 0 || lin == 1
    str.e1elm = constant.str.e1;
    str.e2elm = constant.str.e2;
    str.e3elm = constant.str.e3;
else
    str.e1elm = statics.str.e1;
    str.e2elm = statics.str.e2;
    str.e3elm = statics.str.e3;
end

str.rho = 1;
str.eft = constant.str.EFT;

% DOF 
str.dof = constant.str.dof;

% Structures
str.Nel   = size(str.xyz,1)/3-1; % Number of structural elements
str.Nnode = str.Nel+1;

%% Structural State-Space
cd('ext/PROTEUS/dynamics/structure')
[str] = strorien(str,sens,lin);

if echoflag
    fprintf('\n')
    fprintf('--- Dynamic Structure ---\n')
    fprintf('Calculate mass matrix\n')
end
[str] = dstmam(str,constant,sens,tailflag);

if sens == 1
    if tailflag == 1
        if lin == 1
            dMdt      = str.dMdA*crossmod.dmAdt+str.dMdQ*crossmod.dmQdt+str.dMdI*crossmod.dmIdt;
            str.dMddv = ((str.dMdC)*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv)+dMdt*lampar.dtddv);
            
            dMrsdt      = str.dMrsdA*crossmod.dmAdt+str.dMrsdQ*crossmod.dmQdt+str.dMrsdI*crossmod.dmIdt;
            str.dMrsddv = ((str.dMrsdC)*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv)+dMrsdt*lampar.dtddv);
        else
            dMdp = str.dMdR*(str.dRde1*statics.sens.dp_e1+str.dRde2*statics.sens.dp_e2+str.dRde3*statics.sens.dp_e3);
            dMdt = dMdp*statics.sens.dpdt+str.dMdA*crossmod.dmAdt+str.dMdQ*crossmod.dmQdt+str.dMdI*crossmod.dmIdt;
            str.dMddv = ((dMdp*statics.sens.dC_p+str.dMdC)*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv)+dMdt*lampar.dtddv);
            
            dMrsdp = str.dMrsdR*(str.dRde1*statics.sens.dp_e1+str.dRde2*statics.sens.dp_e2+str.dRde3*statics.sens.dp_e3);
            dMrsdt = dMrsdp*statics.sens.dpdt+str.dMrsdA*crossmod.dmAdt+str.dMrsdQ*crossmod.dmQdt+str.dMrsdI*crossmod.dmIdt;
            str.dMrsddv = ((dMrsdp*statics.sens.dC_p+str.dMrsdC)*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv)+dMrsdt*lampar.dtddv);
        end
    end
end

if echoflag
    fprintf('Assemble structural state space\n')
end
[str] = dstrss_parfor(str,sens,tailflag);
cd(curdir)

%% Output
dynamics = str;

end
