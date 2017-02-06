% function varargout = analysis_beam(dv,ders,trim,outflag,constant,varargin)
function [crossmod,statics3,dynamics] = analysis_beam(constant,lin,trim,grav,varargin)

% Derivative information not required, since no optimisation

% if length(varargin) == 1
%     normal = varargin{1};
%     normflag = 1;
% else
%     normflag = 0;
% end

curdir = cd;

%% Cross-sectional Modeller
tic
cd('cross_mod')
crossmod = cross_mod_beam(constant);
cd(curdir)
toc

%% Include gravity
if grav == 1
    fext = constant.fext;
    nfext = length(fext.type);
    fext.type{nfext+1} = 'Structural mass';
    count = 1;
    for i = 1:constant.str.Ns
        Mmag(i) = crossmod.mA(i)*norm(constant.str.xyz(3*i+(1:3))-constant.str.xyz(3*(i-1)+(1:3)));
        Fmag = Mmag(i)*9.81;
        fext.magnitude{nfext+1}(count,:) = [0,0,-Fmag,0,0,0];
        fext.location{nfext+1}(count,:) = [(constant.str.xyz(3*i+(1:3))+constant.str.xyz(3*(i-1)+(1:3)))/2]';
        fext.location{nfext+1}(count,:) = fext.location{2}(count,:) + (constant.str.R0(3*(i-1)+(1:3),:)*[0;crossmod.mQ(i,2)/crossmod.mA(i);crossmod.mQ(i,1)/crossmod.mA(i)])'; % Think about positive definition of cg(1) and cg(2)
        fext.follower{nfext+1}(1,count)=0;
        count = count+1;
    end
    cd('input')
    constant.fext=fext_inp(constant.str,fext);
    cd(curdir)
end
%% Statics
tic
cd('statics')
cd('Kernel')

statics.str.C = crossmod.C;
statics.str.p = zeros(constant.str.Ndof,1);

[statics2,exitflag2] = statae(constant,statics,0,constant.aero.V,constant.aero.alpha0,constant.aero.rho,0,0,lin);

if trim==1
    [statics2,exitflag2] = statae(constant,statics2,constant.aero.V-1,constant.aero.V,constant.aero.alpha0,constant.aero.rho,trim,0,lin);
else
    statics2.str.alpha = constant.aero.alpha0;
end

cd('../Postprocessor')
statics3 = mr(constant,statics2,0,trim);

cd(curdir)
toc


%% Dynamics
tic
cd('dynamics')
[dynamics] = dynamic_beam(constant,[],crossmod,statics3,0,1,lin,[]);
cd(curdir)
toc

%{
%% Strain computation

tic
cd('strain')
[strain] = strain_comp(constant,lampar,crossmod,statics3,ders);
cd(curdir)
toc

%% Objective and constraints
% Select stiffness matrix sensitivities
tic
if ders == 1
    cd('postproc')
    dC_Ks = selK(statics3.sens.dC_Ks,constant.str.Ns);
    cd(curdir)
end
toc
% keyboard
if normflag == 1
    if feas==1
        obj = norm(squeeze(out.y(out.crit_lc,end-2:end,out.tmaxind(out.crit_lc)))'+statics3.str.Mr)./abs(normal.obj);
    else
        obj = norm(statics3.str.Mr)./abs(normal.obj);
    end
    
    g = -[lampar.g1;lampar.g2;lampar.g3;lampar.g4;out.eigval;constant.lam.exmax-strain.ex;constant.lam.gmax-strain.gamma;constant.general.alpha_max-statics3.str.alpha;constant.general.alpha_max+statics3.str.alpha]./abs(normal.g);
    if ders == 1
        dg = [lampar.dg1;lampar.dg2;lampar.dg3;lampar.dg4];
        
        dl = (out.dldxyz*statics3.sens.dp_xdef*statics3.sens.dC_p+...
            out.dldc0*statics3.sens.dp_c0*statics3.sens.dC_p+...
            out.dldC+...
            out.dldK*dC_Ks+...
            out.dlde1*statics3.sens.dp_e1*statics3.sens.dC_p+...
            out.dlde2*statics3.sens.dp_e2*statics3.sens.dC_p+...
            out.dlde3*statics3.sens.dp_e3*statics3.sens.dC_p+...
            out.dlda*statics3.sens.dC_a)*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv);
        
        dg = -[dg;dl;-strain.dexddv;-strain.dgammaddv;-statics3.sens.dC_a*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv);statics3.sens.dC_a*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv)]./abs(normal.g(:,ones(1,length(dv))));
        if feas==1
            dobj = (out.dydxyz(end-2:end,:)*statics3.sens.dp_xdef*statics3.sens.dC_p+...
                out.dydc0(end-2:end,:)*statics3.sens.dp_c0*statics3.sens.dC_p+...
                out.dydC(end-2:end,:)+...
                out.dydK(end-2:end,:)*dC_Ks+...
                out.dyde1(end-2:end,:)*statics3.sens.dp_e1*statics3.sens.dC_p+...
                out.dyde2(end-2:end,:)*statics3.sens.dp_e2*statics3.sens.dC_p+...
                out.dyde3(end-2:end,:)*statics3.sens.dp_e3*statics3.sens.dC_p+...
                out.dyda(end-2:end,:)*statics3.sens.dC_a+statics3.sens.dC_Mr)*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv);
            dobj = dnorm(squeeze(out.y(out.crit_lc,end-2:end,out.tmaxind(out.crit_lc)))'+statics3.str.Mr,dobj)./abs(normal.obj(:,ones(1,length(dv))));
        else
            dobj = (statics3.sens.dC_Mr)*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv);
            dobj = dnorm(statics3.str.Mr,dobj)./abs(normal.obj(:,ones(1,length(dv))));
        end
    end
else
    
    if feas==1
%         keyboard
        obj = norm(squeeze(out.y(out.crit_lc,end-2:end,out.tmaxind(out.crit_lc)))'+statics3.str.Mr);
    else
        obj = norm(statics3.str.Mr);
    end
    g = -[lampar.g1;lampar.g2;lampar.g3;lampar.g4;out.eigval;constant.lam.exmax-strain.ex;constant.lam.gmax-strain.gamma;constant.general.alpha_max-statics3.str.alpha;constant.general.alpha_max+statics3.str.alpha];
    if ders == 1
        dg = [lampar.dg1;lampar.dg2;lampar.dg3;lampar.dg4;-statics3.sens.dC_a*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv);statics3.sens.dC_a*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv)];
        
        dl = (out.dldxyz*statics3.sens.dp_xdef*statics3.sens.dC_p+...
            out.dldc0*statics3.sens.dp_c0*statics3.sens.dC_p+...
            out.dldC+...
            out.dldK*dC_Ks+...
            out.dlde1*statics3.sens.dp_e1*statics3.sens.dC_p+...
            out.dlde2*statics3.sens.dp_e2*statics3.sens.dC_p+...
            out.dlde3*statics3.sens.dp_e3*statics3.sens.dC_p+...
            out.dlda*statics3.sens.dC_a)*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv);
        
        dg = -[dg;dl;-strain.dexddv;-strain.dgammaddv;-statics3.sens.dC_a*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv);statics3.sens.dC_a*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv)];
        if feas==1
            dobj = (out.dydxyz(end-2:end,:)*statics3.sens.dp_xdef*statics3.sens.dC_p+...
                out.dydc0(end-2:end,:)*statics3.sens.dp_c0*statics3.sens.dC_p+...
                out.dydC(end-2:end,:)+...
                out.dydK(end-2:end,:)*dC_Ks+...
                out.dyde1(end-2:end,:)*statics3.sens.dp_e1*statics3.sens.dC_p+...
                out.dyde2(end-2:end,:)*statics3.sens.dp_e2*statics3.sens.dC_p+...
                out.dyde3(end-2:end,:)*statics3.sens.dp_e3*statics3.sens.dC_p+...
                out.dyda(end-2:end,:)*statics3.sens.dC_a+statics3.sens.dC_Mr)*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv);
            dobj = dnorm(squeeze(out.y(out.crit_lc,end-2:end,out.tmaxind(out.crit_lc)))'+statics3.str.Mr,dobj);
        else
            dobj = (statics3.sens.dC_Mr)*(crossmod.dCdA*lampar.dAddv+crossmod.dCdD*lampar.dDddv);
            dobj = dnorm(statics3.str.Mr,dobj);
        end
    end
end
obj
max(g)

if outflag == 1
    varargout{1} = statics3;
    dynamics.str = str;
    dynamics.aero = aero;
    dynamics.coupl = coupl;
    dynamics.out = out;
    varargout{2} = dynamics;
    varargout{3} = strain;
elseif normflag == 0
    normal.obj = obj;
    normal.g = g;
    varargout{1} = normal;
else
    varargout{1} = obj;
    varargout{2} = g;
    
    if ders==1
        varargout{3} = dobj';
        varargout{4} = dg;
    end
end


function [dnorma] = dnorm(a,da)

dnorma = 1/2/norm(a)*(sum(2*a(:,ones(1,size(da,2))).*da));


%% % Wing plot
% cd('postproc')
% wingplot(reshape(xlz',[],1),c,reshape(c0',[],1),xref.*c,aero.xlz)
% cd(curdir)

%}