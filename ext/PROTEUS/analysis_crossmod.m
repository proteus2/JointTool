function [varargout] = analysis_crossmod(constant,crossmod,dispi,alphai,ANALYSISTYPE,trim,trimdef,grav,ders,outflag)

curdir = cd;

if ANALYSISTYPE == 1 || ANALYSISTYPE == 3
    tailflag = 1;
else
    tailflag = 0;
end
if ANALYSISTYPE == 2 || ANALYSISTYPE == 3
    morphflag = 1;
else
    morphflag = 0;
end

% Structural mass
if ders == 1 && tailflag == 1
    dMmagdt = zeros(constant.str.Ns,length(constant.lam.ID));
end
for i = 1:constant.str.Ns
    Mmag(i) = crossmod.mA(i)*norm(constant.str.xyz(3*i+(1:3))-constant.str.xyz(3*(i-1)+(1:3)));
    if ders == 1 && tailflag == 1
        dMmagdt(i,:) = crossmod.dmAdt(i,:)*norm(constant.str.xyz(3*i+(1:3))-constant.str.xyz(3*(i-1)+(1:3)));
    end
end

mass = sum(Mmag);
constant.general.mass = mass;  % Store Mass
if ders == 1 && tailflag == 1
    dmassdt = sum(dMmagdt);
end

% Include gravity
if grav == 1
    fext = constant.fext;
    
    if isfield(fext,'type')
        nfext = length(fext.type);
    else
        nfext = 0;
    end
    
    fext.type{nfext+1} = 'Structural mass'; % Note, don't change this name, this is used to update gravitational forces in case of span extension
    
    count = 1;
    if ders == 1 && tailflag == 1
        fext.dmagnitudedt{nfext+1} = zeros(6*constant.str.Ns,length(constant.lam.ID));
        fext.dlocationdt{nfext+1} = zeros(3*constant.str.Ns,length(constant.lam.ID));
    end
    for i = 1:constant.str.Ns
        Fmag = Mmag(i)*9.81;
        fext.magnitude{nfext+1}(count,:) = constant.general.nz*[0,0,-Fmag,0,0,0];
        fext.location{nfext+1}(count,:)  = [(constant.str.xyz(3*i+(1:3))+constant.str.xyz(3*(i-1)+(1:3)))/2]';
        fext.location{nfext+1}(count,:)  = fext.location{nfext+1}(count,:) + (constant.str.R0(3*(i-1)+(1:3),:)*[0;crossmod.mQ(i,2)/crossmod.mA(i);crossmod.mQ(i,1)/crossmod.mA(i)])'; % Think about positive definition of cg(1) and cg(2)
        fext.follower{nfext+1}(1,count)  = 0;
        fext.alphaflag{nfext+1}(1,count)  = 1;
        if ders == 1 && tailflag == 1
            fext.dmagnitudedt{nfext+1}(6*(count-1)+3,:) = -constant.general.nz*9.81*dMmagdt(i,:);
            fext.dlocationdt{nfext+1}(3*(count-1)+(1:3),:) = constant.str.R0(3*(i-1)+(1:3),:)*[zeros(1,length(constant.lam.ID));(crossmod.mA(i)*crossmod.dmQdt(2*(i-1)+2,:)-crossmod.mQ(i,2)*crossmod.dmAdt(i,:))/crossmod.mA(i)^2;(crossmod.mA(i)*crossmod.dmQdt(2*(i-1)+1,:)-crossmod.mQ(i,1)*crossmod.dmAdt(i,:))/crossmod.mA(i)^2];
        end
        count = count+1;
    end
    cd([constant.curdir,'/ext/PROTEUS/InputFiles'])
    if ders == 1 && tailflag == 1
        constant.fext=dfext_inp(constant.str,fext,ders);
    else
        constant.fext=dfext_inp(constant.str,fext,0);
    end
    cd(constant.curdir)
end

% Calculate Non-structural masses
nonstrmass = 0;
for i=1:length(constant.lumped.type)
    nonstrmass = nonstrmass+sum(constant.lumped.mass{i});
end

constant.general.nonstrmass = nonstrmass; % store total non-structural mass

% Define trim mass
if trimdef == 1
    statics.str.weight = constant.general.nz*(constant.general.weight + 9.81*mass + 9.81*nonstrmass);
    
    if ders == 1 && tailflag == 1
        statics.str.dweightdt = constant.general.nz*9.81*dmassdt;
    end
elseif trimdef == 2
    statics.str.weight = constant.general.nz*constant.general.weight;
    if ders == 1 && tailflag == 1
        statics.str.dweightdt = zeros(1,length(constant.lam.ID));
    end
end

% Structural Analysis
cd('ext/PROTEUS/statics')
cd('Kernel')

statics.str.C= crossmod.C;

statics.str.p = dispi; 

if grav == 1
    statics.cross.mA = crossmod.mA;
    if ders == 1 && tailflag == 1
        statics.cross.dmAdt = crossmod.dmAdt;
    end
end

if trim == 1
    statics.str.alpha = alphai;
else
    statics.str.alpha = constant.aero.alpha0;
end

statics2 = calculateStiffness(constant,statics,ders,tailflag);

cd('../Postprocessor')
statics3 = mr(constant,statics2,ders,tailflag,morphflag);
cd(curdir)

% Clear memory
clear statics2
if grav == 1
    statics3 = rmfield(statics3,'cross');
end

% Quick fix
% Buckling and strain computation should be implemented as external
% function in the class strModel. (ADJUST A.S.A.P.)
% if constant.general.romflag ~= 2
%     
%     if constant.cross.numel ~= 0
%         % Strain computation
%         cd('ext/PROTEUS/strain')
%         [straindat,strainraw] = strain_comp_single_gust_fail_index(constant,crossmod,statics3,[],ders,gustflag,tailflag);
%         strain.exmax = straindat.exmax;
%         strain.exmin = straindat.exmin;
%         strain.gammamax = straindat.gammamax;
%         strain.rcrit = straindat.rcrit;
%         strain.exmax_full{1} = straindat.exmax_full;
%         strain.exmin_full{1} = straindat.exmin_full;
%         strain.gamma_full{1} = straindat.gamma_full;
%         strain.r_full{1} = straindat.r_full;
%         strain.locemax{1} = straindat.locemax;
%         strain.locemin{1} = straindat.locemin;
%         strain.locgmax{1} = straindat.locgmax;
%         strain.locrcrit{1} = straindat.locrcrit;
%         strain.node1{1} = straindat.node1;
%         strain.node2{1} = straindat.node2;
%         fprintf('\n Strain ;')
%         cd(curdir)
%         
%         % Buckling computation
%         if constant.opt.BucklConst
%             if feas == 0
%                 buckl.r = 2*ones(8*constant.opt.constraint.BucklPerLam*length(constant.lam.ID),1);
%             else
%                 cd('buckling')
%                 constant2.buckl = constant.buckl;
%                 constant2.str = constant.str;
%                 constant2.lam = constant.lam;
%                 constant2.opt.constraint.BucklPerLam = constant.opt.constraint.BucklPerLam;
%                 Numcross = crossmod.Numcross;
%                 parfor i = 1
%                     [buckldat{i}] = buckl_comp_single_gust(constant2,lampar,Numcross,strainraw,ders,gustflag,tailflag);
%                 end
%                 buckl.r = buckldat{1}.r;
%                 buckl.rfull{1} = buckldat{1}.rfull;
%                 buckl.pan{1} = buckldat{1}.pan;
%                 buckl.sec{1} = buckldat{1}.sec;
%                 buckl.cross{1} = buckldat{1}.cross;
%                 fprintf('\n Buckling ;')
%                 cd(curdir)
%             end
%         end
%         
%         % Objective and constraints
%         if feas == 0
%             g.exmax = repmat(strain.exmax,ngust,1);
%             g.exmin = repmat(strain.exmin,ngust,1);
%             g.gammamax = repmat(strain.gammamax,ngust,1);
%             g.strainrcrit = repmat(strain.rcrit,ngust,1);
%             if constant.opt.BucklConst
%                 g.buckl = repmat(buckl.r,ngust,1);
%             end
%         else
%             g.exmax = strain.exmax;
%             g.exmin = strain.exmin;
%             g.gammamax = strain.gammamax;
%             g.strainrcrit = strain.rcrit;
%             if constant.opt.BucklConst
%                 g.buckl = buckl.r;
%             end
%         end
%         
%         obj.mass  = mass;
%         
%     else
%         strain = [];
%         buckl = [];
%     end
% end

% if constant.opt.BucklConst
% 	constant.buckl = rmfield(constant.buckl,'stiff');
% end

if outflag == 1
    varargout{1} = constant;
    varargout{2} = statics3;
%     if constant.general.romflag ~= 2
%         varargout{3} = strain;
%         if constant.opt.BucklConst
%             varargout{4} = buckl;
%         end
%     end
else
    varargout{1} = statics3.str.p;
    varargout{2} = statics3.str.alpha;
    varargout{3} = obj;
    varargout{4} = g;
    
    if ders==1
        varargout{5} = dobj;
        varargout{6} = dg;
    end
end

end

