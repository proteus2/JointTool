function [strain,strainraw] = strain_comp(constant,crossmod,statics,floc,ders,gustflag,tailflag,varargin)

if length(varargin) == 1
    morph = varargin{1};
    morphflag = 1;
else
    morphflag = 0;
end

for i=1:constant.str.Ns
    for ncross=1:length(crossmod.Ccross{i})
        if morphflag == 1
            [constr{crossmod.Numcross{i}(ncross),1},exmaxfull{crossmod.Numcross{i}(ncross),1},exminfull{crossmod.Numcross{i}(ncross),1},gammafull{crossmod.Numcross{i}(ncross),1},rfull{crossmod.Numcross{i}(ncross),1},strainraw{crossmod.Numcross{i}(ncross),1}] = ...
                straincomp(constant,crossmod,statics,floc,i,ncross,ders,gustflag,tailflag,morph);
        else
            [constr{crossmod.Numcross{i}(ncross),1},exmaxfull{crossmod.Numcross{i}(ncross),1},exminfull{crossmod.Numcross{i}(ncross),1},gammafull{crossmod.Numcross{i}(ncross),1},rfull{crossmod.Numcross{i}(ncross),1},strainraw{crossmod.Numcross{i}(ncross),1}] = ...
                straincomp(constant,crossmod,statics,floc,i,ncross,ders,gustflag,tailflag);
        end
    end
end

strain.locemax = [];
strain.locemin = [];
strain.locgmax = [];
strain.locrcrit = [];
strain.exmax = [];
strain.exmin = [];
strain.gammamax = [];
strain.rcrit = [];
if ders == 1
    if tailflag == 1
        strain.dexmaxdA = [];
        strain.dexmaxdD = [];
        strain.dexmaxdt = [];
        
        strain.dexmindA = [];
        strain.dexmindD = [];
        strain.dexmindt = [];
        
        strain.dgammamaxdA = [];
        strain.dgammamaxdD = [];
        strain.dgammamaxdt = [];
        
        strain.drcritdA = [];
        strain.drcritdD = [];
        strain.drcritdt = [];
        
        
        if gustflag == 1
            strain.dexmaxddv = [];
            strain.dexminddv = [];
            strain.dgammamaxddv = [];
            strain.drcritddv = [];
        end
    end
    if morphflag == 1
        if morph.camber
            strain.dexmaxdparam = [];
            strain.dexmindparam = [];
            strain.dgammamaxdparam = [];
            strain.drcritdparam = [];
        end
        if morph.twist
            strain.dexmaxdphi = [];
            strain.dexmindphi = [];
            strain.dgammamaxdphi = [];
            strain.drcritdphi = [];
        end
        if morph.shear
            strain.dexmaxdpsi = [];
            strain.dexmindpsi = [];
            strain.dgammamaxdpsi = [];
            strain.drcritdpsi = [];
        end
        if morph.span
            strain.dexmaxdext = [];
            strain.dexmindext = [];
            strain.dgammamaxdext = [];
            strain.drcritdext = [];
        end
        if morph.fold
            strain.dexmaxdtheta = [];
            strain.dexmindtheta = [];
            strain.dgammamaxdtheta = [];
            strain.drcritdtheta = [];
        end
    end
end


for i=1:constant.str.Ns
    for ncross=1:length(crossmod.Ccross{i})
        Numcross = crossmod.Numcross{i}(ncross);
        % Find the actual physical location of the strain data
        strain.locemax(end+(1:size(constr{Numcross,1}.locemax,1)),:) = (constant.cross.yzlocal{i}{ncross}(crossmod.Elm2D{i}{ncross}(constr{Numcross,1}.locemax(:,1),3),:)+constant.cross.yzlocal{i}{ncross}(crossmod.Elm2D{i}{ncross}(constr{Numcross,1}.locemax(:,1),4),:))/2;
        strain.locgmax(end+(1:size(constr{Numcross,1}.locgmax,1)),:) = (constant.cross.yzlocal{i}{ncross}(crossmod.Elm2D{i}{ncross}(constr{Numcross,1}.locgmax(:,1),3),:)+constant.cross.yzlocal{i}{ncross}(crossmod.Elm2D{i}{ncross}(constr{Numcross,1}.locgmax(:,1),4),:))/2;
        strain.locrcrit(end+(1:size(constr{Numcross,1}.locrcrit,1)),:) = (constant.cross.yzlocal{i}{ncross}(crossmod.Elm2D{i}{ncross}(constr{Numcross,1}.locrcrit(:,1),3),:)+constant.cross.yzlocal{i}{ncross}(crossmod.Elm2D{i}{ncross}(constr{Numcross,1}.locrcrit(:,1),4),:))/2;
        strain.locemin(end+(1:size(constr{Numcross,1}.locemin,1)),:) = (constant.cross.yzlocal{i}{ncross}(crossmod.Elm2D{i}{ncross}(constr{Numcross,1}.locemin(:,1),3),:)+constant.cross.yzlocal{i}{ncross}(crossmod.Elm2D{i}{ncross}(constr{Numcross,1}.locemin(:,1),4),:))/2;
        
        strain.exmax(end+(1:size(constr{Numcross,1}.locemax,1)),1) = constr{Numcross,1}.exmax(:,1);
        strain.gammamax(end+(1:size(constr{Numcross,1}.locgmax,1)),1) = constr{Numcross,1}.gammamax(:,1);
        strain.rcrit(end+(1:size(constr{Numcross,1}.locrcrit,1)),1) = constr{Numcross,1}.rcrit(:,1);
        strain.exmin(end+(1:size(constr{Numcross,1}.locemin,1)),1) = constr{Numcross,1}.exmin(:,1);
        
        if ders == 1
            if tailflag == 1
                strain.dexmaxdA(end+(1:size(constr{Numcross,1}.locemax,1)),:) = constr{Numcross,1}.dexmaxdA;
                strain.dexmaxdD(end+(1:size(constr{Numcross,1}.locemax,1)),:) = constr{Numcross,1}.dexmaxdD;
                strain.dexmaxdt(end+(1:size(constr{Numcross,1}.locemax,1)),:) = constr{Numcross,1}.dexmaxdt;
                
                strain.dgammamaxdA(end+(1:size(constr{Numcross,1}.locgmax,1)),:) = constr{Numcross,1}.dgammamaxdA;
                strain.dgammamaxdD(end+(1:size(constr{Numcross,1}.locgmax,1)),:) = constr{Numcross,1}.dgammamaxdD;
                strain.dgammamaxdt(end+(1:size(constr{Numcross,1}.locgmax,1)),:) = constr{Numcross,1}.dgammamaxdt;
                
                strain.drcritdA(end+(1:size(constr{Numcross,1}.locrcrit,1)),:) = constr{Numcross,1}.drcritdA;
                strain.drcritdD(end+(1:size(constr{Numcross,1}.locrcrit,1)),:) = constr{Numcross,1}.drcritdD;
                strain.drcritdt(end+(1:size(constr{Numcross,1}.locrcrit,1)),:) = constr{Numcross,1}.drcritdt;
                
                
                strain.dexmindA(end+(1:size(constr{Numcross,1}.locemin,1)),:) = constr{Numcross,1}.dexmindA;
                strain.dexmindD(end+(1:size(constr{Numcross,1}.locemin,1)),:) = constr{Numcross,1}.dexmindD;
                strain.dexmindt(end+(1:size(constr{Numcross,1}.locemin,1)),:) = constr{Numcross,1}.dexmindt;
                
                if gustflag == 1
                    strain.dexmaxddv(end+(1:size(constr{Numcross,1}.locemax,1)),:) = constr{Numcross,1}.dexmaxddv;
                    strain.dgammamaxddv(end+(1:size(constr{Numcross,1}.locgmax,1)),:) = constr{Numcross,1}.dgammamaxddv;
                    strain.drcritddv(end+(1:size(constr{Numcross,1}.locrcrit,1)),:) = constr{Numcross,1}.drcritddv;
                    strain.dexminddv(end+(1:size(constr{Numcross,1}.locemin,1)),:) = constr{Numcross,1}.dexminddv;
                end
            end
            if morphflag == 1
                if morph.camber
                    strain.dexmaxdparam(end+(1:size(constr{Numcross,1}.locemax,1)),:) = constr{Numcross,1}.dexmaxdparam;
                    strain.dgammamaxdparam(end+(1:size(constr{Numcross,1}.locgmax,1)),:) = constr{Numcross,1}.dgammamaxdparam;
                    strain.dexmindparam(end+(1:size(constr{Numcross,1}.locemin,1)),:) = constr{Numcross,1}.dexmindparam;
                    strain.drcritdparam(end+(1:size(constr{Numcross,1}.locrcrit,1)),:) = constr{Numcross,1}.drcritdparam;
                end
                if morph.twist
                    strain.dexmaxdphi(end+(1:size(constr{Numcross,1}.locemax,1)),:) = constr{Numcross,1}.dexmaxdphi;
                    strain.dgammamaxdphi(end+(1:size(constr{Numcross,1}.locgmax,1)),:) = constr{Numcross,1}.dgammamaxdphi;
                    strain.dexmindphi(end+(1:size(constr{Numcross,1}.locemin,1)),:) = constr{Numcross,1}.dexmindphi;
                    strain.drcritdphi(end+(1:size(constr{Numcross,1}.locrcrit,1)),:) = constr{Numcross,1}.drcritdphi;
                end
                if morph.shear
                    strain.dexmaxdpsi(end+(1:size(constr{Numcross,1}.locemax,1)),:) = constr{Numcross,1}.dexmaxdpsi;
                    strain.dgammamaxdpsi(end+(1:size(constr{Numcross,1}.locgmax,1)),:) = constr{Numcross,1}.dgammamaxdpsi;
                    strain.dexmindpsi(end+(1:size(constr{Numcross,1}.locemin,1)),:) = constr{Numcross,1}.dexmindpsi;
                    strain.drcritdpsi(end+(1:size(constr{Numcross,1}.locrcrit,1)),:) = constr{Numcross,1}.drcritdpsi;
                end
                if morph.span
                    strain.dexmaxdext(end+(1:size(constr{Numcross,1}.locemax,1)),:) = constr{Numcross,1}.dexmaxdext;
                    strain.dgammamaxdext(end+(1:size(constr{Numcross,1}.locgmax,1)),:) = constr{Numcross,1}.dgammamaxdext;
                    strain.dexmindext(end+(1:size(constr{Numcross,1}.locemin,1)),:) = constr{Numcross,1}.dexmindext;
                    strain.drcritdext(end+(1:size(constr{Numcross,1}.locrcrit,1)),:) = constr{Numcross,1}.drcritdext;
                end
                if morph.fold
                    strain.dexmaxdtheta(end+(1:size(constr{Numcross,1}.locemax,1)),:) = constr{Numcross,1}.dexmaxdtheta;
                    strain.dgammamaxdtheta(end+(1:size(constr{Numcross,1}.locgmax,1)),:) = constr{Numcross,1}.dgammamaxdtheta;
                    strain.dexmindtheta(end+(1:size(constr{Numcross,1}.locemin,1)),:) = constr{Numcross,1}.dexmindtheta;
                    strain.drcritdtheta(end+(1:size(constr{Numcross,1}.locrcrit,1)),:) = constr{Numcross,1}.drcritdtheta;
                end
            end
        end
        
        strain.node1{i}{ncross} = constant.cross.yzlocal{i}{ncross}(crossmod.Elm2D{i}{ncross}(:,3),:);
        strain.node2{i}{ncross} = constant.cross.yzlocal{i}{ncross}(crossmod.Elm2D{i}{ncross}(:,4),:);
    end
end

strain.exmax_full = exmaxfull;
strain.exmin_full = exminfull;
strain.gamma_full = gammafull;
strain.r_full = rfull;

function [constr,exmaxmat,exminmat,gammamat,rmat,strainraw] = straincomp(constant,crossmod,statics,floc,i,ncross,ders,gustflag,tailflag,varargin)

if length(varargin) == 1
    morph = varargin{1};
    morphflag = 1;
else
    morphflag = 0;
end

nlamskin = length(crossmod.laminates{i}{ncross}(crossmod.laminate_type{i}{ncross}==1));
exmax = zeros(4*nlamskin,1);
exmin = zeros(4*nlamskin,1);
gammamax = zeros(8*nlamskin,1);
rcrit = zeros(8*nlamskin,1);
if ders == 1
    if tailflag == 1
        dexmaxdA = zeros(numel(exmax),6*length(constant.lam.ID));
        dexmaxdD = zeros(numel(exmax),6*length(constant.lam.ID));
        dexmaxdt = zeros(numel(exmax),length(constant.lam.ID));
        
        dgammamaxdA = zeros(numel(gammamax),6*length(constant.lam.ID));
        dgammamaxdD = zeros(numel(gammamax),6*length(constant.lam.ID));
        dgammamaxdt = zeros(numel(gammamax),length(constant.lam.ID));
        
        dexmindA = zeros(numel(exmin),6*length(constant.lam.ID));
        dexmindD = zeros(numel(exmin),6*length(constant.lam.ID));
        dexmindt = zeros(numel(exmin),length(constant.lam.ID));
        
        drcritdA = zeros(numel(rcrit),6*length(constant.lam.ID));
        drcritdD = zeros(numel(rcrit),6*length(constant.lam.ID));
        drcritdt = zeros(numel(rcrit),length(constant.lam.ID));
        
        if gustflag == 1
            dexmaxddv = zeros(numel(exmax),9*length(constant.lam.ID));
            dexminddv = zeros(numel(exmin),9*length(constant.lam.ID));
            dgammamaxddv = zeros(numel(gammamax),9*length(constant.lam.ID));
            drcritddv = zeros(numel(rcrit),9*length(constant.lam.ID));
        end
    end
    if morphflag == 1
        if morph.camber
            dexmaxdparam = zeros(numel(exmax),length(constant.morph.camber.loc));
            dexmindparam = zeros(numel(exmin),length(constant.morph.camber.loc));
            dgammamaxdparam = zeros(numel(gammamax),length(constant.morph.camber.loc));
            drcritdparam = zeros(numel(rcrit),length(constant.morph.camber.loc));
        end
        if morph.fold
            dexmaxdtheta = zeros(numel(exmax),length(constant.morph.fold.sec));
            dexmindtheta = zeros(numel(exmin),length(constant.morph.fold.sec));
            dgammamaxdtheta = zeros(numel(gammamax),length(constant.morph.fold.sec));
            drcritdtheta = zeros(numel(rcrit),length(constant.morph.fold.sec));
        end
        if morph.shear
            dexmaxdpsi = zeros(numel(exmax),length(constant.morph.shear.sec));
            dexmindpsi = zeros(numel(exmin),length(constant.morph.shear.sec));
            dgammamaxdpsi = zeros(numel(gammamax),length(constant.morph.shear.sec));
            drcritdpsi = zeros(numel(rcrit),length(constant.morph.shear.sec));
        end
        if morph.twist
            dexmaxdphi = zeros(numel(exmax),length(constant.morph.twist.sec));
            dexmindphi = zeros(numel(exmin),length(constant.morph.twist.sec));
            dgammamaxdphi = zeros(numel(gammamax),length(constant.morph.twist.sec));
            drcritdphi = zeros(numel(rcrit),length(constant.morph.twist.sec));
        end
        if morph.span
            dexmaxdext = zeros(numel(exmax),length(constant.morph.span.sec));
            dexmindext = zeros(numel(exmin),length(constant.morph.span.sec));
            dgammamaxdext = zeros(numel(gammamax),length(constant.morph.span.sec));
            drcritdext = zeros(numel(rcrit),length(constant.morph.span.sec));
        end
    end
end

if gustflag == 1
    rei=statics.str.re((i-1)*12+(1:12))+floc.Fl((i-1)*12+(1:12));
else
    rei=statics.str.re((i-1)*12+(1:12));
end

Ci = crossmod.C((i-1)*6+(1:6),:);
Ccrossj = crossmod.Ccross{i}{ncross};

if ders == 1 && tailflag == 1
    dCidA = crossmod.dCdA(36*(i-1)+(1:36),:);
    dCidD = crossmod.dCdD(36*(i-1)+(1:36),:);
    dCcrossjdA = crossmod.dCcrossdA{i}{ncross};
    dCcrossjdD = crossmod.dCcrossdD{i}{ncross};
end

StrainF = crossmod.StrainF{i}{ncross};
if ders == 1 && tailflag == 1
    dStraindmat = crossmod.dStraindmat{i}{ncross};
end
for j=1:2
    if j==1
        rej = -Ccrossj*(Ci\rei(1:6));
    else
        rej = Ccrossj*(Ci\rei(7:12));
    end
    
    Strain = StrainF*rej;
    
    exmat = Strain(1:6:end);
    gamma = Strain(3:6:end);
    eymat = Strain(2:6:end);
    
    exout(:,j) = exmat;
    eyout(:,j) = eymat;
    gammaout(:,j) = gamma;
    
    ecenter = 0.5*(exmat+eymat);
    eradius = sqrt((0.5*(exmat-eymat)).^2+(0.5*gamma).^2);
    
    exmaxmat(:,j) = ecenter+eradius;
    exminmat(:,j) = ecenter-eradius;
    gammamat(:,j) = 2*eradius;
    
    a0 = constant.strain.C0;
    a1 = constant.strain.C1*exmaxmat(:,j)+constant.strain.C2*exminmat(:,j);
    a2 = constant.strain.C11*exmaxmat(:,j).^2+constant.strain.C22*exminmat(:,j).^2+2*constant.strain.C12*exmaxmat(:,j).*exminmat(:,j);

    rmat1 = (-a1+sqrt(a1.^2-4*a0.*a2))./(2*a0);
    rmat2 = (-a1-sqrt(a1.^2-4*a0.*a2))./(2*a0);
    
    [rmat(:,j),indr]=max([rmat1,rmat2],[],2);
    
    for k=1:crossmod.nlam(i,ncross)
        % Check whether the laminate is a skin laminate or a stringer
        % laminate
        if crossmod.laminate_type{i}{ncross}(k)==1
            lamkloc = (crossmod.Elm2D{i}{ncross}(:,2)==k);
            Strainex = exmaxmat(:,j);
            Strainex = sort(Strainex(lamkloc),'descend');
            Strainexmax = Strainex(1:2);
            emaxloc(1) = find(exmaxmat(:,j)==Strainexmax(1),1);
            emaxloc(2) = find(exmaxmat(:,j)==Strainexmax(2),1);
            [emaxloc,ind] = sort(emaxloc,'ascend');
            locemax(4*(k-1)+2*(j-1)+(1:2),1) = emaxloc;
            exmax(4*(k-1)+2*(j-1)+(1:2),1) = Strainexmax(ind);
            
            Strainex = exminmat(:,j);
            Strainex = sort(Strainex(lamkloc),'ascend');
            Strainexmin = Strainex(1:2);
            eminloc(1) = find(exminmat(:,j)==Strainexmin(1),1);
            eminloc(2) = find(exminmat(:,j)==Strainexmin(2),1);
            [eminloc,ind] = sort(eminloc,'ascend');
            locemin(4*(k-1)+2*(j-1)+(1:2),1) = eminloc;
            exmin(4*(k-1)+2*(j-1)+(1:2),1) = Strainexmin(ind);
            
            try
                Straingamma = gammamat(:,j);
                Straingamma = sort(Straingamma(lamkloc),'descend');
                Straingammamax = Straingamma(1:4);
                gammamaxloc(1) = find(gammamat(:,j)==Straingammamax(1),1);
                gammamaxloc(2) = find(gammamat(:,j)==Straingammamax(2),1);
                gammamaxloc(3) = find(gammamat(:,j)==Straingammamax(3),1);
                gammamaxloc(4) = find(gammamat(:,j)==Straingammamax(4),1);
                [gammamaxloc,ind] = sort(gammamaxloc,'ascend');
                
                locgmax(8*(k-1)+4*(j-1)+(1:4),1) = gammamaxloc;
                gammamax(8*(k-1)+4*(j-1)+(1:4),1) = Straingammamax(ind);
            catch
                error('Please Increase cross.numel')
            end
            
            Strainr = rmat(:,j);
            Strainr = sort(Strainr(lamkloc),'descend');
            Strainr = Strainr(1:4);
            rloc(1) = find(rmat(:,j)==Strainr(1),1);
            rloc(2) = find(rmat(:,j)==Strainr(2),1);
            rloc(3) = find(rmat(:,j)==Strainr(3),1);
            rloc(4) = find(rmat(:,j)==Strainr(4),1);
            [rloc,ind] = sort(rloc,'ascend');
            locrcrit(8*(k-1)+4*(j-1)+(1:4),1) = rloc;
            rcrit(8*(k-1)+4*(j-1)+(1:4),1) = Strainr(ind);
        end
    end
    
    if ders==1
        nex = size(exmaxmat,1);

        da1dexmax = constant.strain.C1*eye(nex);
        da1dexmin = constant.strain.C2*eye(nex);
        da2dexmax = 2*constant.strain.C11*diag(exmaxmat(:,j))+2*constant.strain.C12*diag(exminmat(:,j));
        da2dexmin = 2*constant.strain.C12*diag(exmaxmat(:,j))+2*constant.strain.C22*diag(exminmat(:,j));
        
        drs1dexmax = 1/(2*a0)*(-da1dexmax+1/2./sqrt(a1(:,ones(nex,1)).^2-4*a0*a2(:,ones(nex,1))).*(2*a1(:,ones(nex,1)).*da1dexmax-4*a0*da2dexmax));
        drs1dexmin = 1/(2*a0)*(-da1dexmin+1/2./sqrt(a1(:,ones(nex,1)).^2-4*a0*a2(:,ones(nex,1))).*(2*a1(:,ones(nex,1)).*da1dexmin-4*a0*da2dexmin));
        drs2dexmax = 1/(2*a0)*(-da1dexmax-1/2./sqrt(a1(:,ones(nex,1)).^2-4*a0*a2(:,ones(nex,1))).*(2*a1(:,ones(nex,1)).*da1dexmax-4*a0*da2dexmax));
        drs2dexmin = 1/(2*a0)*(-da1dexmin-1/2./sqrt(a1(:,ones(nex,1)).^2-4*a0*a2(:,ones(nex,1))).*(2*a1(:,ones(nex,1)).*da1dexmin-4*a0*da2dexmin));

            
        if tailflag == 1
            dStraindA = zeros(size(Strain,1),6*length(constant.lam.ID));
            dStraindD = zeros(size(Strain,1),6*length(constant.lam.ID));
            dredt = Ccrossj*(Ci\statics.sens.dredt((i-1)*12+(j-1)*6+(1:6),:));
            
            if j == 1
                dredt = -dredt;
            end
            dStraindt = StrainF*dredt;
            % Conversion between DECAT_V1 and this code to select the correct
            % sensitivies
            conv = [1,4,6,2,5,3];
            
            for m = 1:length(crossmod.laminates{i}{ncross})
                % Check if laminate is actually a skin laminate or a stringer
                % laminate
                if crossmod.laminate_type{i}{ncross}(m)==1
                    for l = 1:6
                        dStraindA(:,(crossmod.laminates{i}{ncross}(m)-1)*6+l) = dStraindmat(:,6*18*(m-1)+6*(conv(l)-1)+(1:6))*rej;
                        
                        dStraindD(:,(crossmod.laminates{i}{ncross}(m)-1)*6+l) = dStraindmat(:,6*18*(m-1)+6*(12+conv(l)-1)+(1:6))*rej;
                    end
                end
            end
            

            dredA = Ccrossj*(Ci\statics.sens.dC_re((i-1)*12+(j-1)*6+(1:6),:)*crossmod.dCdA);
            for nA = 1:size(dCidA,2)
                dCcross = reshape(dCcrossjdA(:,nA),6,6)';
                dC = reshape(dCidA(:,nA),6,6)';
                dredA(:,nA) = dredA(:,nA)+dCcross*(Ci\rei(6*(j-1)+(1:6)))-Ccrossj*(Ci\dC)*(Ci\rei(6*(j-1)+(1:6)));
            end

            if j == 1
                dredA = -dredA;
            end
            
            dStraindA = dStraindA+StrainF*dredA;
            
            dredD = Ccrossj*(Ci\statics.sens.dC_re((i-1)*12+(j-1)*6+(1:6),:)*crossmod.dCdD);
            for nD = 1:size(dCidD,2)
                dCcross = reshape(dCcrossjdD(:,nD),6,6)';
                dC = reshape(dCidD(:,nD),6,6)';
                dredD(:,nD) = dredD(:,nD)+dCcross*(Ci\rei(6*(j-1)+(1:6)))-Ccrossj*(Ci\dC)*(Ci\rei(6*(j-1)+(1:6)));
            end
            
            if j == 1
                dredD = -dredD;
            end
            
            dStraindD = dStraindD+StrainF*dredD;
            
            if gustflag == 1
                dreddv = Ccrossj*(Ci\floc.dFlddv(((i-1)*12+(j-1)*6)+1:((i-1)*12+(j)*6),:));
                
                if j == 1
                    dreddv = -dreddv;
                end
                
                dStrainddv = StrainF*dreddv;
            end
            
            dexmatdA = dStraindA(1:6:end,:);
            dgammadA = dStraindA(3:6:end,:);
            deymatdA = dStraindA(2:6:end,:);
            
            dexmatdD = dStraindD(1:6:end,:);
            dgammadD = dStraindD(3:6:end,:);
            deymatdD = dStraindD(2:6:end,:);
            
            dexmatdt = dStraindt(1:6:end,:);
            dgammadt = dStraindt(3:6:end,:);
            deymatdt = dStraindt(2:6:end,:);
            
            if gustflag == 1
                dexmatddv = dStrainddv(1:6:end,:);
                dgammaddv = dStrainddv(3:6:end,:);
                deymatddv = dStrainddv(2:6:end,:);
            end
            
            dexoutdA{j} = dStraindA(1:6:end,:);
            dgammaoutdA{j} = dStraindA(3:6:end,:);
            deyoutdA{j} = dStraindA(2:6:end,:);
            
            dexoutdD{j} = dStraindD(1:6:end,:);
            dgammaoutdD{j} = dStraindD(3:6:end,:);
            deyoutdD{j} = dStraindD(2:6:end,:);
            
            dexoutdt{j} = dStraindt(1:6:end,:);
            dgammaoutdt{j} = dStraindt(3:6:end,:);
            deyoutdt{j} = dStraindt(2:6:end,:);
            
            if gustflag == 1
                dexoutddv{j} = dStrainddv(1:6:end,:);
                dgammaoutddv{j} = dStrainddv(3:6:end,:);
                deyoutddv{j} = dStrainddv(2:6:end,:);
            end
            
            decenterdA = 0.5*(dexmatdA+deymatdA);
            decenterdD = 0.5*(dexmatdD+deymatdD);
            decenterdt = 0.5*(dexmatdt+deymatdt);
            
            if gustflag == 1
                decenterddv = 0.5*(dexmatddv+deymatddv);
            end
            
            ndA = size(dexmatdA,2);
            ndD = size(dexmatdD,2);
            ndt = size(dexmatdt,2);
            if gustflag == 1
                nddv = size(dexmatddv,2);
            end
            
            deradiusdA = 1./(2*eradius(:,ones(ndA,1))).*...
                0.5.*((exmat(:,ones(ndA,1))-eymat(:,ones(ndA,1))).*(dexmatdA-deymatdA)+...
                gamma(:,ones(ndA,1)).*dgammadA);
            
            deradiusdD = 1./(2*eradius(:,ones(ndD,1))).*...
                0.5.*((exmat(:,ones(ndD,1))-eymat(:,ones(ndD,1))).*(dexmatdD-deymatdD)+...
                gamma(:,ones(ndD,1)).*dgammadD);
            
            deradiusdt = 1./(2*eradius(:,ones(ndt,1))).*...
                0.5.*((exmat(:,ones(ndt,1))-eymat(:,ones(ndt,1))).*(dexmatdt-deymatdt)+...
                gamma(:,ones(ndt,1)).*dgammadt);
            
            if gustflag == 1
                deradiusddv = 1./(2*eradius(:,ones(nddv,1))).*...
                    0.5.*((exmat(:,ones(nddv,1))-eymat(:,ones(nddv,1))).*(dexmatddv-deymatddv)+...
                    gamma(:,ones(nddv,1)).*dgammaddv);
            end
            dexmaxmatdA = decenterdA+deradiusdA;
            dexminmatdA = decenterdA-deradiusdA;
            dgammamatdA = 2*deradiusdA;
            
            dexmaxmatdD = decenterdD+deradiusdD;
            dexminmatdD = decenterdD-deradiusdD;
            dgammamatdD = 2*deradiusdD;
            
            dexmaxmatdt = decenterdt+deradiusdt;
            dexminmatdt = decenterdt-deradiusdt;
            dgammamatdt = 2*deradiusdt;
            
            if gustflag == 1
                dexmaxmatddv = decenterddv+deradiusddv;
                dexminmatddv = decenterddv-deradiusddv;
                dgammamatddv = 2*deradiusddv;
            end
            
            drmatdA(indr==1,:) = drs1dexmax(indr==1,:)*dexmaxmatdA+drs1dexmin(indr==1,:)*dexminmatdA;
            drmatdA(indr==2,:) = drs2dexmax(indr==2,:)*dexmaxmatdA+drs2dexmin(indr==2,:)*dexminmatdA;
            
            drmatdD(indr==1,:) = drs1dexmax(indr==1,:)*dexmaxmatdD+drs1dexmin(indr==1,:)*dexminmatdD;
            drmatdD(indr==2,:) = drs2dexmax(indr==2,:)*dexmaxmatdD+drs2dexmin(indr==2,:)*dexminmatdD;
            
            drmatdt(indr==1,:) = drs1dexmax(indr==1,:)*dexmaxmatdt+drs1dexmin(indr==1,:)*dexminmatdt;
            drmatdt(indr==2,:) = drs2dexmax(indr==2,:)*dexmaxmatdt+drs2dexmin(indr==2,:)*dexminmatdt;
            
            if gustflag == 1
                drmatddv(indr==1,:) = drs1dexmax(indr==1,:)*dexmaxmatddv+drs1dexmin(indr==1,:)*dexminmatddv;
                drmatddv(indr==2,:) = drs2dexmax(indr==2,:)*dexmaxmatddv+drs2dexmin(indr==2,:)*dexminmatddv;
            end
            
            for k=1:crossmod.nlam(i,ncross)
                if crossmod.laminate_type{i}{ncross}(k)==1
                    dexmaxdA((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexmaxmatdA(locemax(4*(k-1)+2*(j-1)+(1:2),1),:);
                    dexmaxdD((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexmaxmatdD(locemax(4*(k-1)+2*(j-1)+(1:2),1),:);
                    dexmaxdt((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexmaxmatdt(locemax(4*(k-1)+2*(j-1)+(1:2),1),:);
                    
                    dgammamaxdA((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = dgammamatdA(locgmax(8*(k-1)+4*(j-1)+(1:4),1),:);
                    dgammamaxdD((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = dgammamatdD(locgmax(8*(k-1)+4*(j-1)+(1:4),1),:);
                    dgammamaxdt((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = dgammamatdt(locgmax(8*(k-1)+4*(j-1)+(1:4),1),:);
                    
                    dexmindA((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexminmatdA(locemin(4*(k-1)+2*(j-1)+(1:2),1),:);
                    dexmindD((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexminmatdD(locemin(4*(k-1)+2*(j-1)+(1:2),1),:);
                    dexmindt((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexminmatdt(locemin(4*(k-1)+2*(j-1)+(1:2),1),:);
                    
                    drcritdA((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = drmatdA(locrcrit(8*(k-1)+4*(j-1)+(1:4),1),:);
                    drcritdD((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = drmatdD(locrcrit(8*(k-1)+4*(j-1)+(1:4),1),:);
                    drcritdt((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = drmatdt(locrcrit(8*(k-1)+4*(j-1)+(1:4),1),:);
                    
                    
                    if gustflag == 1
                        dexmaxddv((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexmaxmatddv(locemax(4*(k-1)+2*(j-1)+(1:2),1),:);
                        dgammamaxddv((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = dgammamatddv(locgmax(8*(k-1)+4*(j-1)+(1:4),1),:);
                        drcritddv((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = drmatddv(locrcrit(8*(k-1)+4*(j-1)+(1:4),1),:);
                        dexminddv((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexminmatddv(locemin(4*(k-1)+2*(j-1)+(1:2),1),:);
                    end
                end
            end
        end
        if morphflag == 1
            if morph.camber
                
                if gustflag == 1
                    dredparam = Ccrossj*(Ci\(statics.sens.dredparam((i-1)*12+(j-1)*6+(1:6),:)+floc.dFldparam(((i-1)*12+(j-1)*6)+1:((i-1)*12+(j)*6),:)));
                else
                    dredparam = Ccrossj*(Ci\statics.sens.dredparam((i-1)*12+(j-1)*6+(1:6),:));
                end
                if j == 1
                    dredparam = -dredparam;
                end
                
                dStraindparam = StrainF*dredparam;
                
                dexmatdparam = dStraindparam(1:6:end,:);
                dgammadparam = dStraindparam(3:6:end,:);
                deymatdparam = dStraindparam(2:6:end,:);
                
                dexoutdparam{j} = dStraindparam(1:6:end,:);
                dgammaoutdparam{j} = dStraindparam(3:6:end,:);
                deyoutdparam{j} = dStraindparam(2:6:end,:);
                
                decenterdparam = 0.5*(dexmatdparam+deymatdparam);
                
                ndparam = size(dexmatdparam,2);
                
                deradiusdparam = 1./(2*eradius(:,ones(ndparam,1))).*...
                    0.5.*((exmat(:,ones(ndparam,1))-eymat(:,ones(ndparam,1))).*(dexmatdparam-deymatdparam)+...
                    gamma(:,ones(ndparam,1)).*dgammadparam);
                
                dexmaxmatdparam = decenterdparam+deradiusdparam;
                dexminmatdparam = decenterdparam-deradiusdparam;
                dgammamatdparam = 2*deradiusdparam;
                
                drmatdparam(indr==1,:) = drs1dexmax(indr==1,:)*dexmaxmatdparam+drs1dexmin(indr==1,:)*dexminmatdparam;
                drmatdparam(indr==2,:) = drs2dexmax(indr==2,:)*dexmaxmatdparam+drs2dexmin(indr==2,:)*dexminmatdparam;

                for k=1:crossmod.nlam(i,ncross)
                    if crossmod.laminate_type{i}{ncross}(k)==1
                        dexmaxdparam((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexmaxmatdparam(locemax(4*(k-1)+2*(j-1)+(1:2),1),:);
                        dgammamaxdparam((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = dgammamatdparam(locgmax(8*(k-1)+4*(j-1)+(1:4),1),:);
                        dexmindparam((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexminmatdparam(locemin(4*(k-1)+2*(j-1)+(1:2),1),:);
                        drcritdparam((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = drmatdparam(locrcrit(8*(k-1)+4*(j-1)+(1:4),1),:);
                    end
                end
            end
            if morph.twist
                
                if gustflag == 1
                    dredphi = Ccrossj*(Ci\(statics.sens.dredphi((i-1)*12+(j-1)*6+(1:6),:)+floc.dFldphi(((i-1)*12+(j-1)*6)+1:((i-1)*12+(j)*6),:)));
                else
                    dredphi = Ccrossj*(Ci\statics.sens.dredphi((i-1)*12+(j-1)*6+(1:6),:));
                end
                if j == 1
                    dredphi = -dredphi;
                end
                
                dStraindphi = StrainF*dredphi;
                
                dexmatdphi = dStraindphi(1:6:end,:);
                dgammadphi = dStraindphi(3:6:end,:);
                deymatdphi = dStraindphi(2:6:end,:);
                
                dexoutdphi{j} = dStraindphi(1:6:end,:);
                dgammaoutdphi{j} = dStraindphi(3:6:end,:);
                deyoutdphi{j} = dStraindphi(2:6:end,:);
                
                decenterdphi = 0.5*(dexmatdphi+deymatdphi);
                
                ndphi = size(dexmatdphi,2);
                
                deradiusdphi = 1./(2*eradius(:,ones(ndphi,1))).*...
                    0.5.*((exmat(:,ones(ndphi,1))-eymat(:,ones(ndphi,1))).*(dexmatdphi-deymatdphi)+...
                    gamma(:,ones(ndphi,1)).*dgammadphi);
                
                dexmaxmatdphi = decenterdphi+deradiusdphi;
                dexminmatdphi = decenterdphi-deradiusdphi;
                dgammamatdphi = 2*deradiusdphi;
                
                drmatdphi(indr==1,:) = drs1dexmax(indr==1,:)*dexmaxmatdphi+drs1dexmin(indr==1,:)*dexminmatdphi;
                drmatdphi(indr==2,:) = drs2dexmax(indr==2,:)*dexmaxmatdphi+drs2dexmin(indr==2,:)*dexminmatdphi;
                
                for k=1:crossmod.nlam(i,ncross)
                    if crossmod.laminate_type{i}{ncross}(k)==1
                        dexmaxdphi((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexmaxmatdphi(locemax(4*(k-1)+2*(j-1)+(1:2),1),:);
                        dgammamaxdphi((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = dgammamatdphi(locgmax(8*(k-1)+4*(j-1)+(1:4),1),:);
                        dexmindphi((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexminmatdphi(locemin(4*(k-1)+2*(j-1)+(1:2),1),:);
                        drcritdphi((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = drmatdphi(locrcrit(8*(k-1)+4*(j-1)+(1:4),1),:);
                    end
                end
            end
            if morph.shear
                
                if gustflag == 1
                    dredpsi = Ccrossj*(Ci\(statics.sens.dredpsi((i-1)*12+(j-1)*6+(1:6),:)+floc.dFldpsi(((i-1)*12+(j-1)*6)+1:((i-1)*12+(j)*6),:)));
                else
                    dredpsi = Ccrossj*(Ci\statics.sens.dredpsi((i-1)*12+(j-1)*6+(1:6),:));
                end
                if j == 1
                    dredpsi = -dredpsi;
                end
                
                dStraindpsi = StrainF*dredpsi;
                
                dexmatdpsi = dStraindpsi(1:6:end,:);
                dgammadpsi = dStraindpsi(3:6:end,:);
                deymatdpsi = dStraindpsi(2:6:end,:);
                
                dexoutdpsi{j} = dStraindpsi(1:6:end,:);
                dgammaoutdpsi{j} = dStraindpsi(3:6:end,:);
                deyoutdpsi{j} = dStraindpsi(2:6:end,:);
                
                decenterdpsi = 0.5*(dexmatdpsi+deymatdpsi);
                
                ndpsi = size(dexmatdpsi,2);
                
                deradiusdpsi = 1./(2*eradius(:,ones(ndpsi,1))).*...
                    0.5.*((exmat(:,ones(ndpsi,1))-eymat(:,ones(ndpsi,1))).*(dexmatdpsi-deymatdpsi)+...
                    gamma(:,ones(ndpsi,1)).*dgammadpsi);
                
                dexmaxmatdpsi = decenterdpsi+deradiusdpsi;
                dexminmatdpsi = decenterdpsi-deradiusdpsi;
                dgammamatdpsi = 2*deradiusdpsi;
                
                drmatdpsi(indr==1,:) = drs1dexmax(indr==1,:)*dexmaxmatdpsi+drs1dexmin(indr==1,:)*dexminmatdpsi;
                drmatdpsi(indr==2,:) = drs2dexmax(indr==2,:)*dexmaxmatdpsi+drs2dexmin(indr==2,:)*dexminmatdpsi;
                
                for k=1:crossmod.nlam(i,ncross)
                    if crossmod.laminate_type{i}{ncross}(k)==1
                        dexmaxdpsi((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexmaxmatdpsi(locemax(4*(k-1)+2*(j-1)+(1:2),1),:);
                        dgammamaxdpsi((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = dgammamatdpsi(locgmax(8*(k-1)+4*(j-1)+(1:4),1),:);
                        dexmindpsi((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexminmatdpsi(locemin(4*(k-1)+2*(j-1)+(1:2),1),:);
                        drcritdpsi((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = drmatdpsi(locrcrit(8*(k-1)+4*(j-1)+(1:4),1),:);
                    end
                end
            end
            if morph.span
                
                if gustflag == 1
                    dredext = Ccrossj*(Ci\(statics.sens.dredext((i-1)*12+(j-1)*6+(1:6),:)+floc.dFldext(((i-1)*12+(j-1)*6)+1:((i-1)*12+(j)*6),:)));
                else
                    dredext = Ccrossj*(Ci\statics.sens.dredext((i-1)*12+(j-1)*6+(1:6),:));
                end
                if j == 1
                    dredext = -dredext;
                end
                
                dStraindext = StrainF*dredext;
                
                dexmatdext = dStraindext(1:6:end,:);
                dgammadext = dStraindext(3:6:end,:);
                deymatdext = dStraindext(2:6:end,:);
                
                dexoutdext{j} = dStraindext(1:6:end,:);
                dgammaoutdext{j} = dStraindext(3:6:end,:);
                deyoutdext{j} = dStraindext(2:6:end,:);
                
                decenterdext = 0.5*(dexmatdext+deymatdext);
                
                ndext = size(dexmatdext,2);
                
                deradiusdext = 1./(2*eradius(:,ones(ndext,1))).*...
                    0.5.*((exmat(:,ones(ndext,1))-eymat(:,ones(ndext,1))).*(dexmatdext-deymatdext)+...
                    gamma(:,ones(ndext,1)).*dgammadext);
                
                dexmaxmatdext = decenterdext+deradiusdext;
                dexminmatdext = decenterdext-deradiusdext;
                dgammamatdext = 2*deradiusdext;
                
                drmatdext(indr==1,:) = drs1dexmax(indr==1,:)*dexmaxmatdext+drs1dexmin(indr==1,:)*dexminmatdext;
                drmatdext(indr==2,:) = drs2dexmax(indr==2,:)*dexmaxmatdext+drs2dexmin(indr==2,:)*dexminmatdext;
                
                for k=1:crossmod.nlam(i,ncross)
                    if crossmod.laminate_type{i}{ncross}(k)==1
                        dexmaxdext((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexmaxmatdext(locemax(4*(k-1)+2*(j-1)+(1:2),1),:);
                        dgammamaxdext((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = dgammamatdext(locgmax(8*(k-1)+4*(j-1)+(1:4),1),:);
                        dexmindext((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexminmatdext(locemin(4*(k-1)+2*(j-1)+(1:2),1),:);
                        drcritdext((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = drmatdext(locrcrit(8*(k-1)+4*(j-1)+(1:4),1),:);
                    end
                end
            end
            if morph.fold
                
                if gustflag == 1
                    dredtheta = Ccrossj*(Ci\(statics.sens.dredtheta((i-1)*12+(j-1)*6+(1:6),:)+floc.dFldtheta(((i-1)*12+(j-1)*6)+1:((i-1)*12+(j)*6),:)));
                else
                    dredtheta = Ccrossj*(Ci\statics.sens.dredtheta((i-1)*12+(j-1)*6+(1:6),:));
                end
                if j == 1
                    dredtheta = -dredtheta;
                end
                
                dStraindtheta = StrainF*dredtheta;
                
                dexmatdtheta = dStraindtheta(1:6:end,:);
                dgammadtheta = dStraindtheta(3:6:end,:);
                deymatdtheta = dStraindtheta(2:6:end,:);
                
                dexoutdtheta{j} = dStraindtheta(1:6:end,:);
                dgammaoutdtheta{j} = dStraindtheta(3:6:end,:);
                deyoutdtheta{j} = dStraindtheta(2:6:end,:);
                
                decenterdtheta = 0.5*(dexmatdtheta+deymatdtheta);
                
                ndtheta = size(dexmatdtheta,2);
                
                deradiusdtheta = 1./(2*eradius(:,ones(ndtheta,1))).*...
                    0.5.*((exmat(:,ones(ndtheta,1))-eymat(:,ones(ndtheta,1))).*(dexmatdtheta-deymatdtheta)+...
                    gamma(:,ones(ndtheta,1)).*dgammadtheta);
                
                dexmaxmatdtheta = decenterdtheta+deradiusdtheta;
                dexminmatdtheta = decenterdtheta-deradiusdtheta;
                dgammamatdtheta = 2*deradiusdtheta;
                
                drmatdtheta(indr==1,:) = drs1dexmax(indr==1,:)*dexmaxmatdtheta+drs1dexmin(indr==1,:)*dexminmatdtheta;
                drmatdtheta(indr==2,:) = drs2dexmax(indr==2,:)*dexmaxmatdtheta+drs2dexmin(indr==2,:)*dexminmatdtheta;
                
                for k=1:crossmod.nlam(i,ncross)
                    if crossmod.laminate_type{i}{ncross}(k)==1
                        dexmaxdtheta((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexmaxmatdtheta(locemax(4*(k-1)+2*(j-1)+(1:2),1),:);
                        dgammamaxdtheta((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = dgammamatdtheta(locgmax(8*(k-1)+4*(j-1)+(1:4),1),:);
                        dexmindtheta((4*(k-1)+2*(j-1))+1:(4*(k-1)+2*(j)),:) = dexminmatdtheta(locemin(4*(k-1)+2*(j-1)+(1:2),1),:);
                        drcritdtheta((8*(k-1)+4*(j-1))+1:(8*(k-1)+4*(j)),:) = drmatdtheta(locrcrit(8*(k-1)+4*(j-1)+(1:4),1),:);
                    end
                end
            end
        end
    end
end



% Filter out stringer strain results
stringerlam = find(crossmod.laminate_type{i}{ncross}==2);
stringerelm= [];
for j=1:length(stringerlam)
    stringerelm = [stringerelm;find(crossmod.Elm2D{i}{ncross}(:,2)==stringerlam(j))];
end
exmaxmat(stringerelm,:) = [];
exminmat(stringerelm,:) = [];
gammamat(stringerelm,:) = [];


strainraw.exout = exout;
strainraw.eyout = eyout;
strainraw.gammaout = gammaout;
if ders == 1
    if tailflag == 1
        strainraw.dexoutdA = dexoutdA;
        strainraw.deyoutdA = deyoutdA;
        strainraw.dgammaoutdA = dgammaoutdA;
        
        strainraw.dexoutdD = dexoutdD;
        strainraw.deyoutdD = deyoutdD;
        strainraw.dgammaoutdD = dgammaoutdD;
        
        strainraw.dexoutdt = dexoutdt;
        strainraw.deyoutdt = deyoutdt;
        strainraw.dgammaoutdt = dgammaoutdt;
        
        if gustflag == 1
            strainraw.dexoutddv = dexoutddv;
            strainraw.deyoutddv = deyoutddv;
            strainraw.dgammaoutddv = dgammaoutddv;
        end
    end
    if morphflag == 1
        if morph.camber
            strainraw.dexoutdparam = dexoutdparam;
            strainraw.deyoutdparam = deyoutdparam;
            strainraw.dgammaoutdparam = dgammaoutdparam;
        end
        if morph.twist
            strainraw.dexoutdphi = dexoutdphi;
            strainraw.deyoutdphi = deyoutdphi;
            strainraw.dgammaoutdphi = dgammaoutdphi;
        end
        if morph.shear
            strainraw.dexoutdpsi = dexoutdpsi;
            strainraw.deyoutdpsi = deyoutdpsi;
            strainraw.dgammaoutdpsi = dgammaoutdpsi;
        end
        if morph.span
            strainraw.dexoutdext = dexoutdext;
            strainraw.deyoutdext = deyoutdext;
            strainraw.dgammaoutdext = dgammaoutdext;
        end
        if morph.fold
            strainraw.dexoutdtheta = dexoutdtheta;
            strainraw.deyoutdtheta = deyoutdtheta;
            strainraw.dgammaoutdtheta = dgammaoutdtheta;
        end
    end
end
constr.exmax = exmax;
constr.exmin = exmin;
constr.gammamax = gammamax;
constr.rcrit = rcrit;
constr.locemax = locemax;
constr.locemin = locemin;
constr.locgmax = locgmax;
constr.locrcrit = locrcrit;
if ders == 1
    if tailflag == 1
        constr.dexmaxdA = dexmaxdA;
        constr.dexmindA = dexmindA;
        constr.dgammamaxdA = dgammamaxdA;
        constr.drcritdA = drcritdA;
        constr.dexmaxdD = dexmaxdD;
        constr.dexmindD = dexmindD;
        constr.dgammamaxdD = dgammamaxdD;
        constr.drcritdD = drcritdD;
        constr.dexmaxdt = dexmaxdt;
        constr.dexmindt = dexmindt;
        constr.dgammamaxdt = dgammamaxdt;
        constr.drcritdt = drcritdt;
        if gustflag == 1
            constr.dexmaxddv = dexmaxddv;
            constr.dexminddv = dexminddv;
            constr.dgammamaxddv = dgammamaxddv;
            constr.drcritddv = drcritddv;
        end
    end
    if morphflag == 1
        if morph.camber
            constr.dexmaxdparam = dexmaxdparam;
            constr.dexmindparam = dexmindparam;
            constr.dgammamaxdparam = dgammamaxdparam;
            constr.drcritdparam = drcritdparam;
        end
        if morph.twist
            constr.dexmaxdphi = dexmaxdphi;
            constr.dexmindphi = dexmindphi;
            constr.dgammamaxdphi = dgammamaxdphi;
            constr.drcritdphi = drcritdphi;
        end
        if morph.shear
            constr.dexmaxdpsi = dexmaxdpsi;
            constr.dexmindpsi = dexmindpsi;
            constr.dgammamaxdpsi = dgammamaxdpsi;
            constr.drcritdpsi = drcritdpsi;
        end
        if morph.span
            constr.dexmaxdext = dexmaxdext;
            constr.dexmindext = dexmindext;
            constr.dgammamaxdext = dgammamaxdext;
            constr.drcritdext = drcritdext;
        end
        if morph.fold
            constr.dexmaxdtheta = dexmaxdtheta;
            constr.dexmindtheta = dexmindtheta;
            constr.dgammamaxdtheta = dgammamaxdtheta;
            constr.drcritdtheta = drcritdtheta;
        end
    end
end



