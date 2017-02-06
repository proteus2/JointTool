cd('cross_mod')
[crossmod,stringer] = cross_mod(constant,lampar,stringer,1,ders,tailflag,morphflag);
cd(curdir)


ddv = 1e-5*(0.5-1*rand(length(dv),1));
lampar2 = lampar;

dA = lampar.dAddv*ddv;
dD = 0*lampar.dDddv*ddv;
dt = lampar.dtddv*ddv;

for i=1:size(lampar.A,1)/3
    dAmat = [dA(6*(i-1)+1),dA(6*(i-1)+2),dA(6*(i-1)+3)
        dA(6*(i-1)+2),dA(6*(i-1)+4),dA(6*(i-1)+5)
        dA(6*(i-1)+3),dA(6*(i-1)+5),dA(6*(i-1)+6)];
    
    dDmat = [dD(6*(i-1)+1),dD(6*(i-1)+2),dD(6*(i-1)+3)
        dD(6*(i-1)+2),dD(6*(i-1)+4),dD(6*(i-1)+5)
        dD(6*(i-1)+3),dD(6*(i-1)+5),dD(6*(i-1)+6)];
    lampar2.A(3*(i-1)+(1:3),:) = lampar.A(3*(i-1)+(1:3),:) + dAmat;
    lampar2.D(3*(i-1)+(1:3),:) = lampar.D(3*(i-1)+(1:3),:) + dDmat;
end

lampar2.t = lampar.t + dt';

cd('cross_mod')
[crossmod2,stringer] = cross_mod(constant,lampar2,stringer,1,0,tailflag,morphflag);
cd(curdir)

dC1 = reshape((crossmod2.C-crossmod.C)',[],1);
dC2 = crossmod.dCdA*dA+crossmod.dCdD*dD;

[dC1,dC2,(dC2-dC1)./dC1]

dC1 = reshape((crossmod2.mA-crossmod.mA)',[],1);
dC2 = crossmod.dmAdt*dt;

[dC1,dC2,(dC2-dC1)./dC1]

dC1 = reshape((crossmod2.mQ-crossmod.mQ)',[],1);
dC2 = crossmod.dmQdt*dt;

[dC1,dC2,(dC2-dC1)./dC1]

dC1 = reshape((crossmod2.mI-crossmod.mI)',[],1);
dC2 = crossmod.dmIdt*dt;

[dC1,dC2,(dC2-dC1)./dC1]

dC1 = reshape((crossmod2.Ccross{5}{2}-crossmod.Ccross{5}{2})',[],1);
dC2 = crossmod.dCcrossdA{5}{2}*dA+crossmod.dCcrossdD{5}{2}*dD;

[dC1,dC2,(dC2-dC1)./dC1]

rej = 100*rand(6,1);
i = 1;
ncross = 1;
conv = [1,4,6,2,5,3];
Strain = crossmod.StrainF{i}{ncross}*rej;
Strain2 = crossmod2.StrainF{i}{ncross}*rej;
dStraindA = zeros(size(Strain,1),6*length(constant.lam.ID));
dStraindD = zeros(size(Strain,1),6*length(constant.lam.ID));
for m = 1:length(crossmod.laminates{i}{ncross})
    % Check if laminate is actually a skin laminate or a stringer
    % laminate
    if crossmod.laminate_type{i}{ncross}(m)==1
        for l = 1:6
            dStraindA(:,(crossmod.laminates{i}{ncross}(m)-1)*6+l) = crossmod.dStraindmat{i}{ncross}(:,6*18*(m-1)+6*(conv(l)-1)+(1:6))*rej;
            dStraindD(:,(crossmod.laminates{i}{ncross}(m)-1)*6+l) = crossmod.dStraindmat{i}{ncross}(:,6*18*(m-1)+6*(12+conv(l)-1)+(1:6))*rej;
        end
    end
end

dC1 = reshape((Strain2-Strain)',[],1);
dC2 = dStraindA*dA+dStraindD*dD;

[dC1,dC2,(dC2-dC1)./dC1]
