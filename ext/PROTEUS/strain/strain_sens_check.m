

ders = 1;

if morphflag == 1
    [strain,strainraw] = strain_comp_single_gust_fail_index(constant,crossmod,statics3,dynamics.floc{1},ders,gustflag,tailflag,morph);
else
    [strain,strainraw] = strain_comp_single_gust_fail_index(constant,crossmod,statics3,dynamics.floc{1},ders,gustflag,tailflag);
end
fprintf('\n Strain 1;')

%%
ddv = 1e-5*(0.5-1*rand(length(dv),1));
% dparam = 1e-5*(0.5-1*rand(length(constant.morph.camber.param),1)).*constant.morph.camber.loc;
% dphi = 1e-5*(0.5-1*rand(length(constant.morph.twist.angle),1)).*constant.morph.twist.sec;
% dpsi = 1e-5*(0.5-1*rand(length(constant.morph.shear.angle),1)).*constant.morph.shear.sec;
% dext = 1e-5*(0.5-1*rand(length(constant.morph.span.ext),1)).*constant.morph.span.sec;
% dtheta = 1e-5*(0.5-1*rand(length(constant.morph.fold.angle),1)).*constant.morph.fold.sec;


statics32=statics3;
crossmod2 = crossmod;
dynamics2 = dynamics;

dA = lampar.dAddv*ddv;
dD = lampar.dDddv*ddv;
dt = lampar.dtddv*ddv;

dC = crossmod.dCdA*dA+crossmod.dCdD*dD;
dre = statics3.sens.dC_re*dC+statics3.sens.dredt*dt;%+statics3.sens.dredparam*dparam+statics3.sens.dredphi*dphi+statics3.sens.dredpsi*dpsi+statics3.sens.dredtheta*dtheta+statics3.sens.dredext*dext;
dFl = dynamics.floc{1}.dFlddv*ddv;%+dynamics.floc.dFlmaxdparam*dparam+dynamics.floc.dFlmaxdphi*dphi+dynamics.floc.dFlmaxdpsi*dpsi+dynamics.floc.dFlmaxdtheta*dtheta+dynamics.floc.dFlmaxdext*dext;
% dFlmin = dynamics.floc.dFlminddv*ddv+dynamics.floc.dFlmindparam*dparam+dynamics.floc.dFlmindphi*dphi+dynamics.floc.dFlmindpsi*dpsi+dynamics.floc.dFlmindtheta*dtheta+dynamics.floc.dFlmindext*dext;

crossmod2.C = crossmod.C+reshape(dC,6,[])';
statics32.str.re = statics3.str.re+dre;
dynamics2.floc{1}.Fl = dynamics.floc{1}.Fl+dFl;
% dynamics2.floc.Flmin = dynamics.floc.Flmin+reshape(dFlmin,4,[])';

conv = [1,4,6,2,5,3];
for i=1:length(crossmod.Ccross)
    for j =1:length(crossmod.Ccross{i})
        dCcross = crossmod.dCcrossdA{i}{j}*dA+crossmod.dCcrossdD{i}{j}*dD;
        crossmod2.Ccross{i}{j} = crossmod.Ccross{i}{j}+reshape(dCcross,6,6)';
        for m = 1:length(crossmod.laminates{i}{j})
            % Check if laminate is actually a skin laminate or a stringer
            % laminate
            if crossmod.laminate_type{i}{j}(m)==1
                for l = 1:6
                    dStrain = crossmod.dStraindmat{i}{j}(:,6*18*(m-1)+6*(conv(l)-1)+(1:6))*dA((crossmod.laminates{i}{j}(m)-1)*6+l)+crossmod.dStraindmat{i}{j}(:,6*18*(m-1)+6*(12+conv(l)-1)+(1:6))*dD((crossmod.laminates{i}{j}(m)-1)*6+l);
                    crossmod2.StrainF{i}{j} = crossmod2.StrainF{i}{j}+dStrain;
                end
            end
        end
    end
end

ders = 0;

if morphflag == 1
    [strain2,strainraw2] = strain_comp_single_gust_fail_index(constant,crossmod2,statics32,dynamics2.floc{1},ders,gustflag,tailflag,morph);
else
    [strain2,strainraw2] = strain_comp_single_gust_fail_index(constant,crossmod2,statics32,dynamics2.floc{1},ders,gustflag,tailflag);
end
fprintf('\n Strain 2;')

%%

out.f=strain.rcrit;
out2.f = strain2.rcrit;
out.dfdA = strain.drcritdA;
out.dfdD = strain.drcritdD;
out.dfdt = strain.drcritdt;
out.dfddv = strain.drcritddv;
% out.dfdparam = strain.dgammamaxdparam;
% out.dfdphi = strain.dgammamaxdphi;
% out.dfdpsi = strain.dgammamaxdpsi;
% out.dfdext = strain.dgammamaxdext;
% out.dfdtheta = strain.dgammamaxdtheta;

% out.f=strain.exmax;
% out2.f = strain2.exmax;
% out.dfdA = strain.dexmaxdA;
% out.dfdD = strain.dexmaxdD;
% out.dfdt = strain.dexmaxdt;
% out.dfddv = strain.dexmaxddv;
% out.dfdparam = strain.dexmaxdparam;
% out.dfdphi = strain.dexmaxdphi;
% out.dfdpsi = strain.dexmaxdpsi;
% out.dfdext = strain.dexmaxdext;
% out.dfdtheta = strain.dexmaxdtheta;

% out.f=strain.exmin;
% out2.f = strain2.exmin;
% out.dfdA = strain.dexmindA;
% out.dfdD = strain.dexmindD;
% out.dfdt = strain.dexmindt;
% out.dfddv = strain.dexminddv;
% out.dfdparam = strain.dexmindparam;
% out.dfdphi = strain.dexmindphi;
% out.dfdpsi = strain.dexmindpsi;
% out.dfdext = strain.dexmindext;
% out.dfdtheta = strain.dexmindtheta;

% out.f=strainraw{12}.gammaout(:,14);
% out2.f = strainraw2{12}.gammaout(:,14);
% out.dfdA = strainraw{12}.dgammaoutdA{14};
% out.dfdD = strainraw{12}.dgammaoutdD{14};
% out.dfdt = strainraw{12}.dgammaoutdt{14};
% out.dfddv = strainraw{12}.dgammaoutddv{14};
% out.dfdparam = strainraw{12}.dgammaoutdparam{14};
% out.dfdphi = strainraw{12}.dgammaoutdphi{14};
% out.dfdpsi = strainraw{12}.dgammaoutdpsi{14};
% out.dfdext = strainraw{12}.dgammaoutdext{14};
% out.dfdtheta = strainraw{12}.dgammaoutdtheta{14};

% out.f=strainraw{12}.exout(:,2);
% out2.f = strainraw2{12}.exout(:,2);
% out.dfdA = strainraw{12}.dexoutdA{2};
% out.dfdD = strainraw{12}.dexoutdD{2};
% out.dfdt = strainraw{12}.dexoutdt{2};
% out.dfddv = strainraw{12}.dexoutddv{2};
% out.dfdparam = strainraw{12}.dexoutdparam{2};
% out.dfdphi = strainraw{12}.dexoutdphi{2};
% out.dfdpsi = strainraw{12}.dexoutdpsi{2};
% out.dfdext = strainraw{12}.dexoutdext{2};
% out.dfdtheta = strainraw{12}.dexoutdtheta{2};

% out.f=strainraw{12}.eyout(:,2);
% out2.f = strainraw2{12}.eyout(:,2);
% out.dfdA = strainraw{12}.deyoutdA{2};
% out.dfdD = strainraw{12}.deyoutdD{2};
% out.dfdt = strainraw{12}.deyoutdt{2};
% out.dfddv = strainraw{12}.deyoutddv{2};
% out.dfdparam = strainraw{12}.deyoutdparam{2};
% out.dfdphi = strainraw{12}.deyoutdphi{2};
% out.dfdpsi = strainraw{12}.deyoutdpsi{2};
% out.dfdext = strainraw{12}.deyoutdext{2};
% out.dfdtheta = strainraw{12}.deyoutdtheta{2};


E1 = reshape((out.f)',[],1);
E2 = reshape((out2.f)',[],1); 
dE1 = full(out.dfdA*dA+out.dfdD*dD+out.dfdt*dt+out.dfddv*ddv);%+out.dfdparam*dparam+out.dfdphi*dphi+out.dfdpsi*dpsi+out.dfdext*dext+out.dfdtheta*dtheta);
% dE1 = full(out.dfdA*dA+out.dfdD*dD);
% dE1 = full(out.dfddv*ddv);
% dE1 = full(out.dfdC*dC+out.dfdt*dt);
dE2 = reshape((out2.f-out.f)',[],1);
% E1(dE2==0) = [];
% dE1(dE2==0) = [];
% dE2(dE2==0) = [];
% [~,idx] = sort(abs((dE2-dE1)./dE1),'ascend');
idx = (1:length(dE1));
[E1(idx),dE2(idx),dE1(idx),(dE2(idx)-dE1(idx))./dE1(idx)]

max(abs((dE2-dE1)./dE1))