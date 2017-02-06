% clear all
close all
% clc
% format short e
% 
% cd ..
% load sens_test_buckl_rand.mat
% cd ('buckling')

ders = 1;
tic
[buckl] = buckl_comp(constant,lampar,crossmod,strainraw,ders,gustflag,tailflag,morph);
toc

%%

ddv = 1e-5*(0.5-1*rand(length(dv),1));
dparam = 1e-5*(0.5-1*rand(length(constant.morph.camber.param),1)).*constant.morph.camber.loc;
dphi = 1e-5*(0.5-1*rand(length(constant.morph.twist.angle),1)).*constant.morph.twist.sec;
dpsi = 1e-5*(0.5-1*rand(length(constant.morph.shear.angle),1)).*constant.morph.shear.sec;
dext = 1e-5*(0.5-1*rand(length(constant.morph.span.ext),1)).*constant.morph.span.sec;
dtheta = 1e-5*(0.5-1*rand(length(constant.morph.fold.angle),1)).*constant.morph.fold.sec;

row    = [1,1,1,2,2,3];     % Row number of independent A and D matrices elements
column = [1,2,3,2,3,3];     % Column number of independent A and D matrices elements

% Perturb A through ddv
dA = lampar.dAddv*ddv;
dD = lampar.dDddv*ddv;
dt = lampar.dtddv*ddv;


lampar2 = lampar;
for i=1:size(lampar.A,1)/3
   dAloc = dA(6*(i-1)+(1:6));
   dAmat = [dAloc(1),dAloc(2),dAloc(3);
            dAloc(2),dAloc(4),dAloc(5);
            dAloc(3),dAloc(5),dAloc(6)];
   lampar2.A(3*(i-1)+(1:3),1:3) = lampar.A(3*(i-1)+(1:3),1:3) + dAmat;
   
   dDloc = dD(6*(i-1)+(1:6));
   dDmat  = [dDloc(1),dDloc(2),dDloc(3);
            dDloc(2),dDloc(4),dDloc(5);
            dDloc(3),dDloc(5),dDloc(6)];
   lampar2.D(3*(i-1)+(1:3),1:3) = lampar.D(3*(i-1)+(1:3),1:3) + dDmat;
end
lampar2.t = lampar.t + dt';

strainraw2=strainraw;

for i=1:length(strainraw)
    for j=1:size(strainraw{1}.exout,2)
        strainraw2{i}.exout(:,j) = strainraw{i}.exout(:,j)+strainraw{i}.dexoutdA{j}*dA+strainraw{i}.dexoutdD{j}*dD+strainraw{i}.dexoutdt{j}*dt+strainraw{i}.dexoutddv{j}*ddv...
            +strainraw{i}.dexoutdparam{j}*dparam+strainraw{i}.dexoutdphi{j}*dphi+strainraw{i}.dexoutdpsi{j}*dpsi+strainraw{i}.dexoutdext{j}*dext+strainraw{i}.dexoutdtheta{j}*dtheta;
        strainraw2{i}.eyout(:,j) = strainraw{i}.eyout(:,j)+strainraw{i}.deyoutdA{j}*dA+strainraw{i}.deyoutdD{j}*dD+strainraw{i}.deyoutdt{j}*dt+strainraw{i}.deyoutddv{j}*ddv...
            +strainraw{i}.deyoutdparam{j}*dparam+strainraw{i}.deyoutdphi{j}*dphi+strainraw{i}.deyoutdpsi{j}*dpsi+strainraw{i}.deyoutdext{j}*dext+strainraw{i}.deyoutdtheta{j}*dtheta;
        strainraw2{i}.gammaout(:,j) = strainraw{i}.gammaout(:,j)+strainraw{i}.dgammaoutdA{j}*dA+strainraw{i}.dgammaoutdD{j}*dD+strainraw{i}.dgammaoutdt{j}*dt+strainraw{i}.dgammaoutddv{j}*ddv...
            +strainraw{i}.dgammaoutdparam{j}*dparam+strainraw{i}.dgammaoutdphi{j}*dphi+strainraw{i}.dgammaoutdpsi{j}*dpsi+strainraw{i}.dgammaoutdext{j}*dext+strainraw{i}.dgammaoutdtheta{j}*dtheta;
    end
end

ders = 0;
tic
[buckl2] =  buckl_comp(constant,lampar2,crossmod,strainraw2,ders,gustflag,tailflag,morph);
toc

f1 = reshape((buckl.r)',[],1);
f2 = reshape((buckl2.r)',[],1);
df1 = reshape((buckl2.r-buckl.r)',[],1);
df2 = buckl.drdA*dA+buckl.drdD*dD+buckl.drdt*dt+buckl.drddv*ddv+buckl.drdparam*dparam+buckl.drdphi*dphi+buckl.drdpsi*dpsi+buckl.drdext*dext+buckl.drdtheta*dtheta;

[(1:length(df1))',f1,f2,df1,df2,(df2-df1)./df1]

figure
plot((df2-df1)./df1)




