function [constant] = gravupd(constant,statics,indgrav,ders,tailflag)

%% Structural mass
if ders == 1 && tailflag == 1
    dMmagdt = zeros(constant.str.Ns,length(constant.lam.ID));
end
for i = 1:constant.str.Ns
    Mmag(i) = statics.cross.mA(i)*norm(statics.morph.span.xyz(3*i+(1:3),end)-statics.morph.span.xyz(3*(i-1)+(1:3),end));
    dMmagdext(i,:) = statics.cross.mA(i)*dnorm(statics.morph.span.xyz(3*i+(1:3),end)-statics.morph.span.xyz(3*(i-1)+(1:3),end),statics.morph.span.dxyzdext(3*i+(1:3),:)-statics.morph.span.dxyzdext(3*(i-1)+(1:3),:));

    if ders == 1 && tailflag == 1
        dMmagdt(i,:) = statics.cross.dmAdt(i,:)*norm(statics.morph.span.xyz(3*i+(1:3),end)-statics.morph.span.xyz(3*(i-1)+(1:3),end));
    end
end

%% Gravity

if ders == 1 && tailflag == 1
    constant.fext.dmagnitudedt{indgrav} = zeros(6*constant.str.Ns,length(constant.lam.ID));
end

constant.fext.dmagnitudedext{indgrav} = zeros(6*constant.str.Ns,constant.str.Ns);

for i = 1:constant.str.Ns
    Fmag = Mmag(i)*9.81;
    constant.fext.magnitude{indgrav}(i,:) = constant.general.nz*[0,0,-Fmag,0,0,0];
    constant.fext.dmagnitudedext{indgrav}(6*(i-1)+3,:) = -constant.general.nz*9.81*dMmagdext(i,:);

    if ders == 1 && tailflag == 1
        constant.fext.dmagnitudedt{indgrav}(6*(i-1)+3,:) = -constant.general.nz*9.81*dMmagdt(i,:);
    end
end

function [dnorma] = dnorm(a,da)
    dnorma = 1/2/norm(a)*(sum(2*a(:,ones(1,size(da,2))).*da));
