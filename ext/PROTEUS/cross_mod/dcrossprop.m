function [mA,mQ2,mQ3,mI11,mI22,mI23,mI33,dmAdt,dmQ2dt,dmQ3dt,dmI11dt,dmI22dt,dmI23dt,dmI33dt] = dcrossprop(xy,t,elmloc,rho,ders)

Nelmloc = size(elmloc,1);

SquareXY = (xy(elmloc(:,2),:)-xy(elmloc(:,1),:)).^2;
xy_norm  = (sqrt(sum( SquareXY ,2)));
Aloc     = t.*xy_norm;                              % correct
if ders == 1
    dmAdt = rho.*xy_norm; % vectorised           % correct
end
alpha = atan((xy(elmloc(:,2),2)-xy(elmloc(:,1),2))./(xy(elmloc(:,2),1)-xy(elmloc(:,1),1)));

mask = alpha==0; %& [(xy(elmloc(:,2),1)-xy(elmloc(:,1),1))<0];
alpha(mask) = eps;

tx(sin(alpha)~=0,:)  = t./sin(alpha(sin(alpha)~=0)); % if error, then one alpha == 0
dtx(sin(alpha)~=0,:) = 1./sin(alpha(sin(alpha)~=0));

ty(cos(alpha)~=0,:)  = t./cos(alpha(sin(alpha)~=0));
dty(cos(alpha)~=0,:) = 1./cos(alpha(sin(alpha)~=0));

I22loc = 1/12*tx.*(xy(elmloc(:,2),2)-xy(elmloc(:,1),2)).^3; % correct
I33loc = 1/12*ty.*(xy(elmloc(:,2),1)-xy(elmloc(:,1),1)).^3; % correct
I23loc = 1/12*tx.*(xy(elmloc(:,2),1)-xy(elmloc(:,1),1)).*(xy(elmloc(:,2),2)-xy(elmloc(:,1),2)).^2; % correct

loc23 = (xy(elmloc(:,2),:)+xy(elmloc(:,1),:))/2;

if ders == 1
    partdI22loc = 1/12*dtx.*(xy(elmloc(:,2),2)-xy(elmloc(:,1),2)).^3; % correct
    partdI33loc = 1/12*dty.*(xy(elmloc(:,2),1)-xy(elmloc(:,1),1)).^3; % correct
    partdI23loc = 1/12*dtx.*(xy(elmloc(:,2),1)-xy(elmloc(:,1),1)).*(xy(elmloc(:,2),2)-xy(elmloc(:,1),2)).^2;  % correct
end

mA   = sum(rho.*Aloc);
mQ2  = sum(rho.*loc23(:,2).*Aloc);
mQ3  = sum(rho.*loc23(:,1).*Aloc);
mI22 = sum(rho.*I22loc+rho.*Aloc.*(loc23(:,2)).^2);
mI33 = sum(rho.*I33loc+rho.*Aloc.*(loc23(:,1)).^2);
mI23 = sum(rho.*I23loc+rho.*Aloc.*(loc23(:,1)).*(loc23(:,2)));
mI11 = mI22+mI33;

if ders == 1
    dmQ2dt = loc23(:,2).*dmAdt;
    dmQ3dt = loc23(:,1).*dmAdt;
    dmI22dt = rho.*partdI22loc+dmAdt.*(loc23(:,2)).^2;
    dmI33dt = rho.*partdI33loc+dmAdt.*(loc23(:,1)).^2;
    dmI23dt = rho.*partdI23loc+dmAdt.*(loc23(:,1)).*(loc23(:,2));
    dmI11dt = dmI22dt+dmI33dt;
else
    dmAdt   = [];
    dmQ2dt  = [];
    dmQ3dt  = [];
    dmI22dt = [];
    dmI33dt = [];
    dmI23dt = [];
    dmI11dt = [];
end