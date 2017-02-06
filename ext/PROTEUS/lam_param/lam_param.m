function [lampar] = lam_param(mat,matID,t,lam,sens)

% lamparvec = lam2lampar(Lam);
[As,Ds,dAdV,dDdV] = dADcalc(mat.E11(matID),mat.E22(matID),mat.nu12(matID),mat.G12(matID),lam,sens);

A = As*t;
D = Ds*t^3/12;

lampar = lampar_constr(lam,sens);

lampar.A = A;
lampar.D = D;
if sens == 1
    lampar.dAddv = dAdV*t;
    lampar.dDddv = dDdV*t^3/12;
    
    dAdtfull = reshape(As',[],1);
    lampar.dAddv(:,9) = [dAdtfull(1,:);dAdtfull(2,:);dAdtfull(3,:);dAdtfull(5,:);dAdtfull(6,:);dAdtfull(9,:)];
    
    dDdtfull = reshape((Ds*t^2/4)',[],1);
    lampar.dDddv(:,9) = [dDdtfull(1,:);dDdtfull(2,:);dDdtfull(3,:);dDdtfull(5,:);dDdtfull(6,:);dDdtfull(9,:)];
end