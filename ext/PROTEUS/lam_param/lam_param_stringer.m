function [lampar] = lam_param_stringer(mat,matID,EA,h,lam)

% lamparvec = lam2lampar(Lam);
[As,Ds,~,~] = dADcalc(mat.E11(matID),mat.E22(matID),mat.nu12(matID),mat.G12(matID),lam,0);

a = As\eye(3);

t = a(1,1)*EA/h;

A = As*t;
D = Ds*t^3/12;

lampar.A = A;
lampar.D = D;
lampar.t = t;
