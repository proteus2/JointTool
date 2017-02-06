function [constant] = prepr_strain(constant)

F11 = 1/constant.mat.Xt(1)/constant.mat.Xc(1);
F22 = 1/constant.mat.Yt(1)/constant.mat.Yc(1);
F66 = 1/constant.mat.S(1)^2;
F12 = -1/2/sqrt(constant.mat.Xt(1)*constant.mat.Xc(1)*constant.mat.Yt(1)*constant.mat.Yc(1));
F1 = 1/constant.mat.Xt(1)-1/constant.mat.Xc(1);
F2 = 1/constant.mat.Yt(1)-1/constant.mat.Yc(1);

nu21 = constant.mat.nu12(1)/constant.mat.E11(1)*constant.mat.E22(1);
Dnom = 1/(1-constant.mat.nu12(1)*nu21);
Q11 = constant.mat.E11(1)*Dnom;
Q22 = constant.mat.E22(1)*Dnom;
Q12 = constant.mat.nu12(1)*constant.mat.E22(1)*Dnom;
Q66 = constant.mat.G12(1);

G11 = Q11^2*F11+Q12^2*F22+2*F12*Q11*Q12;
G22 = Q12^2*F11+Q22^2*F22+2*F12*Q12*Q22;
G12 = Q11*Q12*F11+Q12*Q22*F22+F12*Q12^2+F12*Q11*Q22;
G66 = 4*Q66^2*F66;
G1 = Q11*F1+Q12*F2;
G2 = Q12*F1+Q22*F2;

u1 = G22-G66/2;
u2 = G66/2;
u3 = 2*G12-2*G22+G66;
u4 = G11-2*G12+G22-G66;
u5 = G2;
u6 = G1-G2;

constant.strain.C0(1) = -(1/4)*u6^2/u4-1;
constant.strain.C1(1) = -(1/2)*u3*u6/u4+u5;
constant.strain.C2(1) = -(1/2)*u3*u6/u4+u5;
constant.strain.C11(1) = -(1/4)*u3^2/u4+u2+u1;
constant.strain.C12(1) = u1-(1/4)*u3^2/u4;
constant.strain.C22(1) = -(1/4)*u3^2/u4+u2+u1;