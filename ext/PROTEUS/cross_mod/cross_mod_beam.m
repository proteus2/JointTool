function [crossmod] = cross_mod_beam(constant)

C=zeros(6,6);
C(1,1)=constant.cross.prop.EA;
C(2,2)=constant.cross.prop.kGA;
C(3,3)=constant.cross.prop.kGA;
% C(1,1)=2e7*constant.cross.mat.E*constant.cross.prop.A;
% C(2,2)=2e9*constant.cross.prop.k*constant.cross.mat.G*constant.cross.prop.A;
% C(3,3)=2e9*constant.cross.prop.k*constant.cross.mat.G*constant.cross.prop.A;
C(4,4)=constant.cross.prop.GJ;
C(5,5)=constant.cross.prop.EI22;
C(5,6)=constant.cross.prop.EI23;
C(6,5)=constant.cross.prop.EI23;
C(6,6)=constant.cross.prop.EI33;

for j=1:constant.str.Ns % Loop over the number of structural elements
    crossmod.mI((j-1)*3+(1:3),1:3) = [constant.cross.prop.mI11,0,0;0,constant.cross.prop.mI22,constant.cross.prop.mI23;0,constant.cross.prop.mI23,constant.cross.prop.mI33];
    crossmod.mQ(j,:) = [constant.cross.prop.mQ2,constant.cross.prop.mQ3];
    crossmod.mA(j,1) = constant.cross.prop.mA;
    crossmod.C((j-1)*6+(1:6),1:6)=C;
end
