function [A,D,dA,dD]=dADcalc(E1,E2,nu12,G12,LamParVec,sens)

% Determine the reduced laminate stiffnesses
nu21 = nu12/E1*E2;
Dnom = 1/(1-nu12*nu21);
Q11 = E1*Dnom;
Q22 = E2*Dnom;
Q12 = nu12*E2*Dnom;
Q66 = G12;

% Determine the lamination invariants
U1 = (3*Q11+3*Q22+2*Q12+4*Q66)/8;
U2 = (Q11-Q22)/2;
U3 = (Q11+Q22-2*Q12-4*Q66)/8;
U4 = (Q11+Q22+6*Q12-4*Q66)/8;
U5 = (Q11+Q22-2*Q12+4*Q66)/8;

% Determine the stiffness invariants
C0 = [U1 U4 0;U4 U1 0;0 0 U5];
C1 = [U2 0 0;0 -U2 0;0 0 0];
C2 = [0 0 U2/2;0 0 U2/2;U2/2 U2/2 0];
C3 = [U3 -U3 0;-U3 U3 0;0 0 -U3];
C4 = [0 0 U3;0 0 -U3;U3 -U3 0];

%*****************************
% Element material properties
%*****************************
% -- Upper/lower layer between spars --
A = (C0+LamParVec(1)*C1+LamParVec(2)*C2+LamParVec(3)*C3+LamParVec(4)*C4);
D = (C0+LamParVec(5)*C1+LamParVec(6)*C2+LamParVec(7)*C3+LamParVec(8)*C4);

if sens == 1

    dAfull = zeros(numel(A),numel(LamParVec));
    dAfull(:,1:4) = [reshape(C1',[],1) reshape(C2',[],1) reshape(C3',[],1) reshape(C4',[],1)];
    
    dA = [dAfull(1,:);dAfull(2,:);dAfull(3,:);dAfull(5,:);dAfull(6,:);dAfull(9,:)];
    
    dDfull = zeros(numel(D),numel(LamParVec));
    dDfull(:,5:8) = [reshape(C1',[],1) reshape(C2',[],1) reshape(C3',[],1) reshape(C4',[],1)];
    
    dD = [dDfull(1,:);dDfull(2,:);dDfull(3,:);dDfull(5,:);dDfull(6,:);dDfull(9,:)];
else
    dA = [];
    dD = [];
end