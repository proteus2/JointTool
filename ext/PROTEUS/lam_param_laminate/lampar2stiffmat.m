function A = lampar2stiffmat(mat,v0A,v1A,v2A,v3A,v4A,outflag)
E11  = mat.E11;
E22  = mat.E22;
G12  = mat.G12;
nu12 = mat.nu12;
nu21 = nu12*E22/E11;

% Constitutive matrix components
Q11 = E11/(1-nu12*nu21);
Q22 = E22/(1-nu12*nu21);
Q12 = nu12*E22/(1-nu12*nu21);
Q66 = G12;

% Material invariants
U1=1/8*(3*Q11+3*Q22+2*Q12+4*Q66);
U2=1/2*(Q11-Q22);
U3=1/8*(Q11+Q22-2*Q12-4*Q66);
U4=1/8*(Q11+Q22+6*Q12-4*Q66);
U5=1/8*(Q11+Q22-2*Q12+4*Q66);

% Stiffness matrix
A0 = [U1, U4, 0 ;
      U4, U1, 0 ;
      0,  0,  U5];
A1 = [U2, 0,  0;
      0, -U2, 0; 
      0,  0,  0]*v1A;
A2 = [0,    0,    U2/2; 
      0,    0,    U2/2;
      U2/2, U2/2, 0 ]*v2A;
A3 = [U3, -U3, 0;
     -U3,  U3, 0;
      0,   0,  -U3]*v3A;  % should be - U3 here (was an error)
A4 = [0,   0,  U3;
      0,   0, -U3; 
      U3, -U3, 0 ]*v4A;

  keyboard
if outflag==1
    A = v0A*(A0 + A1 + A2 + A3 + A4);
elseif outflag==0
    A = A0 + A1 + A2 + A3 + A4;
else 
    err('Reference to non-existent outflag. Select right flag when calling lampar2stiffmat.m')
end


