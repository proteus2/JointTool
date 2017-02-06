function [S6,DS6,Strain,DStrain,DStrain_eps] = DTSM(Elm,Node,Prop,eps,Cons_elm,strainflag,sens)
% FACAS: Fast Asymtotically correct Cross-sectional Analysis Software
%this function calculates the cross-sectional stiffness matrix
%input: Elm: [elm_id prop_id id_node1 id_node2 dof_node1 dof_node2]
%       Node: Node locations [y,z]
%       Prop: Structured array with the ABD-matrix
%       eps: 1D strain measures: [E11;Gamma12;Gamma13;K1;K2;K3]
%       Cons_elm: Element to be constraint in warping solution, can be any
%       element
%       Strainflag: 1, strains; 0, no strains
%       Sens: 1, sensitivities; 0, no sensitivities
%output S6: the Timoshenko stiffness matrix
%       DS6: The sensitivities of the timoshenko stiffness matrix wrt to
%       the cross-sectional properties
%       Each column: [dS11;dS12;....]
%       where the columns are first ordered along prop_id and then
%       according to location in the ABD-matrix:
%       [11,12,13,22,23,33,44,45,46,55,56,66]
%       Strain: Cross-sectional strains per cross-sectional element,
%       ordered according to: [Ex1;Gamma1;K11;K21;Ex2;Gamma2;...]
%       DStrain: Direct sensitivities of the strains wrt to the
%       cross-sectional properties, similar to DS6
%       DStrain_eps: Sensitivities of the strains wrt to a change in 1D
%       strain measures
%% assemble: F,H0,K00,H1,K10,K11
Ne = size(Elm,1); %number of elements
dof_max = max(max(Elm(:,5:end))); %maximum dof

[F,H0,K00,H1,K10,K11,R] = assemble(dof_max,Elm,Node,Prop); %note K00 is singular so 1 node needs to be fixed

Re = R(:,[1,5,6,4]); % rigid body modes related to the euler forces
Rn = R(:,[1,2,3,4]); % the rigid body modes that span the null space of K00

dof_tot = 1:1:dof_max; %total dof
dof_clamped = Elm(Cons_elm,5:8); %fixed dof
dof_free = dof_tot(find(ismember(dof_tot,dof_clamped)==0))'; % free dof

%% first order approximation of warping w0 and Euler stiffness matrix
V0 = zeros(size(H0));
V0(dof_free,:) = -K00(dof_free,dof_free)\H0(dof_free,:);

Se = F + H0'*V0;
Se = (Se+Se')/2; %Symmetrize the euler stiffness


%% second order approximation of warping w1 and the Timoshenko stiffness matrix
V1  = zeros(size(H0));
H1b = H1+K10*V0;
V1(dof_free,:) = K00(dof_free,dof_free)\(H1b(dof_free,:));
D = V1'*H1b;

P=H0'*V1;
Ce = Se\eye(4); %inverse of euler stiffness
E = [0 0 -1 0;0 1 0 0];
Ie = zeros(4,6);
Ie(1,1)= 1;Ie(2,5)= 1;Ie(3,6)= 1;Ie(4,4)= 1;
Is = zeros(2,6);
Is(1,2)= 1;Is(2,3)= 1;

C6 = Ie'*Ce*Ie+Is'*E*Ce*(V1'*H1b+P'*Ce*P)*Ce*E'*Is-(Ie'*Ce*P*Ce*E'*Is+Is'*E*Ce*P'*Ce*Ie);
S6 = C6\eye(6);
S6 = (S6+S6')/2; %Symmetrize the timoshenko stiffness

%% calculation of the strain distribution across the cross section
[Strain_var,Curv_var] = strain_var(Elm,Node,Prop,V0,V1,P,Se,S6,Ie,Is,E,eps);

Strain(1:4:4*Ne,1) = Strain_var(:,1);
Strain(2:4:4*Ne,1) = Strain_var(:,2);
Strain(3:4:4*Ne,1) = Curv_var(:,1);
Strain(4:4:4*Ne,1) = Curv_var(:,2);

%% calculation of sensitivity
if sens==1
    Cgen = generator;
    %DProp = Prop;
    V2 = zeros(size(V0));
    KV0 = K10'*V0;
    V2(dof_free,:) = K00(dof_free,dof_free)\(KV0(dof_free,:));
    V3 = zeros(size(V1));
    KV1 = K10'*V1;
    V3(5:end,:) = -K00(dof_free,dof_free)\(KV1(dof_free,:));
    
    invK00 = K00(dof_free,dof_free)\eye(size(K00(dof_free,dof_free)));
    DGamma_eps = zeros(6);
    for el = 1:Ne
        %calculate average element properties
        [B0_elm_tot{el},G_elm_tot{el}] = elmvar_av(Elm(el,:),Node,Prop);
        elmdof = Elm(el,5:end); %local dof of element:el
        
        V0_elm = V0(elmdof,:);
        V1_elm = V1(elmdof,:);
        
        Ie = zeros(4,6);
        Ie(1,1)= 1;Ie(2,5)= 1;Ie(3,6)= 1;Ie(4,4)= 1;
        
        DGamma_eps = (G_elm_tot{el}+B0_elm_tot{el}*V0_elm)*inv(Se)*Ie*S6+(B0_elm_tot{el}*V1_elm-(G_elm_tot{el}+B0_elm_tot{el}*V0_elm)*Ce*P)*Ce*E'*Is*S6;
        DStrain_eps(4*(el-1)+(1:4),:) = DGamma_eps([1,3,4,6],:);
        
        Amat{el} = (G_elm_tot{el}+B0_elm_tot{el}*V0_elm);
        Bmat{el} = (B0_elm_tot{el}*V1_elm-(G_elm_tot{el}+B0_elm_tot{el}*V0_elm)*Ce*P);
        Cmat{el} = Ce*E'*Is*S6;
    end
    for pr = 1:length(Prop) %run over the properties
        el_ac = find(Elm(:,2)==pr); %active elements with property from prop(pr)
        for cij = 1:12
            DProp{pr} = Cgen{cij}; %derivative of that elm
            DS6_dummy = zeros(6);
            DStrain_dummy = zeros(size(Strain_var));
            DCurv_dummy = zeros(size(Curv_var));
%             
%             DStrain_dummy = zeros(size(Strain_var,1),6);
%             DCurv_dummy = zeros(size(Curv_var,1),6);
            for el=1:size(el_ac,1)
%                 [DS6_elm,DStrain_var,DCurv_var] = elm_sensitivity(dof_max,Elm,Elm(el_ac(el),:),Node,Prop,DProp,P,K10,K00,V0,V1,Ce,E,Ie,Is,S6,D,dof_free,eps,strainflag,Deps{pr,cij}); %sensitivity for one element
                [DS6_elm,DStrain_var_x,DCurv_var_x] = elm_sensitivity(dof_max,Elm,Elm(el_ac(el),:),Node,DProp,P,K10,invK00,V0,V1,V2,V3,Ce,E,Ie,Is,S6,D,B0_elm_tot,G_elm_tot,Amat,Bmat,Cmat,dof_free,eps,strainflag); %sensitivity for one element
                
                DS6_dummy =DS6_dummy+DS6_elm;
%                 if el == 1
%                     DStrain_dummy = DStrain_dummy+DStrain_var_x+DStrain_var_eps;
%                     DCurv_dummy = DCurv_dummy+DCurv_var_x+DCurv_var_eps;
%                 else
                    DStrain_dummy = DStrain_dummy+DStrain_var_x;
                    DCurv_dummy = DCurv_dummy+DCurv_var_x;
%                 end
            end
            DS6(:,(pr-1)*12+cij) = reshape(DS6_dummy',[],1);
            DStrain(1:4:4*Ne,(pr-1)*12+cij) = DStrain_dummy(:,1);
            DStrain(2:4:4*Ne,(pr-1)*12+cij) = DStrain_dummy(:,2);
            DStrain(3:4:4*Ne,(pr-1)*12+cij) = DCurv_dummy(:,1);
            DStrain(4:4:4*Ne,(pr-1)*12+cij) = DCurv_dummy(:,2);
        end
    end
else
    DS6=0;
    DStrain=0;
    DCurv=0;
end

end

function [DS6,DStrain_var_x,DCurv_var_x] = elm_sensitivity(dof_max,Elm,Elm_ac,Node,DProp,P,K10,invK00,V0,V1,V2,V3,Ce,E,Ie,Is,S6,D,B0_elm_tot,G_elm_tot,Amat,Bmat,Cmat,dof_free,eps,strainflag) %sensitivity for one element

%% assemble: DF,DH0,DK00,DH1,DK10,DK11
   [DF,DH0,DK00,DH1,DK10] = assemble(dof_max,Elm_ac,Node,DProp); %note K00 is singular so 1 node needs to be fixed

%% Sensitivity of Euler stiffness matrix  
   DSe = DF + DH0'*V0+V0'*(DH0+DK00*V0);
   DSe = (DSe+DSe')/2; 
   DCe = -Ce*DSe*Ce; % the euler compliance matrix 
   
%% Sensitivity of Timoshenko stiffness matrix
    DP = DH0'*V1+V0'*DK00*V1-V0'*(DH1+DK10*V0)+V2'*(DH0+DK00*V0);
    D1 = V1'*(DH1+DK10*V0)+V3'*(DH0+DK00*V0);
    DD = D1'+D1-V1'*DK00*V1;
         
    DCs = E*DCe'*(D+P'*Ce*P)*Ce*E'+...
               E*Ce'*(DD+DP'*Ce*P+P'*DCe*P+P'*Ce*DP)*Ce*E'+...
               E*Ce'*(D+P'*Ce*P)*DCe*E';
    DCes = -(DCe*P*Ce+Ce*DP*Ce+Ce*P*DCe)*E';
     
    %% derivative of C6 wrt ABD 
     DC6 = Ie'*DCe*Ie + Is'*DCs*Is + Ie'*DCes*Is + Is'*DCes'*Ie;
     DC6 = 1/2*(DC6'+DC6); % ensure the symmetry of the matrix
     %derivative of S6 wrt ABD
     DS6 = -S6*DC6*S6;
     DS6 = (DS6+DS6')/2; %ensure symmetricity
     
     %% derivative of the first and second warping solution (DV0,DV1)
     DV0 = zeros(size(V0));
     DV1 = zeros(size(V0));
     DHV0 = -(DH0+DK00*V0);
     DV0(dof_free,:) = invK00*DHV0(dof_free,:);
     DHV1 = DH1+DK10*V0+K10*DV0-DK00*V1;
     DV1(dof_free,:) = invK00*DHV1(dof_free,:);
    
     %% derivative of strain variation
     if strainflag == 1
        Ne = size(Elm,1); %number of elements
        eps_e = eps([1,5,6,4]); %euler strains
        Dmat = (DCe*P+Ce*DP);
        Emat = (DCe*E'*Is*S6+Ce*E'*Is*DS6);
        for el = 1:Ne
            elmdof = Elm(el,5:end); %local dof of element:el
            
            %calculate average element properties
            B0_elm = B0_elm_tot{el};
            G_elm = G_elm_tot{el};
%             [B0_elm,G_elm] = elmvar_av(Elm(el,:),Node,Prop);
            V0_elm = V0(elmdof,:);
            DV0_elm  = DV0(elmdof,:);
            V1_elm = V1(elmdof,:);
            DV1_elm  = DV1(elmdof,:);
%             Gamma_e = (G_elm+B0_elm*V0_elm)*inv(Se)*Ie*S6*eps;
%             Gamma_s = (B0_elm*V1_elm-(G_elm+B0_elm*V0_elm)*Ce*P)*Ce*E'*Is*S6*eps;
%             Gamma = Gamma_e+Gamma_s;
%             
            Ie = zeros(4,6);
            Ie(1,1)= 1;Ie(2,5)= 1;Ie(3,6)= 1;Ie(4,4)= 1;
            
%             DGamma_eps = (G_elm+B0_elm*V0_elm)*Ie+(B0_elm*V1_elm-(G_elm+B0_elm*V0_elm)*Ce*P)*Ce*E'*Is*S6;

            DGamma_x = B0_elm*DV0_elm*Ce*Ie*S6+Amat{el}*DCe*Ie*S6+Amat{el}*Ce*Ie*DS6+(B0_elm*DV1_elm-B0_elm*DV0_elm*Ce*P-Amat{el}*Dmat)*Cmat{el}+Bmat{el}*Emat;
            
            DGamma_x = DGamma_x*eps;
%             DGamma_eps = DGamma_eps_tot{el}*Deps;
            
            DStrain_var_x(el,:) = [DGamma_x(1),DGamma_x(3)];
            DCurv_var_x(el,:) = [DGamma_x(4),DGamma_x(6)];
            
%             DStrain_var_eps(el,:) = [DGamma_eps(1),DGamma_eps(3)];
%             DCurv_var_eps(el,:) = [DGamma_eps(4),DGamma_eps(6)];
        end
        
        
     else
        DStrain_var_x = 0;
        DCurv_var_x = 0;
        DStrain_var_eps = 0;
        DCurv_var_eps = 0;
     end
end

function[Strain_var,Curv_var] = strain_var(Elm,Node,Prop,V0,V1,P,Se,S6,Ie,Is,E,eps)

Ne = size(Elm,1); %number of elements
eps_e = eps([1,5,6,4]); %euler strains
for el = 1:Ne
 elmdof = Elm(el,5:end); %local dof of element:el
    
 %calculate average element properties
 [B0_elm,G_elm] = elmvar_av(Elm(el,:),Node,Prop);
 V0_elm = V0(elmdof,:);
 V1_elm = V1(elmdof,:);

 Gamma_e = (G_elm+B0_elm*V0_elm)*inv(Se)*Ie*S6*eps;
 Gamma_s = (B0_elm*V1_elm-(G_elm+B0_elm*V0_elm)*inv(Se)*P)*inv(Se)*E'*Is*S6*eps;
 Gamma = Gamma_e+Gamma_s;
 Strain_var(el,:) = [Gamma(1),Gamma(3)];
 Curv_var(el,:) = [Gamma(4),Gamma(6)];
end

end

function [F,H0,K00,H1,K10,K11,R] = assemble(dof,Elm,Node,Prop)

% Ne = size(Elm,1); %number of elements
% dof = max(max(Elm(:,5:end))); %maximum dof

K00 = sparse(dof,dof);
K11 = sparse(dof,dof);
H0 = zeros(dof,4);
H1 = zeros(dof,4);
K10 = sparse(dof,dof);
F = zeros(4,4);
R = zeros(dof,6); % maxtrix where the columns are the 6 ridig body modes of a cross-section


for el = 1:size(Elm,1)
 elmdof = Elm(el,5:end); %local dof of element:el
    
 %run the function that calculates local properties
 [Fe,H0e,K00e,H1e,K10e,K11e,Re] = elmvar(Elm(el,:),Node,Prop);
    
 %assemble H,K,etc.....
 H0(elmdof,:) = H0(elmdof,:)+H0e;
 K00(elmdof,elmdof) = K00(elmdof,elmdof)+K00e;
 K11(elmdof,elmdof) = K11(elmdof,elmdof)+K11e;
 H1(elmdof,:) = H1(elmdof,:)+H1e;
 K10(elmdof,elmdof) = K10(elmdof,elmdof)+K10e;
 R(elmdof,:) = R(elmdof,:)+Re;
 F = F+Fe;
end
K00 = full(K00);
K11 = full(K11);
K10 = full(K10);
end

function [F,H0,K00,H1,K10,K11,R] = elmvar(Elm,Node,Prop)
% this function calulated the coefficients of element strain energy

%% element properties
C = Prop{Elm(2)}; %material property of an element
y = Node([Elm(3) Elm(4)],1); %y coordinate of the 2 nodes
z = Node([Elm(3) Elm(4)],2); %z coordinate of the 2 nodes

L = sqrt((y(2)-y(1))^2+(z(2)-z(1))^2);% length of the element
ydot = (y(2)-y(1))/L; %y derivative wrt s
zdot = (z(2)-z(1))/L; %z derivative wrt s

%% the rigid body modes of currect element: ri i1,2,..6
  % 8dof
  r1 = [1,0,0,0,1,0,0,0]';
  r2 = [0,ydot,-zdot,0,0,ydot,-zdot,0]';
  r3 = [0,zdot,ydot,0,0,zdot,ydot,0]';
  r4 = [0,y(1)*zdot-z(1)*ydot,z(1)*zdot+y(1)*ydot,-1,0,y(2)*zdot-z(2)*ydot,z(2)*zdot+y(2)*ydot,-1]';
  r5 = [z(1),0,0,0,z(2),0,0,0]';
  r6 = -[y(1),0,0,0,y(2),0,0,0]';
%   10dof
%   r1 = [1,0,0,0,0,1,0,0,0,0]';
%   r2 = [0,ydot,-zdot,0,0,0,ydot,-zdot,0,0]';
%   r3 = [0,zdot,ydot,0,0,0,zdot,ydot,0,0]';
%   r4 = [0,y(1)*zdot-z(1)*ydot,z(1)*zdot+y(1)*ydot,-1,0,0,y(2)*zdot-z(2)*ydot,z(2)*zdot+y(2)*ydot,-1,0]';
%   r5 = [z(1),0,0,0,ydot,z(2),0,0,0,ydot]';
%   r6 = -[y(1),0,0,0,zdot,y(2),0,0,0,zdot]';
  
  Re = [r1 r5 r6 r4];
  
  
%% B0 and B1 and G
  B0i = [0 0 0 0 0 0 0 0;
         0 -1/L 0 0 0 1/L 0 0;
         -1/L 0 0 0 1/L 0 0 0;
         0 0 0 0 0 0 0 0;
         0 0 6/L^2 -4/L 0 0 -6/L^2 -2/L;
         0 0 0 0 0 0 0 0];
     
  B0j = [0 0 0 0 0 0 0 0;
         0 -1/L 0 0 0 1/L 0 0;
         -1/L 0 0 0 1/L 0 0 0;
         0 0 0 0 0 0 0 0;
         0 0 -6/L^2 2/L 0 0 6/L^2 4/L;
         0 0 0 0 0 0 0 0];
     
 B1i = [1 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 1/L 1/2 0 0 -1/L -1/2];
 
 B1j = [0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 1/L -1/2 0 0 -1/L 1/2];
 
Gi = B1i*Re;
Gj = B1j*Re;

%% B0,B1,G using 5 dof localy: U = {u,v,w,ts,tx}
% B0i = [0 0 0 0 0 0 0 0 0 0;
%        0 -1/L 0 0 0 0 1/L 0 0 0;
%        -1/L 0 0 0 0 1/L 0 0 0 0;
%        0 0 0 0 0 0 0 0 0 0;
%        0 0 6/L^2 -4/L 0 0 0 -6/L^2 -2/L 0;
%        0 0 0 0 -1/L 0 0 0 0 1/L];
%    
% B0j = [0 0 0 0 0 0 0 0 0 0;
%        0 -1/L 0 0 0 0 1/L 0 0 0;
%        -1/L 0 0 0 0 1/L 0 0 0 0;
%        0 0 0 0 0 0 0 0 0 0;
%        0 0 -6/L^2 2/L 0 0 0 6/L^2 4/L 0;
%        0 0 0 0 -1/L 0 0 0 0 1/L];
% 
% B1i = [1 0 0 0 0 0 0 0 0 0;
%        0 0 0 0 0 0 0 0 0 0;
%        0 1 0 0 0 0 0 0 0 0;
%        0 0 0 0 1 0 0 0 0 0;
%        0 0 0 0 0 0 0 0 0 0;
%        0 0 1/L 1/2 0 0 -1/L -1/2 0 0];
% B1j = [0 0 0 0 0 1 0 0 0 0;
%        0 0 0 0 0 0 0 0 0 0;
%        0 0 0 0 0 0 1 0 0 0;
%        0 0 0 0 0 0 0 0 0 1;
%        0 0 0 0 0 0 0 0 0 0;
%        0 0 1/L -1/2 0 0 -1/L 1/2 1 0];
%    
% Gi = B1i*Re;
% Gj = B1j*Re;   
%checkrbm = [norm(B0i*r4) norm(B0j*r4)]
%% transformation matrix from local to global
Tf = eye(10);
Tf(2,2) = (y(2)-y(1))/L;
Tf(2,3) = -(z(2)-z(1))/L;
Tf(3,2) = (z(2)-z(1))/L;
Tf(3,3) = (y(2)-y(1))/L;
Tf(7,7) = (y(2)-y(1))/L;
Tf(7,8) = -(z(2)-z(1))/L;
Tf(8,7) = (z(2)-z(1))/L;
Tf(8,8) = (y(2)-y(1))/L;
%T = Tf;
T = Tf([1 2 3 5 6 7 8 10],[1 2 3 5 6 7 8 10]); %Reduced matrix containing the rotation of: u,v,w,theta

%% coefficients of the strain energy expression
F = L*(1/3*Gi'*C*Gi+1/6*(Gi'*C*Gj+Gj'*C*Gi)+1/3*Gj'*C*Gj);
H0 = L*(1/3*B0i'*C*Gi+1/6*B0i'*C*Gj+1/6*B0j'*C*Gi+1/3*B0j'*C*Gj);
K00 =L*(1/3*B0i'*C*B0i+1/6*(B0i'*C*B0j+B0j'*C*B0i)+1/3*B0j'*C*B0j);
H1 = L*(1/3*B1i'*C*Gi+1/6*B1i'*C*Gj+1/6*B1j'*C*Gi+1/3*B1j'*C*Gj);
K10 =L*(1/3*B1i'*C*B0i+1/6*(B1i'*C*B0j+B1j'*C*B0i)+1/3*B1j'*C*B0j);
K11 =L*(1/3*B1i'*C*B1i+1/6*(B1i'*C*B1j+B1j'*C*B1i)+1/3*B1j'*C*B1j);
  
%% transformation to the global(cartesian) coordinate system
H0 = T*H0;
K00 = T*K00*T';
H1 = T*H1;
K10 = T*K10*T';
K11 = T*K11*T';
R = T*[r1 r2 r3 r4 r5 r6]; %matrix of rigid body modes(in global coordinate system)
end

function[B0_elm,G_elm] = elmvar_av(Elm,Node,Prop)
% this function calulated the coefficients of element strain energy

%% element properties
C = Prop{Elm(2)}; %material property of an element
y = Node([Elm(3) Elm(4)],1); %y coordinate of the 2 nodes
z = Node([Elm(3) Elm(4)],2); %z coordinate of the 2 nodes

L = sqrt((y(2)-y(1))^2+(z(2)-z(1))^2);% length of the element
ydot = (y(2)-y(1))/L; %y derivative wrt s
zdot = (z(2)-z(1))/L; %z derivative wrt s

%% the rigid body modes of currect element: ri i1,2,..6
  % 8dof
  r1 = [1,0,0,0,1,0,0,0]';
  r2 = [0,ydot,-zdot,0,0,ydot,-zdot,0]';
  r3 = [0,zdot,ydot,0,0,zdot,ydot,0]';
  r4 = [0,y(1)*zdot-z(1)*ydot,z(1)*zdot+y(1)*ydot,-1,0,y(2)*zdot-z(2)*ydot,z(2)*zdot+y(2)*ydot,-1]';
  r5 = [z(1),0,0,0,z(2),0,0,0]';
  r6 = -[y(1),0,0,0,y(2),0,0,0]';
  
  Re = [r1 r5 r6 r4];
  
  
%% B0 and B1 and G
  B0i = [0 0 0 0 0 0 0 0;
         0 -1/L 0 0 0 1/L 0 0;
         -1/L 0 0 0 1/L 0 0 0;
         0 0 0 0 0 0 0 0;
         0 0 6/L^2 -4/L 0 0 -6/L^2 -2/L;
         0 0 0 0 0 0 0 0];
     
  B0j = [0 0 0 0 0 0 0 0;
         0 -1/L 0 0 0 1/L 0 0;
         -1/L 0 0 0 1/L 0 0 0;
         0 0 0 0 0 0 0 0;
         0 0 -6/L^2 2/L 0 0 6/L^2 4/L;
         0 0 0 0 0 0 0 0];
     
 B1i = [1 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 1/L 1/2 0 0 -1/L -1/2];
 
 B1j = [0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 1/L -1/2 0 0 -1/L 1/2];
 
Gi = B1i*Re;
Gj = B1j*Re;


%% transformation matrix from local to global
Tf = eye(10);
Tf(2,2) = (y(2)-y(1))/L;
Tf(2,3) = -(z(2)-z(1))/L;
Tf(3,2) = (z(2)-z(1))/L;
Tf(3,3) = (y(2)-y(1))/L;
Tf(7,7) = (y(2)-y(1))/L;
Tf(7,8) = -(z(2)-z(1))/L;
Tf(8,7) = (z(2)-z(1))/L;
Tf(8,8) = (y(2)-y(1))/L;
%T = Tf;
T = Tf([1 2 3 5 6 7 8 10],[1 2 3 5 6 7 8 10]); %Reduced matrix containing the rotation of: u,v,w,theta

%% average values of G and B0
G_elm = (Gi+Gj)/2;
B0_elm = (B0i+B0j)/2;
%% transformation to the global(cartesian) coordinate system
B0_elm  = B0_elm*T'; 
end

function [Cgen] = generator


for pa = 1:2
    switch pa
        case 1
            for i=1:6
                switch i
                    case 1
                        Cgen_nom = zeros(6);
                        Cgen_nom(1,1) = 1;
                        Cgen{i} = Cgen_nom;
                    case 4
                        Cgen_nom = zeros(6);
                        Cgen_nom(2,2) = 1;
                        Cgen{i} = Cgen_nom;
                    case 6
                        Cgen_nom = zeros(6);
                        Cgen_nom(3,3) = 1;
                        Cgen{i} = Cgen_nom;
                    case 2
                        Cgen_nom = zeros(6);
                        Cgen_nom(1,2) = 1;Cgen_nom(2,1) = 1;
                        Cgen{i} = Cgen_nom;
                    case 3
                        Cgen_nom = zeros(6);
                        Cgen_nom(1,3) = 1;Cgen_nom(3,1) = 1;
                        Cgen{i} = Cgen_nom;
                    case 5
                        Cgen_nom = zeros(6);
                        Cgen_nom(2,3) = 1;Cgen_nom(3,2) = 1;
                        Cgen{i} = Cgen_nom;
                end
            end
        case 2
            for i=1:6
                switch i
                    case 1
                        Cgen_nom = zeros(6);
                        Cgen_nom(4,4) = 1;
                        Cgen{i+6} = Cgen_nom;
                    case 4
                        Cgen_nom = zeros(6);
                        Cgen_nom(5,5) = 1;
                        Cgen{i+6} = Cgen_nom;
                    case 6
                        Cgen_nom = zeros(6);
                        Cgen_nom(6,6) = 1;
                        Cgen{i+6} = Cgen_nom;
                    case 2
                        Cgen_nom = zeros(6);
                        Cgen_nom(4,5) = 1;Cgen_nom(5,4) = 1;
                        Cgen{i+6} = Cgen_nom;
                    case 3
                        Cgen_nom = zeros(6);
                        Cgen_nom(4,6) = 1;Cgen_nom(6,4) = 1;
                        Cgen{i+6} = Cgen_nom;
                    case 5
                        Cgen_nom = zeros(6);
                        Cgen_nom(5,6) = 1;Cgen_nom(6,5) = 1;
                        Cgen{i+6} = Cgen_nom;
                end
           end
    end
end
 
end
