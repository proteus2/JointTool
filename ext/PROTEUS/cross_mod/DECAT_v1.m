function [S6,Gamma,DS6_mat,DGamma_mat,DS6_geom_y,DS6_geom_z,DGamma_geom_y,DGamma_geom_z] = DECAT_v1(Elm,Node,Prop,flag_strain,flag_sens,flag_geom)
% FACAS: Fast Asymtotically correct Cross-sectional Analysis Software
%this function calculates the cross-sectional stiffness matrix
%input:Elm matric containing the beam elementgs
%      elm.y: [y1 y2] the y coordinate of the element nodes
%      elm.z: [z1 z2] the z coordinate of the element nodes
%      elm.dofi: dof of node i
%      elm.dofj: dof of node j
%output S6,Gamma,DS6,DGamma_eps(the thimoshenko stiffness matrix and strain distribution across the cross-section, and sensitivity results if requested)

%% constants
E = [0 0 -1 0;0 1 0 0];
Ie = zeros(4,6);
Ie(1,1)= 1;Ie(2,5)= 1;Ie(3,6)= 1;Ie(4,4)= 1;
Is = zeros(2,6);
Is(1,2)= 1;Is(2,3)= 1;
ndof = 4; %Dof per node
Cons_elm = 1; % id constrained node to remove rigid body modes
%% assemble: F,H0,K00,H1,K10,R,B0,G
Ne      = size(Elm,1); %number of elements
dof_max = ndof*max(max(Elm(:,3:4))); %maximum dof

[F,H0,K00,H1,K10,B0,G] = assemble(Ne,dof_max,Elm,Node,Prop,Ne); %note K00 is singular so 1 node needs to be fixed

%Re = R(:,[1,5,6,4]); % rigid body modes related to the euler forces
%Rn = R(:,[1,2,3,4]); % the rigid body modes that span the null space of K00


%% contrained dof and free dof
dof_tot       = 1:dof_max; %total dof
node_clamped  = Elm(Cons_elm,3); %constrained node
dof_clamped   = (node_clamped-1)*ndof+1:ndof*node_clamped;%constrained dof
dof_free      = setdiff(dof_tot,dof_clamped); % free dof

%% first order approximation of warping w0 and Euler stiffness matrix
V0 = zeros(size(H0));
V0(dof_free,:) = -K00(dof_free,dof_free)\H0(dof_free,:);

Se = F+H0'*V0; % Euler: modulus/stiffness matrix
Se = (Se+Se')/2; %Symmetrize

%% second order approximation of warping w1 and the Timoshenko stiffness matrix
V1 = zeros(size(H0));
H1b = H1+K10*V0;
V1(dof_free,:) = K00(dof_free,dof_free)\(H1b(dof_free,:));

P=H0'*V1;
Ce =  Se\eye(4);                     %Euler: compliance/flexibility matrix
Cs =  E*Ce*(V1'*H1b+P'*Ce*P)*Ce*E';  %Shear: compliance of the shear stiffness components
Ces = Ce*P*Ce*E';                    %Coupling: euler-shear forces 

C6 = Ie'*Ce*Ie+Is'*Cs*Is-(Ie'*Ces*Is+Is'*Ces'*Ie); % Timoshenko: compliance/flexibility matrix

S6 = C6\eye(6);  % Timoshenko: modulus/stiffness matrix
S6 = (S6+S6')/2; %Symmetrize 

%% Normalized strain variation over cross-section: Gamma = Gamma_hat*epsilon

if(strcmp(flag_strain,'yes'))%check if strain analysis is requested
    Gamma_euler = (G+B0*V0)*Ce;
    Gamma_shear = (B0*V1-Gamma_euler*P)*Ce;
    Gamma       = Gamma_euler*Ie+Gamma_shear*E'*Is; %strain across cross-seaction
else
    Gamma_euler = [];
    Gamma       = [];
end

%% Sensitivity analysis
if(strcmp(flag_sens,'yes'))%check if sensitivity analysis is requested
    DS6_mat     = zeros(21,18*size(Prop,2));
    DS6_geom_y  = zeros(21,size(Node,1));
    DS6_geom_z  = zeros(21,size(Node,1));

    DGamma_mat      = zeros(6*Ne,6*18*size(Prop,2));
    DGamma_geom_y   = zeros(6*Ne,6*size(Node,1));
    DGamma_geom_z   = zeros(6*Ne,6*size(Node,1));
  
    %% adjoint sensitivity papameters
    RV101            = K10'*V1;
    V101             = zeros(size(RV101));
    V101(dof_free,:) = K00(dof_free,dof_free)\RV101(dof_free,:);
    
    RV010            = K10'*V0;
    V010             = zeros(size(RV010));
    V010(dof_free,:) = K00(dof_free,dof_free)\RV010(dof_free,:);
    
    VH0              = zeros(size(H0));
    VH0(dof_free,:)  = K00(dof_free,dof_free)\H0(dof_free,:);
    
    B0t              = B0';
    VB0              = zeros(size(B0t));
    VB0(dof_free,:)  = K00(dof_free,dof_free)\B0t(dof_free,:);
    
    K10t             = K10';
    VK10             = zeros(size(K10t));
    VK10(dof_free,:) = K00(dof_free,dof_free)\K10t(dof_free,:);
    
    %% derivative wrt to laminate stiffness Cij
    DProp = zeros(size(Prop));
    for el=1:size(Elm,1) %loop over the active elements
        for cij=1:18 % loop over the elements of C
            switch cij % set up the dProp matrix
                case {1,2,3,4,5,6} % dA
                    dB = zeros(3,3);dD = zeros(3,3);
                    dA_vec = zeros(6); dA_vec(cij) = 1;
                    dA = symmat2vec(dA_vec,-1,3);%vector to matrix
                    dC = [dA,dB;dB,dD];
                    DProp(:,Elm(el,2)) = symmat2vec(dC,1);
                case{7,8,9,10,11,12} %dB
                    dA = zeros(3,3);dD = zeros(3,3);
                    dB_vec = zeros(6); dB_vec(cij-6) = 1;
                    dB = symmat2vec(dB_vec,-1,3);%vector to matrix
                    dC = [dA,dB;dB,dD];
                    DProp(:,Elm(el,2)) = symmat2vec(dC,1);
                case{13,14,15,16,17,18} %dD
                    dA = zeros(3,3);dB = zeros(3,3);
                    dD_vec = zeros(6); dD_vec(cij-2*6) = 1;
                    dD = symmat2vec(dD_vec,-1,3);%vector to matrix
                    dC = [dA,dB;dB,dD];
                    DProp(:,Elm(el,2)) = symmat2vec(dC,1);
            end
            [DS6_elm,DGamma_elm] =  sensitivity(dof_max,Elm(el,:),Node,DProp,V0,V1,V101,V010,H1b,P,Ce,S6,E,Ie,Is,'mat',[],[],Ne,VH0,VB0,VK10,B0,G,Gamma_euler,flag_strain);
            DS6_mat(:,18*(Elm(el,2)-1)+cij) = DS6_mat(:,18*(Elm(el,2)-1)+cij)+symmat2vec(DS6_elm,1); %matrix to vector
            Nselect = 6*18*(Elm(el,2)-1)+6*(cij-1)+(1:6);
            if(strcmp(flag_strain,'yes'))%check if strain analysis is requested
                DGamma_mat(:,Nselect) = DGamma_mat(:,Nselect)+DGamma_elm;
            end
        end
    end
    if(strcmp(flag_geom,'yes'))%check if sensitivity analysis wrt geometry is requested
        %% derivative wrt to goemetry; i.e node locations
        Nnode = size(Node,1); % # of nodes
        for nd = 1:Nnode %loop over the nodes
            Nselect = 6*(nd-1)+1:6*nd;
            dNode = zeros(size(Node)); %initialization of geometric derivative
            [row,col] = find(Elm(:,[3,4]) == nd);  % find the elements where the node(nd) operates on
            %% derivative wrt y coordinate of Node: nd
            dNode(nd,1) = 1;dNode(nd,2) = 0;
            [DS6_node,DGamma_node] =  sensitivity(dof_max,Elm,Node,Prop,V0,V1,V101,V010,H1b,P,Ce,S6,E,Ie,Is,'geom',flipud(row),dNode,Ne,VH0,VB0,VK10,B0,G,Gamma_euler);
            DS6_geom_y(:,nd) =  symmat2vec(DS6_node,1);
            DGamma_geom_y(:,Nselect) = DGamma_node;
            %% derivative wrt z coordinate of Node: nd
            dNode(nd,1) = 0;dNode(nd,2) = 1;
            [DS6_node,DGamma_node] =  sensitivity(dof_max,Elm,Node,Prop,V0,V1,V101,V010,H1b,P,Ce,S6,E,Ie,Is,'geom',flipud(row),dNode,Ne,VH0,VB0,VK10,B0,G,Gamma_euler);
            DS6_geom_z(:,nd) =  symmat2vec(DS6_node,1);
            DGamma_geom_z(:,Nselect) = DGamma_node;
        end
    end
else
    DS6_mat     = [];
    DS6_geom_y  = [];
    DS6_geom_z  = [];
    DGamma_mat      = [];
    DGamma_geom_y   = [];
    DGamma_geom_z   = [];
end

end

function [F,H0,K00,H1,K10,B0,G,R,L] = assemble(Ne,dof,Elm,Node,Prop,Ne_tot)

ndof = 4; %dof per node

K00 = zeros(dof,dof);
H0  = zeros(dof,4);
H1  = zeros(dof,4);
K10 = zeros(dof,dof);
F   = zeros(4,4);
R   = zeros(dof,6); % maxtrix where the columns are the 6 ridig body modes of a cross-section

B0  = zeros(6*Ne_tot,dof); %B0 matrix of all elements
G   = zeros(6*Ne_tot,4); %G matrix for all elements
L   = 0;

for el = 1:Ne

 elm_node = Elm(el,3:4); %local nodes for element:el
 elmdof = [(elm_node(1)-1)*ndof+1:ndof*elm_node(1),(elm_node(2)-1)*ndof+1:ndof*elm_node(2)];%local dof of element:el
 nelm = 6*(el-1)+1:6*el; 
 %run the function that calculates local properties
 [Fe,H0e,K00e,H1e,K10e,Re,B0e,Ge,Le] = elmvar_reduced(Elm(el,:),Node,Prop(:,Elm(el,2)));
 
 %assemble H,K,etc.....
 H0(elmdof,:)       = H0(elmdof,:)+H0e;
 K00(elmdof,elmdof) = K00(elmdof,elmdof)+K00e;
 H1(elmdof,:)       = H1(elmdof,:)+H1e;
 K10(elmdof,elmdof) = K10(elmdof,elmdof)+K10e;
 F                  = F+Fe;
 L                  = L+Le;
 %strain comp
 R(elmdof,:)        = R(elmdof,:)+Re;
 B0(nelm,elmdof)    = B0e;
 G(nelm,:)          = Ge;
end

end

function [F,H0,K00,H1,K10,R,B0,G,L] = elmvar_reduced(Elm,Node,Prop)
% this function calulates the coefficients of element strain energy

%% element properties
C = symmat2vec(Prop,-1,6);      %material property of an element
y = Node([Elm(3) Elm(4)],1);    %y coordinate of the 2 nodes
z = Node([Elm(3) Elm(4)],2);    %z coordinate of the 2 nodes

L = sqrt(diff(y)^2+diff(z)^2);
ydot = diff(y)/L; %y derivative wrt s
zdot = diff(z)/L; %z derivative wrt s

%% the rigid body modes of currect element: ri i1,2,..6 for 8dof
  r1 = [1,0,0,0,1,0,0,0]';
  r2 = [0,ydot,-zdot,0,0,ydot,-zdot,0]';
  r3 = [0,zdot,ydot,0,0,zdot,ydot,0]';
  r4 = [0,y(1)*zdot-z(1)*ydot,z(1)*zdot+y(1)*ydot,-2,0,y(2)*zdot-z(2)*ydot,z(2)*zdot+y(2)*ydot,-2]';
  r5 = [z(1),0,0,0,z(2),0,0,0]';
  r6 = -[y(1),0,0,0,y(2),0,0,0]';
  
%% B0 and B1 and G
  B0i = [0 0 0 0 0 0 0 0;
         0 -1/L 0 0 0 1/L 0 0;
         -1/L 0 0 0 1/L 0 0 0;
         0 0 0 0 0 0 0 0;
         0 0 -6/L^2 4/L 0 0 6/L^2 2/L;
         0 0 0 0 0 0 0 0];
     
  B0j = [0 0 0 0 0 0 0 0;
         0 -1/L 0 0 0 1/L 0 0;
         -1/L 0 0 0 1/L 0 0 0;
         0 0 0 0 0 0 0 0;
         0 0 6/L^2 -2/L 0 0 -6/L^2 -4/L;
         0 0 0 0 0 0 0 0];
     
 B1i = [1 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 -1/L -1/2 0 0 1/L 1/2];
 
 B1j = [0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 -1/L 1/2 0 0 1/L -1/2];


%% normal notation of the rigid body modes    
Gi = [1 z(1)  -y(1)      0;
      0 0     0         0;
      0 0     0 1*(y(1)*zdot-z(1)*ydot);
      0 ydot  zdot     0;
      0 0     0         0;
      0 0     0        -2];
Gj = [1 z(2)  -y(2)      0;
      0 0     0         0;
      0 0     0 1*(y(2)*zdot-z(2)*ydot);
      0 ydot  zdot     0;
      0 0     0         0;
      0 0     0        -2];

%% transformation matrix from local to global
T = eye(8);
T(2,2) = ydot;
T(2,3) = zdot;
T(3,2) = -zdot;
T(3,3) = ydot;
T(6,6) = ydot;
T(6,7) = zdot;
T(7,6) = -zdot;
T(7,7) = ydot;


%% coefficients of the strain energy expression
F = L*(1/3*Gi'*C*Gi+1/6*(Gi'*C*Gj+Gj'*C*Gi)+1/3*Gj'*C*Gj);
H0 = L*(1/3*B0i'*C*Gi+1/6*B0i'*C*Gj+1/6*B0j'*C*Gi+1/3*B0j'*C*Gj);
K00 =L*(1/3*B0i'*C*B0i+1/6*(B0i'*C*B0j+B0j'*C*B0i)+1/3*B0j'*C*B0j);
H1 = L*(1/3*B1i'*C*Gi+1/6*B1i'*C*Gj+1/6*B1j'*C*Gi+1/3*B1j'*C*Gj);
K10 =L*(1/3*B1i'*C*B0i+1/6*(B1i'*C*B0j+B1j'*C*B0i)+1/3*B1j'*C*B0j);
%K11 =L*(1/3*B1i'*C*B1i+1/6*(B1i'*C*B1j+B1j'*C*B1i)+1/3*B1j'*C*B1j);
  
%% transformation to the global(cartesian) coordinate system
H0 = T*H0;
K00 = T*K00*T';
H1 = T*H1;
K10 = T*K10*T';

%% some additional parameters
R = T*[r1 r2 r3 r4 r5 r6]; %matrix of rigid body modes(in global coordinate system)
B0 = 1/2*(B0i+B0j)*T'; %geometry info strain contribution of warping displacement
G = 1/2*(Gi+Gj); %strain contribution of the beam strains
    
end


function[DS6,DGamma] = sensitivity(dof_max,Elm,Node,Prop,V0,V1,V101,V010,H1b,P,Ce,S6,E,Ie,Is,flag,Nelm,dNode,Ne,VH0,VB0,VK10,B0,G,Ge,flag_strain)

    switch flag
        case {'mat'}
            [DF,DH0,DK00,DH1,DK10,DB0,DG] = assemble(1,dof_max,Elm,Node,Prop,Ne); %derivative of the element cefficients
            DB0 = 0*DB0;DG = 0*DG; %B0 and G are not dependent on stiffness components/material
        case {'geom'}
            [DF,DH0,DK00,DH1,DK10,DB0,DG] = sens_elm_geom(Node,dNode,Elm,Prop,dof_max,Nelm,Ne);
    end
    
    %% derivative of the first order approximation of warping w0 and Euler stiffness matrix
    DSe = DF+V0'*DK00*V0+DH0'*V0+V0'*DH0; % Euler: modulus/stiffness matrix
    DSe = (DSe+DSe')/2; %Symmetrize
                
    %% derivative of the Timoshenko stiffness matrix
    DV1H1b  =V1'*(DH1+DK10*V0)-V101'*(DH0+DK00*V0);
    DP      =DH0'*V1+V0'*DK00*V1+V010'*(DH0+DK00*V0)-V0'*(DH1+DK10*V0); 
    
    DCe = -Ce*DSe*Ce;                                          %Euler: compliance/flexibility matrix
    DCs =  E*DCe*(V1'*H1b+P'*Ce*P)*Ce*E'+...
           E*Ce*(DV1H1b'+DV1H1b-V1'*DK00*V1+DP'*Ce*P+P'*DCe*P+P'*Ce*DP)*Ce*E'+...
           E*Ce*(V1'*H1b+P'*Ce*P)*DCe*E' ;                     %Shear: compliance of the shear stiffness components
    DCes = (DCe*P*Ce+Ce*DP*Ce+Ce*P*DCe)*E';                    %Coupling: euler-shear forces 
  
    DC6 = Ie'*DCe*Ie+Is'*DCs*Is-(Ie'*DCes*Is+Is'*DCes'*Ie); % Timoshenko: compliance/flexibility matrix

    DS6 = -S6*DC6*S6;
    DS6 = (DS6+DS6')/2; %Symmetrize
    
    if(strcmp(flag_strain,'yes'))%check if strain analysis is requested
        %% derivative of strain
        
        % Definition of intermediate matrices for speed
        MatA = (DH1+DK10*V0-DK00*V1);
        MatB = VK10'*(DH0+DK00*V0);
        MatC = Ce*DSe*Ce;
        
        DGe     = (DG+DB0*V0)*Ce-VB0'*(DH0+DK00*V0)*Ce-(G+B0*V0)*MatC;
        DP      = DH0'*V1+VH0'*MatA-VH0'*MatB;
        DGs     = DB0*V1*Ce+VB0'*MatA*Ce-VB0'*MatB*Ce-(DGe*P+Ge*DP)*Ce-(B0*V1-Ge*P)*MatC;
        DGamma  = DGe*Ie+DGs*E'*Is; %strain across cross-seaction
    else
        DGamma = [];
    end

end

function[dF,dH0,dK00,dH1,dK10,dB0,dG] = sens_elm_geom(Node,dNode,Elm,Prop,dof,Nelm,Ne_tot)
 
 ndof = 4; %dof per node
 %% initialization
 dK00 = zeros(dof,dof);
 dH0 = zeros(dof,4);
 dH1 = zeros(dof,4);
 dK10 = zeros(dof,dof);
 dF = zeros(4,4);
 
 dB0     = zeros(6*Ne_tot,dof); %B0 matrix of all elements
 dG      = zeros(6*Ne_tot,4); %G matrix for all elements

 %% run over the active elements
 for el=1:length(Nelm) %loop over the elements where node: nd, operates
    elm_node = Elm(Nelm(el),3:4); %local nodes for element:el
    elmdof = [(elm_node(1)-1)*ndof+1:ndof*elm_node(1),(elm_node(2)-1)*ndof+1:ndof*elm_node(2)];%local dof of element:el
    nelm = 6*(Nelm(el)-1)+1:6*Nelm(el); 
    
   [dFe,dH0e,dK00e,dH1e,dK10e,dB0e,dGe] = sens_elmvar_reduced(Elm(Nelm(el),:),Node,dNode,Prop(:,Elm(Nelm(el),2)));
    
   %assemble H,K,etc.....
   dH0(elmdof,:)        = dH0(elmdof,:)+dH0e;
   dK00(elmdof,elmdof)  = dK00(elmdof,elmdof)+dK00e;
   dH1(elmdof,:)        = dH1(elmdof,:)+dH1e;
   dK10(elmdof,elmdof)  = dK10(elmdof,elmdof)+dK10e;
   dF                   = dF+dFe;
   
   %strain comp
   dB0(nelm,elmdof)     = dB0e;
   dG(nelm,:)           = dGe;
 end

end

function [dF,dH0,dK00,dH1,dK10,dB0,dG] = sens_elmvar_reduced(Elm,Node,dNode,Prop)
%% element properties
C = symmat2vec(Prop,-1,6);      %material property of an element
y = Node([Elm(3) Elm(4)],1);    %y coordinate of the 2 nodes
z = Node([Elm(3) Elm(4)],2);    %z coordinate of the 2 nodes

L = sqrt(diff(y)^2+diff(z)^2);
ydot = diff(y)/L; %y derivative wrt s
zdot = diff(z)/L; %z derivative wrt s

%% derivatives
dy = dNode([Elm(3) Elm(4)],1);    %derivative of y coordinate of the 2 nodes
dz = dNode([Elm(3) Elm(4)],2);    %derivative of z coordinate of the 2 nodes

dL = 1/L*(diff(y)*diff(dy)+diff(z)*diff(dz));
dydot = 1/L*(diff(dy)-dL*ydot);
dzdot = 1/L*(diff(dz)-dL*zdot);


%% B0 and B1 and G
  B0i = [0 0 0 0 0 0 0 0;
         0 -1/L 0 0 0 1/L 0 0;
         -1/L 0 0 0 1/L 0 0 0;
         0 0 0 0 0 0 0 0;
         0 0 -6/L^2 4/L 0 0 6/L^2 2/L;
         0 0 0 0 0 0 0 0];
     
  B0j = [0 0 0 0 0 0 0 0;
         0 -1/L 0 0 0 1/L 0 0;
         -1/L 0 0 0 1/L 0 0 0;
         0 0 0 0 0 0 0 0;
         0 0 6/L^2 -2/L 0 0 -6/L^2 -4/L;
         0 0 0 0 0 0 0 0];
     
 B1i = [1 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 -1/L -1/2 0 0 1/L 1/2];
 
 B1j = [0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 -1/L 1/2 0 0 1/L -1/2];

%% derivative of B0 and B1 and G
  dB0i = dL*[0 0 0 0 0 0 0 0;
         0 1/L^2 0 0 0 -1/L^2 0 0;
         1/L^2 0 0 0 -1/L^2 0 0 0;
         0 0 0 0 0 0 0 0;
         0 0 6*2*1/L^3 -4*1/L^2 0 0 -6*2*1/L^3 -2*1/L^2;
         0 0 0 0 0 0 0 0];
     
  dB0j = dL*[0 0 0 0 0 0 0 0;
         0 1/L^2 0 0 0 -1/L^2 0 0;
         1/L^2 0 0 0 -1/L^2 0 0 0;
         0 0 0 0 0 0 0 0;
         0 0 -6*2*1/L^3 2/L^2 0 0 6*2*1/L^3 4*1/L^2;
         0 0 0 0 0 0 0 0];
     
 dB1i = dL*[0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 1/L^2 0 0 0 -1/L^2 0];
 
 dB1j = dL*[0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 1/L^2 0 0 0 -1/L^2 0];    

%% beam strain projected on the shell strain using the rigid body modes    
Gi = [1 z(1)  -y(1)      0;
      0 0     0         0;
      0 0     0 (y(1)*zdot-z(1)*ydot);
      0 ydot  zdot     0;
      0 0     0         0;
      0 0     0        -2];
Gj = [1 z(2)  -y(2)      0;
      0 0     0         0;
      0 0     0 1*(y(2)*zdot-z(2)*ydot);
      0 ydot  zdot     0;
      0 0     0         0;
      0 0     0        -2];
  
%% derivative of the beam strain projected on the shell strain using the rigid body modes    
dGi = [0 dz(1) -dy(1)    0;
       0 0     0         0;
       0 0     0 (dy(1)*zdot-dz(1)*ydot)+(y(1)*dzdot-z(1)*dydot);
       0 dydot dzdot     0;
       0 0     0         0;
       0 0     0         0];
   
dGj = [0 dz(2)  -dy(2)  0;
      0 0     0         0;
      0 0     0 (dy(2)*zdot-dz(2)*ydot)+(y(2)*dzdot-z(2)*dydot);
      0 dydot  dzdot    0;
      0 0     0         0;
      0 0     0         0];  
      
%% transformation matrix from local to global
T = eye(8);
T(2,2) = ydot;
T(2,3) = zdot;
T(3,2) = -zdot;
T(3,3) = ydot;
T(6,6) = ydot;
T(6,7) = zdot;
T(7,6) = -zdot;
T(7,7) = ydot;

%% derivative of transformation matrix from local to global
dT = zeros(8);
dT(2,2) = dydot;
dT(2,3) = dzdot;
dT(3,2) = -dzdot;
dT(3,3) = dydot;
dT(6,6) = dydot;
dT(6,7) = dzdot;
dT(7,6) = -dzdot;
dT(7,7) = dydot;

%% coefficients of the strain energy expression
H0 = L*(1/3*B0i'*C*Gi+1/6*(B0i'*C*Gj+B0j'*C*Gi)+1/3*B0j'*C*Gj);
H1 = L*(1/3*B1i'*C*Gi+1/6*(B1i'*C*Gj+B1j'*C*Gi)+1/3*B1j'*C*Gj);
K00 =L*(1/3*B0i'*C*B0i+1/6*(B0i'*C*B0j+B0j'*C*B0i)+1/3*B0j'*C*B0j);
K10 =L*(1/3*B1i'*C*B0i+1/6*(B1i'*C*B0j+B1j'*C*B0i)+1/3*B1j'*C*B0j);

%% derivative of the coefficients of the strain energy expression
dF = 1/3*L*(dGi'*C*Gi+Gi'*C*dGi+dGj'*C*Gj+Gj'*C*dGj)+1/6*L*(dGi'*C*Gj+Gi'*C*dGj+dGj'*C*Gi+Gj'*C*dGi)+dL*(1/3*(Gi'*C*Gi+Gj'*C*Gj)+1/6*(Gi'*C*Gj+Gj'*C*Gi));
dH0 = 1/3*L*(dB0i'*C*Gi+B0i'*C*dGi+dB0j'*C*Gj+B0j'*C*dGj)+1/6*L*(dB0i'*C*Gj+B0i'*C*dGj+dB0j'*C*Gi+B0j'*C*dGi)+dL*(1/3*B0i'*C*Gi+1/6*B0i'*C*Gj+1/6*B0j'*C*Gi+1/3*B0j'*C*Gj);
dH1 = 1/3*L*(dB1i'*C*Gi+B1i'*C*dGi+dB1j'*C*Gj+B1j'*C*dGj)+1/6*L*(dB1i'*C*Gj+B1i'*C*dGj+dB1j'*C*Gi+B1j'*C*dGi)+dL*(1/3*B1i'*C*Gi+1/6*B1i'*C*Gj+1/6*B1j'*C*Gi+1/3*B1j'*C*Gj);
dK00 =1/3*L*(dB0i'*C*B0i+B0i'*C*dB0i+dB0j'*C*B0j+B0j'*C*dB0j)+1/6*L*(dB0i'*C*B0j+B0i'*C*dB0j+dB0j'*C*B0i+B0j'*C*dB0i)+dL*(1/3*B0i'*C*B0i+1/6*(B0i'*C*B0j+B0j'*C*B0i)+1/3*B0j'*C*B0j);
dK10 =1/3*L*(dB1i'*C*B0i+B1i'*C*dB0i+dB1j'*C*B0j+B1j'*C*dB0j)+1/6*L*(dB1i'*C*B0j+B1i'*C*dB0j+dB1j'*C*B0i+B1j'*C*dB0i)+dL*(1/3*B1i'*C*B0i+1/6*(B1i'*C*B0j+B1j'*C*B0i)+1/3*B1j'*C*B0j);

%% transformation to the global(cartesian) coordinate system
dH0 = dT*H0+T*dH0;
dH1 = dT*H1+T*dH1;
dK00 = dT*K00*T'+T*dK00*T'+T*K00*dT';
dK10 = dT*K10*T'+T*dK10*T'+T*K10*dT';

%% sensitivity of strain coefficients
dB0 = 1/2*(dB0i+dB0j)*T'+1/2*(B0i+B0j)*dT'; %geometry info strain contribution of warping displacement
dG = 1/2*(dGi+dGj); %strain contribution of the beam strains

end

