function [ft,Kt,re,Ke,pdC1_Kt,pdC1_ft,pdC2_Kt,pdC2_ft,dp_Kt_s,dp_ft_s,pdC1_re,pdC2_re,dp_re,dC1_Ke,dC2_Ke,varargout] = elem3d(x,p,C1,C2,Ro2,ders,tailflag,varargin)

if isempty(varargin) == 0
    morph = varargin{1};
    morph_flag = 1;
else
    morph_flag = 0;
end

if morph_flag == 1
    if isfield(morph,'shear')
        Rpsi = [cos(morph.shear)   -sin(morph.shear)   0;
            sin(morph.shear)    cos(morph.shear)   0;
            0           0        1 ];
        
        Ro3 = Ro2*Rpsi;
    else
        Ro3 = Ro2;
    end
    if isfield(morph,'fold')
        Rtheta = expon(morph.fold*morph.c0);
        Ro = Rtheta*Ro3;
    else
        Ro = Ro3;
    end
else
    Ro = Ro2;
end

tg1 = p(4:6);
tg2 = p(10:12);

Rg1 = expon(tg1);
Rg2 = expon(tg2);

x21 = x(4:6)-x(1:3);

d21 = p(7:9)-p(1:3);

lo  = sqrt(x21'*x21);
l   = sqrt((x21+d21)'*(x21+d21));
u   = l-lo;

% rigid rotation

e1=(x21+d21)/l;

qb1=Rg1*Ro*[0;1;0];
qb2=Rg2*Ro*[0;1;0];

qb=(qb1+qb2)/2;

e3b = vec(e1,qb);
e3  = e3b/norm(e3b);

e2 = vec(e3,e1);

Rr=[e1 e2 e3];

q=Rr'*qb;
q1=Rr'*qb1;

nu=q(1)/q(2);
nu11=q1(1)/q(2);
nu12=q1(2)/q(2);
nu21=2*nu-nu11;
nu22=2-nu12;

%local rotations

Re1=Rr'*Rg1*Ro;
Re2=Rr'*Rg2*Ro;

tl1 = logar(Re1);
tl2 = logar(Re2);

% local force vector and tangent stiffness matrix

if morph_flag == 1
    if isfield(morph,'twist')
        [kl,fl,re,Ke] = locel(u,tl1,tl2,lo,C1,C2,morph);
    else
        [kl,fl,re,Ke] = locel(u,tl1,tl2,lo,C1,C2);
    end
elseif morph_flag == 0
    [kl,fl,re,Ke] = locel(u,tl1,tl2,lo,C1,C2);
end

% transformation to the new local coordinates

De1 = invTs(tl1);
De2 = invTs(tl2);

fe=[fl(1)
    De1'*fl(2:4)
    De2'*fl(5:7)];

O3=zeros(3,3);
O1=zeros(1,3);

H1=[1   O1  O1
    O1' De1 O3
    O1' O3  De2];

Dh1=dinvTs(tl1,fl(2:4))*De1;
Dh2=dinvTs(tl2,fl(5:7))*De2;

Kh=[0   O1  O1
    O1' Dh1 O3
    O1' O3  Dh2];

ke=H1'*kl*H1+Kh;


% transformation to the global coordinates

r=[-e1;0;0;0;e1;0;0;0];

B=[ r'
    -nu/l*e3' (1-nu12/2)*e1'+nu11/2*e2'  nu/l*e3' 1/2*(-nu22*e1'+nu21*e2')
    -e3'/l e2' e3'/l 0 0 0
    e2'/l e3' -e2'/l 0 0 0
    -nu/l*e3' 1/2*(-nu12*e1'+nu11*e2')  nu/l*e3' (1-nu22/2)*e1'+nu21/2*e2'
    -e3'/l 0 0 0 e3'/l e2'
    e2'/l 0 0 0 -e2'/l e3'];

fg=B'*fe;

I3=eye(3);

A=(I3-e1*e1')/l;

Dr=[A  O3 -A  O3
    O3 O3  O3 O3
    -A  O3  A  O3
    O3 O3  O3 O3];

G=[0   0    nu/l  nu12/2  -nu11/2  0  0  0    -nu/l  nu22/2  -nu21/2  0
    0   0    1/l     0        0     0  0  0    -1/l     0        0     0
    0  -1/l  0       0        0     0  0  1/l   0       0        0     0]';

II=[O3 I3 O3 O3
    O3 O3 O3 I3];

P=II-[G';G'];

F=P'*fe(2:7);

sF=[skew(F(1:3))
    skew(F(4:6))
    skew(F(7:9))
    skew(F(10:12))];

EE=[Rr O3 O3 O3
    O3 Rr O3 O3
    O3 O3 Rr O3
    O3 O3 O3 Rr];

nab=[0
    (nu*(fe(2)+fe(5))+fe(3)+fe(6))/l
    (fe(4)+fe(7))/l];

Kg=B'*ke*B+Dr*fe(1)-EE*sF*G'*EE'+EE*G*nab*r';

% transformation to the new global coordinates

Dg1 = Ts(tg1);
Dg2 = Ts(tg2);

ft=[fg(1:3)
    Dg1'*fg(4:6)
    fg(7:9)
    Dg2'*fg(10:12)];

Dk1 = dTs(tg1,fg(4:6));
Dk2 = dTs(tg2,fg(10:12));

H2=[I3 O3  O3 O3
    O3 Dg1 O3 O3
    O3 O3  I3 O3
    O3 O3  O3 Dg2];

Kt = H2'*Kg*H2;

Kt(4:6,4:6) = Kt(4:6,4:6)+Dk1;

Kt(10:12,10:12) = Kt(10:12,10:12)+Dk2;

if morph_flag == 1
    
    %% Calculation of Mpsi
    if isfield(morph,'shear')
        if morph.energy
            dp = eye(12,12);
            
            dptg1 = dp(4:6,:);
            dptg2 = dp(10:12,:);
            
            [~,dpRg1] = expon(tg1,dptg1);
            [~,dpRg2] = expon(tg2,dptg2);
            
            dpd21 = dp(7:9,:)-dp(1:3,:);
            dpl = 1/2/l*(2*(x21+d21)'*dpd21);
            dpu = dpl;
            
            dpe1 = 1/l^2*(dpd21.*l-(x21(:,ones(12,1))+d21(:,ones(12,1))).*dpl(ones(3,1),:));
        end
        dpsiRpsi   = [-sin(morph.shear)   -cos(morph.shear)   0;
            cos(morph.shear)   -sin(morph.shear)   0;
            0            0       0];
        
        if morph.energy
            d2psiRpsi  = [-cos(morph.shear)    sin(morph.shear)   0;
                -sin(morph.shear)   -cos(morph.shear)   0;
                0            0       0];
        end
        if isfield(morph,'fold')
            dpsiRo = Rtheta*Ro2*dpsiRpsi;
            if morph.energy
                d2psiRo = Rtheta*Ro2*d2psiRpsi;
            end
        else
            dpsiRo = Ro2*dpsiRpsi;
            if morph.energy
                d2psiRo = Ro2*d2psiRpsi;
            end
        end
        
        dpsiqb1 = Rg1*dpsiRo*[0;1;0];
        dpsiqb2 = Rg2*dpsiRo*[0;1;0];
        dpsiqb = (dpsiqb1+dpsiqb2)/2;
        
        if morph.energy
            d2psiqb1 = Rg1*d2psiRo*[0;1;0];
            d2psiqb2 = Rg2*d2psiRo*[0;1;0];
            d2psiqb = (d2psiqb1+d2psiqb2)/2;
            for i=1:size(dp,2)
                dpqb1(:,i)     = dpRg1(:,:,i)*Ro*[0;1;0];
                dpqb2(:,i)     = dpRg2(:,:,i)*Ro*[0;1;0];
                dpsidpqb1(:,i) = dpRg1(:,:,i)*dpsiRo*[0;1;0];
                dpsidpqb2(:,i) = dpRg2(:,:,i)*dpsiRo*[0;1;0];
            end
            dpqb     = (dpqb1+dpqb2)/2;
            dpsidpqb = (dpsidpqb1+dpsidpqb2)/2;
        end
        
        if morph.energy
            [~,dpsie3b,d2psie3b]    = vec(e1,qb,zeros(size(e1)),dpsiqb,zeros(size(e1)),d2psiqb);
            [~,dpe3b,~,dpsidpe3b] = vec(e1,qb,dpe1,dpqb,[],[],zeros(size(e1)),dpsiqb,zeros(size(e1,1),size(dpqb,2)),dpsidpqb);
            dpsine3b = 1/2/norm(e3b)*(2*e3b(1)*dpsie3b(1,:)+2*e3b(2)*dpsie3b(2,:)+2*e3b(3)*dpsie3b(3,:));
            d2psine3b = 1/2/norm(e3b)^2*((2*(dpsie3b(1,:)*dpsie3b(1,:)+e3b(1)*d2psie3b(1,:))+2*(dpsie3b(2,:)*...
                dpsie3b(2,:)+e3b(2)*d2psie3b(2,:))+2*(dpsie3b(3,:)*dpsie3b(3,:)+e3b(3)*d2psie3b(3,:)))*norm(e3b)-...
                (2*e3b(1)*dpsie3b(1,:)+2*e3b(2)*dpsie3b(2,:)+2*e3b(3)*dpsie3b(3,:))*dpsine3b);
            dpsie3 = 1/norm(e3b)^2*(dpsie3b*norm(e3b)-e3b.*dpsine3b);
            d2psie3 = 1/norm(e3b)^4*((d2psie3b*norm(e3b)+dpsie3b*dpsine3b-dpsie3b.*dpsine3b-e3b.*d2psine3b)*norm(e3b)^2-...
                (dpsie3b*norm(e3b)-e3b.*dpsine3b)*2*norm(e3b)*dpsine3b);
            dpne3b        = 1/2/norm(e3b)*(2*e3b(1).*dpe3b(1,:)+2*e3b(2).*dpe3b(2,:)+2*e3b(3).*dpe3b(3,:));
            dpsidpne3b    = 1/2/norm(e3b)^2*((2*(dpsie3b(1).*dpe3b(1,:)+e3b(1).*dpsidpe3b(1,:))+2*(dpsie3b(2).*dpe3b(2,:)+e3b(2).*...
                dpsidpe3b(2,:))+2*(dpsie3b(3).*dpe3b(3,:)+e3b(3).*dpsidpe3b(3,:)))*norm(e3b)-(2*e3b(1).*dpe3b(1,:)+2*e3b(2).*...
                dpe3b(2,:)+2*e3b(3).*dpe3b(3,:)).*dpsine3b);
            dpe3          = 1/norm(e3b)^2*[dpe3b(1,:)*norm(e3b)-e3b(1).*dpne3b;dpe3b(2,:)*norm(e3b)-e3b(2).*dpne3b;...
                dpe3b(3,:)*norm(e3b)-e3b(3).*dpne3b];
            dpsidpe3      = 1/norm(e3b)^4*([dpsidpe3b(1,:)*norm(e3b)+dpe3b(1,:)*dpsine3b-dpsie3b(1).*dpne3b-e3b(1).*...
                dpsidpne3b;dpsidpe3b(2,:)*norm(e3b)+dpe3b(2,:)*dpsine3b-dpsie3b(2).*dpne3b-e3b(2).*dpsidpne3b;...
                dpsidpe3b(3,:)*norm(e3b)+dpe3b(3,:)*dpsine3b-dpsie3b(3).*dpne3b-e3b(3).*dpsidpne3b]*norm(e3b)^2-...
                [dpe3b(1,:)*norm(e3b)-e3b(1).*dpne3b;dpe3b(2,:)*norm(e3b)-e3b(2).*dpne3b;dpe3b(3,:)*norm(e3b)-...
                e3b(3).*dpne3b]*2*norm(e3b)*dpsine3b);
            
            [~,dpsie2,d2psie2] = vec(e3,e1,dpsie3,zeros(size(dpsie3)),d2psie3,zeros(size(dpsie3)));
            [~,dpe2,~,dpsidpe2] = vec(e3,e1,dpe3,dpe1,[],[],dpsie3,zeros(size(e1)),dpsidpe3,zeros(size(e3,1),size(dpqb,2)));
        else
            [~,dpsie3b]    = vec(e1,qb,zeros(size(e1)),dpsiqb);
            dpsine3b = 1/2/norm(e3b)*(2*e3b(1)*dpsie3b(1,:)+2*e3b(2)*dpsie3b(2,:)+2*e3b(3)*dpsie3b(3,:));
            dpsie3 = 1/norm(e3b)^2*(dpsie3b*norm(e3b)-e3b.*dpsine3b);
            [~,dpsie2] = vec(e3,e1,dpsie3,zeros(size(dpsie3)));
        end
        
        dpsiRr = [zeros(size(dpsie2)) dpsie2 dpsie3];
        
        if morph.energy
            d2psiRr = [zeros(size(dpsie3)) d2psie2 d2psie3];
            dpRr(:,1,:) = dpe1;
            dpRr(:,2,:) = dpe2;
            dpRr(:,3,:) = dpe3;
            dpsidpRr(:,1,:) = zeros(size(dpsidpe2));
            dpsidpRr(:,2,:) = dpsidpe2;
            dpsidpRr(:,3,:) = dpsidpe3;
        end
        
        dpsiq = dpsiRr'*q+Rr'*dpsiqb;
        dpsiq1 = dpsiRr'*q1+Rr'*dpsiqb1;
        
        dpsinu = 1/q(2)^2*(dpsiq(1)*q(2)-q(1)*dpsiq(2));
        dpsinu11 = 1/q(2)^2*(dpsiq1(1)*q(2)-q1(1)*dpsiq(2));
        dpsinu12 = 1/q(2)^2*(dpsiq1(2)*q(2)-q1(2)*dpsiq(2));
        dpsinu21 = 2*dpsinu-dpsinu11;
        dpsinu22 = -dpsinu12;
        
        dpsiRe1 = dpsiRr'*Rg1*Ro+Rr'*Rg1*dpsiRo;
        dpsiRe2 = dpsiRr'*Rg2*Ro+Rr'*Rg2*dpsiRo;
        
        if morph.energy
            d2psiRe1 = d2psiRr'*Rg1*Ro+dpsiRr'*Rg1*dpsiRo+Rr'*Rg1*d2psiRo+dpsiRr'*Rg1*dpsiRo;
            d2psiRe2 = d2psiRr'*Rg2*Ro+dpsiRr'*Rg2*dpsiRo+Rr'*Rg2*d2psiRo+dpsiRr'*Rg2*dpsiRo;
            for i=1:size(dp,2)
                dpRe1(:,:,i) = dpRr(:,:,i)'*Rg1*Ro+Rr'*dpRg1(:,:,i)*Ro;
                dpRe2(:,:,i) = dpRr(:,:,i)'*Rg2*Ro+Rr'*dpRg2(:,:,i)*Ro;
                dpsidpRe1(:,:,i) = dpsidpRr(:,:,i)'*Rg1*Ro+dpRr(:,:,i)'*Rg1*dpsiRo+...
                    dpsiRr'*dpRg1(:,:,i)*Ro+Rr'*dpRg1(:,:,i)*dpsiRo;
                dpsidpRe2(:,:,i) = dpsidpRr(:,:,i)'*Rg2*Ro+dpRr(:,:,i)'*Rg2*dpsiRo+...
                    dpsiRr'*dpRg2(:,:,i)*Ro+Rr'*dpRg2(:,:,i)*dpsiRo;
            end
        end
        
        if morph.energy
            [~,dpsitl1,d2psitl1]    = logar(Re1,dpsiRe1,d2psiRe1);
            [~,dpsitl2,d2psitl2]    = logar(Re2,dpsiRe2,d2psiRe2);
            [~,dptl1,~,dpsidptl1] = logar(Re1,dpRe1,[],dpsiRe1,dpsidpRe1);
            [~,dptl2,~,dpsidptl2] = logar(Re2,dpRe2,[],dpsiRe2,dpsidpRe2);
        else
            [~,dpsitl1]    = logar(Re1,dpsiRe1);
            [~,dpsitl2]    = logar(Re2,dpsiRe2);
        end
        
        morph.dpsitl1 = dpsitl1;
        morph.dpsitl2 = dpsitl2;
        
        if morph.energy
            morph.d2psitl1  = d2psitl1;
            morph.d2psitl2  = d2psitl2;
            morph.dptl1     = dptl1;
            morph.dptl2     = dptl2;
            morph.dpu       = dpu;
            morph.dpsidptl1 = dpsidptl1;
            morph.dpsidptl2 = dpsidptl2;
        end
        
        [~,~,~,~,morphosh] = locel(u,tl1,tl2,lo,C1,C2,morph);
        morph.dpsire = morphosh.dredpsi;
        dpsifl = morphosh.dfldpsi;
        
        if morph.energy
            GAskin = 0;
            morph.Mpsi = fl'*morphosh.dqedpsi+GAskin*morph.shear;
            morph.dpMpsi   = morphosh.dqedpsi'*morphosh.dpfl+fl'*morphosh.dqedpsidp;
            morph.dpsiMpsi = morphosh.dfldpsi'*morphosh.dqedpsi+fl'*morphosh.d2qedpsi+GAskin;
        end
        
        [~,dpsiDe1] = invTs(tl1,dpsitl1);
        [~,dpsiDe2] = invTs(tl2,dpsitl2);
        %
        dpsife = [dpsifl(1)
            dpsiDe1'*fl(2:4)+De1'*dpsifl(2:4)
            dpsiDe2'*fl(5:7)+De2'*dpsifl(5:7)];
        
        dpsiH1 = [0   O1     O1
            O1' dpsiDe1 O3
            O1' O3     dpsiDe2];
        
        [Dh1_dum,dpsiDh1_dum] = dinvTs(tl1,fl(2:4),dpsitl1,dpsifl(2:4));
        [Dh2_dum,dpsiDh2_dum] = dinvTs(tl2,fl(5:7),dpsitl2,dpsifl(5:7));
        dpsiDh1 = dpsiDh1_dum*De1+Dh1_dum*dpsiDe1;
        dpsiDh2 = dpsiDh2_dum*De2+Dh2_dum*dpsiDe2;
        
        dpsiKh = [0   O1     O1
            O1' dpsiDh1 O3
            O1' O3     dpsiDh2];
        
        dpsike = dpsiH1'*kl*H1+H1'*kl*dpsiH1+dpsiKh;
        
        % transformation to the global coordinates
        
        dpsiB = [ zeros(size(r'))
            -1/l*(dpsinu*e3'+nu*dpsie3') (-dpsinu12/2)*e1'+1/2*(dpsinu11*e2'+nu11*dpsie2')  1/l*(dpsinu*e3'+nu*dpsie3') 1/2*(-dpsinu22*e1'+dpsinu21*e2'+nu21*dpsie2')
            -dpsie3'/l dpsie2' dpsie3'/l 0 0 0
            dpsie2'/l dpsie3' -dpsie2'/l 0 0 0
            -1/l*(dpsinu*e3'+nu*dpsie3') 1/2*(-dpsinu12*e1'+dpsinu11*e2'+nu11*dpsie2')  1/l*(dpsinu*e3'+nu*dpsie3') (-dpsinu22/2)*e1'+1/2*(dpsinu21*e2'+nu21*dpsie2')
            -dpsie3'/l 0 0 0 dpsie3'/l dpsie2'
            dpsie2'/l 0 0 0 -dpsie2'/l dpsie3'];
        
        dpsifg = dpsiB'*fe+B'*dpsife;
        
        dpsiG = [0   0           dpsinu/l  dpsinu12/2  -dpsinu11/2  0  0       0        -dpsinu/l  dpsinu22/2  -dpsinu21/2  0
            0   0           0                   0           0     0  0       0         0                        0        0     0
            0  0         0                       0           0     0  0  0           0                          0        0     0]';
        
        dpsiP = -[dpsiG';dpsiG'];
        
        dpsiF = dpsiP'*fe(2:7)+P'*dpsife(2:7);
        
        [~,dpsisF(1:3,:)]   = skew(F(1:3),dpsiF(1:3));
        [~,dpsisF(4:6,:)]   = skew(F(4:6),dpsiF(4:6));
        [~,dpsisF(7:9,:)]   = skew(F(7:9),dpsiF(7:9));
        [~,dpsisF(10:12,:)] = skew(F(10:12),dpsiF(10:12));
        
        dpsiEE = [dpsiRr O3    O3    O3
            O3    dpsiRr O3    O3
            O3    O3    dpsiRr O3
            O3    O3    O3    dpsiRr];
        
        dpsinab = [0
            1/l*(dpsinu*(fe(2)+fe(5))+nu*(dpsife(2)+dpsife(5))+dpsife(3)+dpsife(6))
            1/l*(dpsife(4)+dpsife(7))];
        
        dpsiKg = dpsiB'*ke*B+B'*dpsike*B+B'*ke*dpsiB+Dr*dpsife(1)-dpsiEE*sF*G'*EE'-EE*dpsisF*G'*EE'-EE*sF*dpsiG'*EE'-EE*sF*G'*dpsiEE'+dpsiEE*G*nab*r'+EE*dpsiG*nab*r'+EE*G*dpsinab*r';
        
        % transformation to the new global coordinates
        
        morph.dpsift = [dpsifg(1:3)
            Dg1'*dpsifg(4:6)
            dpsifg(7:9)
            Dg2'*dpsifg(10:12)];
        
        [~,dpsiDk1] = dTs(tg1,fg(4:6),zeros(3,1),dpsifg(4:6));
        [~,dpsiDk2] = dTs(tg2,fg(10:12),zeros(3,1),dpsifg(10:12));
        
        dpsiKt = H2'*dpsiKg*H2;
        
        dpsiKt(4:6,4:6) = dpsiKt(4:6,4:6)+dpsiDk1;
        
        dpsiKt(10:12,10:12) = dpsiKt(10:12,10:12)+dpsiDk2;
        
        morph.dpsiKt = reshape(dpsiKt',[],1);
    end
    %% Calculation of Mphi
    if isfield(morph,'twist')
        % Note Mphi is dependent upon the local force vector that is
        % directly dependent upon the structural deformations
        if morph.energy
            dp = eye(12,12);
            
            dptg1 = dp(4:6,:);
            dptg2 = dp(10:12,:);
            
            [~,dpRg1] = expon(tg1,dptg1);
            [~,dpRg2] = expon(tg2,dptg2);
            
            dpd21 = dp(7:9,:)-dp(1:3,:);
            dpl = 1/2/l*(2*(x21+d21)'*dpd21);
            dpu = dpl;
            
            dpe1 = 1/l^2*(dpd21.*l-[(x21+d21) (x21+d21) ...
                (x21+d21) (x21+d21) (x21+d21) (x21+d21) ...
                (x21+d21) (x21+d21) (x21+d21) (x21+d21) ...
                (x21+d21) (x21+d21)].*[dpl;dpl;dpl]);
            
            for i=1:size(dp,2)
                dpqb1(:,i)     = dpRg1(:,:,i)*Ro*[0;1;0];
                dpqb2(:,i)     = dpRg2(:,:,i)*Ro*[0;1;0];
            end
            dpqb     = (dpqb1+dpqb2)/2;
            
            [~,dpe3b] = vec(e1,qb,dpe1,dpqb);
            dpne3b        = 1/2/norm(e3b)*(2*e3b(1).*dpe3b(1,:)+2*e3b(2).*dpe3b(2,:)+2*e3b(3).*dpe3b(3,:));
            dpe3          = 1/norm(e3b)^2*[dpe3b(1,:)*norm(e3b)-e3b(1).*dpne3b;dpe3b(2,:)*norm(e3b)-e3b(2).*dpne3b;...
                dpe3b(3,:)*norm(e3b)-e3b(3).*dpne3b];
            [~,dpe2] = vec(e3,e1,dpe3,dpe1);
            
            dpRr(:,1,:) = dpe1;
            dpRr(:,2,:) = dpe2;
            dpRr(:,3,:) = dpe3;
            
            for i=1:size(dp,2)
                dpRe1(:,:,i) = dpRr(:,:,i)'*Rg1*Ro+Rr'*dpRg1(:,:,i)*Ro;
                dpRe2(:,:,i) = dpRr(:,:,i)'*Rg2*Ro+Rr'*dpRg2(:,:,i)*Ro;
            end
            
            [~,dptl1] = logar(Re1,dpRe1);
            [~,dptl2] = logar(Re2,dpRe2);
        end
        
        morphi = morph;
        
        if morph.energy
            morphi.dptl1 = dptl1;
            morphi.dptl2 = dptl2;
            morphi.dpu   = dpu;
        end
        
        [~,~,~,~,morphotw] = locel(u,tl1,tl2,lo,C1,C2,morphi);
        
        % In case of GJ for morphing energy of twist
        if morph.energy
            GJphi = 0;
            morph.Mphi     = fl'*morphotw.dqedphi+GJphi*morph.twist;
            morph.dpMphi   = morphotw.dpfl'*morphotw.dqedphi;
            morph.dphiMphi = morphotw.dfldphi'*morphotw.dqedphi+GJphi;
            if isfield(morph,'shear')
                morph.dpsiMphi = morphosh.dfldpsi'*morphotw.dqedphi;
                morph.dphiMpsi = morphotw.dfldphi'*morphosh.dqedpsi;
            end
        end
        
        morph.dphire = morphotw.dredphi;
        dphifl = morphotw.dfldphi;
        
        dphife = [dphifl(1)
            De1'*dphifl(2:4)
            De2'*dphifl(5:7)];
        
        [~,dphiDh1_dum] = dinvTs(tl1,fl(2:4),zeros(3,1),dphifl(2:4));
        [~,dphiDh2_dum] = dinvTs(tl2,fl(5:7),zeros(3,1),dphifl(5:7));
        dphiDh1 = dphiDh1_dum*De1;
        dphiDh2 = dphiDh2_dum*De2;
        
        dphiKh = [0   O1     O1
            O1' dphiDh1 O3
            O1' O3     dphiDh2];
        
        dphike = dphiKh;
        
        dphifg = B'*dphife;
        
        dphiF = P'*dphife(2:7);
        
        [~,dphisF(1:3,:)]   = skew(F(1:3),dphiF(1:3));
        [~,dphisF(4:6,:)]   = skew(F(4:6),dphiF(4:6));
        [~,dphisF(7:9,:)]   = skew(F(7:9),dphiF(7:9));
        [~,dphisF(10:12,:)] = skew(F(10:12),dphiF(10:12));
        
        dphinab = [0
            1/l*(nu*(dphife(2)+dphife(5))+dphife(3)+dphife(6))
            1/l*(dphife(4)+dphife(7))];
        
        dphiKg = B'*dphike*B+Dr*dphife(1)-EE*dphisF*G'*EE'+EE*G*dphinab*r';
        
        % transformation to the new global coordinates
        
        morph.dphift = [dphifg(1:3)
            Dg1'*dphifg(4:6)
            dphifg(7:9)
            Dg2'*dphifg(10:12)];
        
        [~,dphiDk1] = dTs(tg1,fg(4:6),zeros(3,1),dphifg(4:6));
        [~,dphiDk2] = dTs(tg2,fg(10:12),zeros(3,1),dphifg(10:12));
        
        dphiKt = H2'*dphiKg*H2;
        
        dphiKt(4:6,4:6) = dphiKt(4:6,4:6)+dphiDk1;
        
        dphiKt(10:12,10:12) = dphiKt(10:12,10:12)+dphiDk2;
        
        morph.dphiKt = reshape(dphiKt',[],1);
    end
    %% Calculation of Mtheta
    if isfield(morph,'fold')
        
        if morph.energy
            dp = eye(12,12);
            
            dptg1 = dp(4:6,:);
            dptg2 = dp(10:12,:);
            
            [~,dpRg1] = expon(tg1,dptg1);
            [~,dpRg2] = expon(tg2,dptg2);
            
            dpd21 = dp(7:9,:)-dp(1:3,:);
            dpl = 1/2/l*(2*(x21+d21)'*dpd21);
            dpu = dpl;
            
            dpe1 = 1/l^2*(dpd21.*l-(x21(:,ones(12,1))+d21(:,ones(12,1))).*dpl(ones(3,1),:));
        end
        
        [~,dthetaRtheta,d2thetaRtheta]   = expon(morph.fold*morph.c0,morph.c0);        
        
        if isfield(morph,'shear')
            dthetaRo = dthetaRtheta*Ro2*Rpsi;
            if morph.energy
                d2thetaRo = d2thetaRtheta*Ro2*Rpsi;
                dpsidthetaRo = dthetaRtheta*Ro2*dpsiRpsi;
            end
        else
            dthetaRo = dthetaRtheta*Ro2;
            if morph.energy
                d2thetaRo = d2thetaRtheta*Ro2;
            end
        end
        
        dthetaqb1 = Rg1*dthetaRo*[0;1;0];
        dthetaqb2 = Rg2*dthetaRo*[0;1;0];
        dthetaqb = (dthetaqb1+dthetaqb2)/2;
        
        if morph.energy
            d2thetaqb1 = Rg1*d2thetaRo*[0;1;0];
            d2thetaqb2 = Rg2*d2thetaRo*[0;1;0];
            d2thetaqb = (d2thetaqb1+d2thetaqb2)/2;
            for i=1:size(dp,2)
                dpqb1(:,i)     = dpRg1(:,:,i)*Ro*[0;1;0];
                dpqb2(:,i)     = dpRg2(:,:,i)*Ro*[0;1;0];
                dthetadpqb1(:,i) = dpRg1(:,:,i)*dthetaRo*[0;1;0];
                dthetadpqb2(:,i) = dpRg2(:,:,i)*dthetaRo*[0;1;0];
            end
            dpqb     = (dpqb1+dpqb2)/2;
            dthetadpqb = (dthetadpqb1+dthetadpqb2)/2;
            
            if isfield(morph,'shear')
                dpsidthetaqb1 = Rg1*dpsidthetaRo*[0;1;0];
                dpsidthetaqb2 = Rg2*dpsidthetaRo*[0;1;0];
                dpsidthetaqb = (dpsidthetaqb1+dpsidthetaqb2)/2;
            end
        end
        
        if morph.energy
            [~,dthetae3b,d2thetae3b]    = vec(e1,qb,zeros(size(e1)),dthetaqb,zeros(size(e1)),d2thetaqb);
            [~,dpe3b,~,dthetadpe3b] = vec(e1,qb,dpe1,dpqb,[],[],zeros(size(e1)),dthetaqb,zeros(size(e1,1),size(dpqb,2)),dthetadpqb);
            
            if isfield(morph,'shear')
                [~,~,~,dpsidthetae3b] = vec(e1,qb,zeros(size(e1)),dpsiqb,[],[],zeros(size(e1)),dthetaqb,zeros(size(e1,1),size(dpsiqb,2)),dpsidthetaqb);
            end
            
            dthetane3b = 1/2/norm(e3b)*(2*e3b(1)*dthetae3b(1,:)+2*e3b(2)*dthetae3b(2,:)+2*e3b(3)*dthetae3b(3,:));
            d2thetane3b = 1/2/norm(e3b)^2*((2*(dthetae3b(1,:)*dthetae3b(1,:)+e3b(1)*d2thetae3b(1,:))+2*(dthetae3b(2,:)*...
                dthetae3b(2,:)+e3b(2)*d2thetae3b(2,:))+2*(dthetae3b(3,:)*dthetae3b(3,:)+e3b(3)*d2thetae3b(3,:)))*norm(e3b)-...
                (2*e3b(1)*dthetae3b(1,:)+2*e3b(2)*dthetae3b(2,:)+2*e3b(3)*dthetae3b(3,:))*dthetane3b);
            dthetae3 = 1/norm(e3b)^2*(dthetae3b*norm(e3b)-e3b.*dthetane3b);
            d2thetae3 = 1/norm(e3b)^4*((d2thetae3b*norm(e3b)+dthetae3b*dthetane3b-dthetae3b.*dthetane3b-e3b.*d2thetane3b)*norm(e3b)^2-...
                (dthetae3b*norm(e3b)-e3b.*dthetane3b)*2*norm(e3b)*dthetane3b);
            dpne3b        = 1/2/norm(e3b)*(2*e3b(1).*dpe3b(1,:)+2*e3b(2).*dpe3b(2,:)+2*e3b(3).*dpe3b(3,:));
            dthetadpne3b    = 1/2/norm(e3b)^2*((2*(dthetae3b(1).*dpe3b(1,:)+e3b(1).*dthetadpe3b(1,:))+2*(dthetae3b(2).*dpe3b(2,:)+e3b(2).*...
                dthetadpe3b(2,:))+2*(dthetae3b(3).*dpe3b(3,:)+e3b(3).*dthetadpe3b(3,:)))*norm(e3b)-(2*e3b(1).*dpe3b(1,:)+2*e3b(2).*...
                dpe3b(2,:)+2*e3b(3).*dpe3b(3,:)).*dthetane3b);
            dpe3          = 1/norm(e3b)^2*[dpe3b(1,:)*norm(e3b)-e3b(1).*dpne3b;dpe3b(2,:)*norm(e3b)-e3b(2).*dpne3b;...
                dpe3b(3,:)*norm(e3b)-e3b(3).*dpne3b];
            dthetadpe3      = 1/norm(e3b)^4*([dthetadpe3b(1,:)*norm(e3b)+dpe3b(1,:)*dthetane3b-dthetae3b(1).*dpne3b-e3b(1).*...
                dthetadpne3b;dthetadpe3b(2,:)*norm(e3b)+dpe3b(2,:)*dthetane3b-dthetae3b(2).*dpne3b-e3b(2).*dthetadpne3b;...
                dthetadpe3b(3,:)*norm(e3b)+dpe3b(3,:)*dthetane3b-dthetae3b(3).*dpne3b-e3b(3).*dthetadpne3b]*norm(e3b)^2-...
                [dpe3b(1,:)*norm(e3b)-e3b(1).*dpne3b;dpe3b(2,:)*norm(e3b)-e3b(2).*dpne3b;dpe3b(3,:)*norm(e3b)-...
                e3b(3).*dpne3b]*2*norm(e3b)*dthetane3b);
            
            if isfield(morph,'shear')
                dpsidthetane3b    = 1/2/norm(e3b)^2*((2*(dthetae3b(1).*dpsie3b(1,:)+e3b(1).*dpsidthetae3b(1,:))+2*(dthetae3b(2).*dpsie3b(2,:)+e3b(2).*...
                    dpsidthetae3b(2,:))+2*(dthetae3b(3).*dpsie3b(3,:)+e3b(3).*dpsidthetae3b(3,:)))*norm(e3b)-(2*e3b(1).*dpsie3b(1,:)+2*e3b(2).*...
                    dpsie3b(2,:)+2*e3b(3).*dpsie3b(3,:)).*dthetane3b);
                dpsidthetae3      = 1/norm(e3b)^4*([dpsidthetae3b(1,:)*norm(e3b)+dpsie3b(1,:)*dthetane3b-dthetae3b(1).*dpsine3b-e3b(1).*...
                    dpsidthetane3b;dpsidthetae3b(2,:)*norm(e3b)+dpsie3b(2,:)*dthetane3b-dthetae3b(2).*dpsine3b-e3b(2).*dpsidthetane3b;...
                    dpsidthetae3b(3,:)*norm(e3b)+dpsie3b(3,:)*dthetane3b-dthetae3b(3).*dpsine3b-e3b(3).*dpsidthetane3b]*norm(e3b)^2-...
                    [dpsie3b(1,:)*norm(e3b)-e3b(1).*dpsine3b;dpsie3b(2,:)*norm(e3b)-e3b(2).*dpsine3b;dpsie3b(3,:)*norm(e3b)-...
                    e3b(3).*dpsine3b]*2*norm(e3b)*dthetane3b);
            end
            
            
            [~,dthetae2,d2thetae2] = vec(e3,e1,dthetae3,zeros(size(dthetae3)),d2thetae3,zeros(size(dthetae3)));
            [~,dpe2,~,dthetadpe2] = vec(e3,e1,dpe3,dpe1,[],[],dthetae3,zeros(size(e1)),dthetadpe3,zeros(size(e3,1),size(dpqb,2)));
            if isfield(morph,'shear')
                [~,~,~,dpsidthetae2] = vec(e3,e1,dpsie3,zeros(size(e1)),[],[],dthetae3,zeros(size(e1)),dpsidthetae3,zeros(size(e1)));
            end
        else
            [~,dthetae3b]    = vec(e1,qb,zeros(size(e1)),dthetaqb);
            dthetane3b = 1/2/norm(e3b)*(2*e3b(1)*dthetae3b(1,:)+2*e3b(2)*dthetae3b(2,:)+2*e3b(3)*dthetae3b(3,:));
            dthetae3 = 1/norm(e3b)^2*(dthetae3b*norm(e3b)-e3b.*dthetane3b);
            [~,dthetae2] = vec(e3,e1,dthetae3,zeros(size(dthetae3)));
        end
        
        dthetaRr = [zeros(size(dthetae2)) dthetae2 dthetae3];
        
        if morph.energy
            d2thetaRr = [zeros(size(dthetae3)) d2thetae2 d2thetae3];
            dpRr(:,1,:) = dpe1;
            dpRr(:,2,:) = dpe2;
            dpRr(:,3,:) = dpe3;
            dthetadpRr(:,1,:) = zeros(size(dthetadpe2));
            dthetadpRr(:,2,:) = dthetadpe2;
            dthetadpRr(:,3,:) = dthetadpe3;
            if isfield(morph,'shear')
                dpsidthetaRr = [zeros(size(dpsidthetae3)) dpsidthetae2 dpsidthetae3];
            end
        end
        
        dthetaq = dthetaRr'*q+Rr'*dthetaqb;
        dthetaq1 = dthetaRr'*q1+Rr'*dthetaqb1;
        
        dthetanu = 1/q(2)^2*(dthetaq(1)*q(2)-q(1)*dthetaq(2));
        dthetanu11 = 1/q(2)^2*(dthetaq1(1)*q(2)-q1(1)*dthetaq(2));
        dthetanu12 = 1/q(2)^2*(dthetaq1(2)*q(2)-q1(2)*dthetaq(2));
        dthetanu21 = 2*dthetanu-dthetanu11;
        dthetanu22 = -dthetanu12;
        
        dthetaRe1 = dthetaRr'*Rg1*Ro+Rr'*Rg1*dthetaRo;
        dthetaRe2 = dthetaRr'*Rg2*Ro+Rr'*Rg2*dthetaRo;
        
        if morph.energy
            d2thetaRe1 = d2thetaRr'*Rg1*Ro+dthetaRr'*Rg1*dthetaRo+Rr'*Rg1*d2thetaRo+dthetaRr'*Rg1*dthetaRo;
            d2thetaRe2 = d2thetaRr'*Rg2*Ro+dthetaRr'*Rg2*dthetaRo+Rr'*Rg2*d2thetaRo+dthetaRr'*Rg2*dthetaRo;
            for i=1:size(dp,2)
                dpRe1(:,:,i) = dpRr(:,:,i)'*Rg1*Ro+Rr'*dpRg1(:,:,i)*Ro;
                dpRe2(:,:,i) = dpRr(:,:,i)'*Rg2*Ro+Rr'*dpRg2(:,:,i)*Ro;
                dthetadpRe1(:,:,i) = dthetadpRr(:,:,i)'*Rg1*Ro+dpRr(:,:,i)'*Rg1*dthetaRo+...
                    dthetaRr'*dpRg1(:,:,i)*Ro+Rr'*dpRg1(:,:,i)*dthetaRo;
                dthetadpRe2(:,:,i) = dthetadpRr(:,:,i)'*Rg2*Ro+dpRr(:,:,i)'*Rg2*dthetaRo+...
                    dthetaRr'*dpRg2(:,:,i)*Ro+Rr'*dpRg2(:,:,i)*dthetaRo;
            end
            
            if isfield(morph,'shear')
                dpsidthetaRe1 = dthetaRr'*Rg1*dpsiRo+dpsiRr'*Rg1*dthetaRo+Rr'*Rg1*dpsidthetaRo+dpsidthetaRr'*Rg1*Ro;
                dpsidthetaRe2 = dthetaRr'*Rg2*dpsiRo+dpsiRr'*Rg2*dthetaRo+Rr'*Rg2*dpsidthetaRo+dpsidthetaRr'*Rg2*Ro;
            end
        end
        
        
        if morph.energy
            [~,dthetatl1,d2thetatl1]    = logar(Re1,dthetaRe1,d2thetaRe1);
            [~,dthetatl2,d2thetatl2]    = logar(Re2,dthetaRe2,d2thetaRe2);
            [~,dptl1,~,dthetadptl1] = logar(Re1,dpRe1,[],dthetaRe1,dthetadpRe1);
            [~,dptl2,~,dthetadptl2] = logar(Re2,dpRe2,[],dthetaRe2,dthetadpRe2);
            if isfield(morph,'shear')
                [~,~,~,dpsidthetatl1] = logar(Re1,dpsiRe1,[],dthetaRe1,dpsidthetaRe1);
                [~,~,~,dpsidthetatl2] = logar(Re2,dpsiRe2,[],dthetaRe2,dpsidthetaRe2);
            end
        else
            [~,dthetatl1]    = logar(Re1,dthetaRe1);
            [~,dthetatl2]    = logar(Re2,dthetaRe2);
        end
        morph.dthetatl1 = dthetatl1;
        morph.dthetatl2 = dthetatl2;
        
        if morph.energy
            morph.d2thetatl1  = d2thetatl1;
            morph.d2thetatl2  = d2thetatl2;
            morph.dptl1     = dptl1;
            morph.dptl2     = dptl2;
            morph.dpu       = dpu;
            morph.dthetadptl1 = dthetadptl1;
            morph.dthetadptl2 = dthetadptl2;
            if isfield(morph,'shear')
                morph.dpsitl1     = dpsitl1;
                morph.dpsitl2     = dpsitl2;
                morph.dpsidthetatl1 = dpsidthetatl1;
                morph.dpsidthetatl2 = dpsidthetatl2;
            end
        end
        
        [~,~,~,~,morphoth] = locel(u,tl1,tl2,lo,C1,C2,morph);
        morph.dthetare = morphoth.dredtheta;
        dthetafl = morphoth.dfldtheta;
        
        if morph.energy
            morph.Mtheta = fl'*morphoth.dqedtheta;
            morph.dpMtheta   = morphoth.dqedtheta'*morphoth.dpfl+fl'*morphoth.dqedthetadp;
            morph.dthetaMtheta = morphoth.dfldtheta'*morphoth.dqedtheta+fl'*morphoth.d2qedtheta;
            
            if isfield(morph,'shear')
                morph.dthetaMpsi = morphoth.dfldtheta'*morphosh.dqedpsi + fl'*morphoth.dqedpsidtheta;
                morph.dpsiMtheta = morphosh.dfldpsi'*morphoth.dqedtheta + fl'*morphoth.dqedpsidtheta;
            end
            if isfield(morph,'twist')
                morph.dthetaMphi = morphoth.dfldtheta'*morphotw.dqedphi;
                morph.dphiMtheta = morphotw.dfldphi'*morphoth.dqedtheta;
            end
        end
        
        [~,dthetaDe1] = invTs(tl1,dthetatl1);
        [~,dthetaDe2] = invTs(tl2,dthetatl2);
        %
        dthetafe = [dthetafl(1)
            dthetaDe1'*fl(2:4)+De1'*dthetafl(2:4)
            dthetaDe2'*fl(5:7)+De2'*dthetafl(5:7)];
        
        dthetaH1 = [0   O1     O1
            O1' dthetaDe1 O3
            O1' O3     dthetaDe2];
        
        [Dh1_dum,dthetaDh1_dum] = dinvTs(tl1,fl(2:4),dthetatl1,dthetafl(2:4));
        [Dh2_dum,dthetaDh2_dum] = dinvTs(tl2,fl(5:7),dthetatl2,dthetafl(5:7));
        dthetaDh1 = dthetaDh1_dum*De1+Dh1_dum*dthetaDe1;
        dthetaDh2 = dthetaDh2_dum*De2+Dh2_dum*dthetaDe2;
        
        dthetaKh = [0   O1     O1
            O1' dthetaDh1 O3
            O1' O3     dthetaDh2];
        
        dthetake = dthetaH1'*kl*H1+H1'*kl*dthetaH1+dthetaKh;
        
        % transformation to the global coordinates
        
        dthetaB = [ zeros(size(r'))
            -1/l*(dthetanu*e3'+nu*dthetae3') (-dthetanu12/2)*e1'+1/2*(dthetanu11*e2'+nu11*dthetae2')  1/l*(dthetanu*e3'+nu*dthetae3') 1/2*(-dthetanu22*e1'+dthetanu21*e2'+nu21*dthetae2')
            -dthetae3'/l dthetae2' dthetae3'/l 0 0 0
            dthetae2'/l dthetae3' -dthetae2'/l 0 0 0
            -1/l*(dthetanu*e3'+nu*dthetae3') 1/2*(-dthetanu12*e1'+dthetanu11*e2'+nu11*dthetae2')  1/l*(dthetanu*e3'+nu*dthetae3') (-dthetanu22/2)*e1'+1/2*(dthetanu21*e2'+nu21*dthetae2')
            -dthetae3'/l 0 0 0 dthetae3'/l dthetae2'
            dthetae2'/l 0 0 0 -dthetae2'/l dthetae3'];
        
        dthetafg = dthetaB'*fe+B'*dthetafe;
        
        dthetaG = [0   0           dthetanu/l  dthetanu12/2  -dthetanu11/2  0  0       0        -dthetanu/l  dthetanu22/2  -dthetanu21/2  0
            0   0           0                   0           0     0  0       0         0                        0        0     0
            0  0         0                       0           0     0  0  0           0                          0        0     0]';
        
        dthetaP = -[dthetaG';dthetaG'];
        
        dthetaF = dthetaP'*fe(2:7)+P'*dthetafe(2:7);
        
        [~,dthetasF(1:3,:)]   = skew(F(1:3),dthetaF(1:3));
        [~,dthetasF(4:6,:)]   = skew(F(4:6),dthetaF(4:6));
        [~,dthetasF(7:9,:)]   = skew(F(7:9),dthetaF(7:9));
        [~,dthetasF(10:12,:)] = skew(F(10:12),dthetaF(10:12));
        
        dthetaEE = [dthetaRr O3    O3    O3
            O3    dthetaRr O3    O3
            O3    O3    dthetaRr O3
            O3    O3    O3    dthetaRr];
        
        dthetanab = [0
            1/l*(dthetanu*(fe(2)+fe(5))+nu*(dthetafe(2)+dthetafe(5))+dthetafe(3)+dthetafe(6))
            1/l*(dthetafe(4)+dthetafe(7))];
        
        dthetaKg = dthetaB'*ke*B+B'*dthetake*B+B'*ke*dthetaB+Dr*dthetafe(1)-dthetaEE*sF*G'*EE'-EE*dthetasF*G'*EE'-EE*sF*dthetaG'*EE'-EE*sF*G'*dthetaEE'+dthetaEE*G*nab*r'+EE*dthetaG*nab*r'+EE*G*dthetanab*r';
        
        % transformation to the new global coordinates
        
        morph.dthetaft = [dthetafg(1:3)
            Dg1'*dthetafg(4:6)
            dthetafg(7:9)
            Dg2'*dthetafg(10:12)];
        
        [~,dthetaDk1] = dTs(tg1,fg(4:6),zeros(3,1),dthetafg(4:6));
        [~,dthetaDk2] = dTs(tg2,fg(10:12),zeros(3,1),dthetafg(10:12));
        
        dthetaKt = H2'*dthetaKg*H2;
        
        dthetaKt(4:6,4:6) = dthetaKt(4:6,4:6)+dthetaDk1;
        
        dthetaKt(10:12,10:12) = dthetaKt(10:12,10:12)+dthetaDk2;
        
        morph.dthetaKt = reshape(dthetaKt',[],1);
    end
    
    
    if isfield(morph,'span')
        dx = eye(6);
        morph.dx_Kt = sparse(numel(Kt),numel(x));
        morph.dx_ft = sparse(numel(ft),numel(x));
        morph.dx_Kl = sparse(numel(Ke),numel(x));
        for i=1:6
            
            dx_x21 = dx(4:6,i)-dx(1:3,i);
            
            dx_lo = 1/2/lo*((dx_x21)'*(x21)+(x21)'*(dx_x21));
            dx_l = 1/2/l*((dx_x21)'*(x21+d21)+(x21+d21)'*(dx_x21));
            
            dx_u = dx_l-dx_lo;
            
            % rigid rotation
            
            dx_e1 = 1/l^2*((dx_x21)*l-(x21+d21)*dx_l);
            
            [~,dx_e3b] = vec(e1,qb,dx_e1,zeros(3,1));
            dx_ne3b      = 1/2/norm(e3b)*(2*e3b(1)*dx_e3b(1)+2*e3b(2)*dx_e3b(2)+2*e3b(3)*dx_e3b(3));
            dx_e3        = 1/norm(e3b)^2*(dx_e3b*norm(e3b)-e3b*dx_ne3b);
            
            [~,dx_e2]  = vec(e3,e1,dx_e3,dx_e1);
            
            dx_Rr = [dx_e1 dx_e2 dx_e3];
            
            dx_q  = dx_Rr'*qb;
            dx_q1 = dx_Rr'*qb1;
            
            dx_nu   = 1/q(2)^2*(dx_q(1)*q(2)-q(1)*dx_q(2));
            dx_nu11 = 1/q(2)^2*(dx_q1(1)*q(2)-q1(1)*dx_q(2));
            dx_nu12 = 1/q(2)^2*(dx_q1(2)*q(2)-q1(2)*dx_q(2));
            dx_nu21 = 2*dx_nu-dx_nu11;
            dx_nu22 = -dx_nu12;
            
            %local rotations
            
            dx_Re1 = (dx_Rr'*Rg1)*Ro;
            dx_Re2 = (dx_Rr'*Rg2)*Ro;
            
            [~,dx_tl1] = logar(Re1,dx_Re1);
            [~,dx_tl2] = logar(Re2,dx_Re2);
            
            % local force vector and tangent stiffness matrix
            
            [pdx_fl,pdx_rei] = dplocel(dx_u,dx_tl1,dx_tl2,lo,C1,C2);
            
            if isfield(morph,'twist')
                [dlo_kl,dlo_fl,dlo_re,dlo_Ke] = dllocel(u,tl1,tl2,lo,C1,C2,morph);
            else
                [dlo_kl,dlo_fl,dlo_re,dlo_Ke] = dllocel(u,tl1,tl2,lo,C1,C2);
            end
            
            dx_fl = pdx_fl+dlo_fl*dx_lo;
            dx_rei = pdx_rei+dlo_re*dx_lo;
            dx_kl = dlo_kl*dx_lo;
            
            dx_Ke = dlo_Ke*dx_lo;
            morph.dx_Kl(:,i) = reshape(dx_Ke',[],1);
            morph.dx_re(:,i) = dx_rei;
            
            if morph.energy
                if isfield(morph,'twist')
                    morph.dxMphi(:,i) = dx_fl'*morphotw.dqedphi;
                end
                if isfield(morph,'shear')
                    morph.dxMpsi(:,i) = dx_fl'*morphosh.dqedpsi;
                end
                if isfield(morph,'fold')
                    morph.dxMtheta(:,i) = dx_fl'*morphoth.dqedtheta;
                end
            end
            
            % transformation to the new local coordinates
            
            [~,dx_De1] = invTs(tl1,dx_tl1);
            [~,dx_De2] = invTs(tl2,dx_tl2);
            
            dx_fe=[dx_fl(1)
                dx_De1'*fl(2:4)+De1'*dx_fl(2:4)
                dx_De2'*fl(5:7)+De2'*dx_fl(5:7)];
            
            dx_H1 = [0   O1     O1
                O1' dx_De1 O3
                O1' O3     dx_De2];
            
            [Dh1_dum,dx_Dh1_dum] = dinvTs(tl1,fl(2:4),dx_tl1,dx_fl(2:4));
            [Dh2_dum,dx_Dh2_dum] = dinvTs(tl2,fl(5:7),dx_tl2,dx_fl(5:7));
            dx_Dh1 = dx_Dh1_dum*De1+Dh1_dum*dx_De1;
            dx_Dh2 = dx_Dh2_dum*De2+Dh2_dum*dx_De2;
            
            dx_Kh = [0   O1     O1
                O1' dx_Dh1 O3
                O1' O3     dx_Dh2];
            
            dx_ke = dx_H1'*kl*H1+H1'*dx_kl*H1+H1'*kl*dx_H1+dx_Kh;
            
            % transformation to the global coordinates
            
            dx_r = [-dx_e1;0;0;0;dx_e1;0;0;0];
            
            dx_B = [ dx_r'
                -(1/l^2*(dx_nu*l-nu*dx_l)*e3'+nu/l*dx_e3'), (-dx_nu12/2)*e1'+(1-nu12/2)*dx_e1'+dx_nu11/2*e2'+nu11/2*dx_e2',  1/l^2*(dx_nu*l-nu*dx_l)*e3'+nu/l*dx_e3', 1/2*(-(dx_nu22*e1'+nu22*dx_e1')+dx_nu21*e2'+nu21*dx_e2')
                -1/l^2*(dx_e3'*l-e3'*dx_l), dx_e2', 1/l^2*(dx_e3'*l-e3'*dx_l), 0, 0, 0
                1/l^2*(dx_e2'*l-e2'*dx_l), dx_e3', -1/l^2*(dx_e2'*l-e2'*dx_l), 0, 0, 0
                -(1/l^2*(dx_nu*l-nu*dx_l)*e3'+nu/l*dx_e3'), 1/2*(-dx_nu12*e1'-nu12*dx_e1'+dx_nu11*e2'+nu11*dx_e2') , (1/l^2*(dx_nu*l-nu*dx_l))*e3'+nu/l*dx_e3', (-dx_nu22/2)*e1'+(1-nu22/2)*dx_e1'+dx_nu21/2*e2'+nu21/2*dx_e2'
                -1/l^2*(dx_e3'*l-e3'*dx_l), 0, 0 ,0 ,1/l^2*(dx_e3'*l-e3'*dx_l), dx_e2'
                1/l^2*(dx_e2'*l-e2'*dx_l), 0, 0, 0, -1/l^2*(dx_e2'*l-e2'*dx_l), dx_e3'];
            
            dx_fg = dx_B'*fe+B'*dx_fe;
            
            dx_A = 1/l^2*((-dx_e1*e1'-e1*dx_e1')*l-(I3-e1*e1')*dx_l);
            
            dx_Dr = [dx_A  O3 -dx_A  O3
                O3    O3  O3    O3
                -dx_A  O3  dx_A  O3
                O3    O3  O3    O3];
            
            dx_G = [0   0           1/l^2*(dx_nu*l-nu*dx_l)  dx_nu12/2  -dx_nu11/2  0  0       0        -1/l^2*(dx_nu*l-nu*dx_l)  dx_nu22/2  -dx_nu21/2  0
                0   0           -1/l^2*dx_l                   0           0     0  0       0         1/l^2*dx_l                        0        0     0
                0  1/l^2*dx_l         0                       0           0     0  0  -1/l^2*dx_l           0                          0        0     0]';
            
            dx_P = -[dx_G';dx_G'];
            
            dx_F = dx_P'*fe(2:7)+P'*dx_fe(2:7);
            
            [~,dx_sF(1:3,:)]   = skew(F(1:3),dx_F(1:3));
            [~,dx_sF(4:6,:)]   = skew(F(4:6),dx_F(4:6));
            [~,dx_sF(7:9,:)]   = skew(F(7:9),dx_F(7:9));
            [~,dx_sF(10:12,:)] = skew(F(10:12),dx_F(10:12));
            
            dx_EE = [dx_Rr O3    O3    O3
                O3    dx_Rr O3    O3
                O3    O3    dx_Rr O3
                O3    O3    O3    dx_Rr];
            
            dx_nab = [0
                1/l^2*((dx_nu*(fe(2)+fe(5))+nu*(dx_fe(2)+dx_fe(5))+dx_fe(3)+dx_fe(6))*l-(nu*(fe(2)+fe(5))+fe(3)+fe(6))*dx_l)
                1/l^2*((dx_fe(4)+dx_fe(7))*l-(fe(4)+fe(7))*dx_l)];
            
            dx_Kg = dx_B'*ke*B+B'*dx_ke*B+B'*ke*dx_B+dx_Dr*fe(1)+Dr*dx_fe(1)-dx_EE*sF*G'*EE'-EE*dx_sF*G'*EE'-EE*sF*dx_G'*EE'-EE*sF*G'*dx_EE'+dx_EE*G*nab*r'+EE*dx_G*nab*r'+EE*G*dx_nab*r'+EE*G*nab*dx_r';
            
            % transformation to the new global coordinates
            dx_ft = [dx_fg(1:3)
                Dg1'*dx_fg(4:6)
                dx_fg(7:9)
                Dg2'*dx_fg(10:12)];
            
            morph.dx_ft(:,i) = dx_ft;
            
            [~,dx_Dk1] = dTs(tg1,fg(4:6),zeros(3,1),dx_fg(4:6));
            [~,dx_Dk2] = dTs(tg2,fg(10:12),zeros(3,1),dx_fg(10:12));
            
            dx_Kt = H2'*dx_Kg*H2;
            
            dx_Kt(4:6,4:6) = dx_Kt(4:6,4:6)+dx_Dk1;
            
            dx_Kt(10:12,10:12) = dx_Kt(10:12,10:12)+dx_Dk2;
            
            morph.dx_Kt(:,i) = reshape(dx_Kt',[],1);
        end
    end
    
    
    varargout{1} = morph;
end

if ders==1
    
    %% Derivative of Kt and ft wrt p
    dp = eye(12);
    dp_Kt_s = sparse(numel(Kt),numel(p));
    dp_ft_s = sparse(numel(ft),numel(p));
    
    for i=1:12
        dp_tg1 = dp(4:6,i);
        dp_tg2 = dp(10:12,i);
        
        [~,dp_Rg1] = expon(tg1,dp_tg1);
        [~,dp_Rg2] = expon(tg2,dp_tg2);
        
        dp_d21 = dp(7:9,i)-dp(1:3,i);
        
        dp_l = 1/2/l*((dp_d21)'*(x21+d21)+(x21+d21)'*(dp_d21));
        dp_u = dp_l;
        
        % rigid rotation
        
        dp_e1 = 1/l^2*((dp_d21)*l-(x21+d21)*dp_l);
        
        dp_qb1 = dp_Rg1*Ro*[0;1;0];
        dp_qb2 = dp_Rg2*Ro*[0;1;0];
        dp_qb  = (dp_qb1+dp_qb2)/2;
        
        [~,dp_e3b] = vec(e1,qb,dp_e1,dp_qb);
        dp_ne3b      = 1/2/norm(e3b)*(2*e3b(1)*dp_e3b(1)+2*e3b(2)*dp_e3b(2)+2*e3b(3)*dp_e3b(3));
        dp_e3        = 1/norm(e3b)^2*(dp_e3b*norm(e3b)-e3b*dp_ne3b);
        
        [~,dp_e2]  = vec(e3,e1,dp_e3,dp_e1);
        
        dp_Rr = [dp_e1 dp_e2 dp_e3];
        
        dp_q  = dp_Rr'*qb+Rr'*dp_qb;
        dp_q1 = dp_Rr'*qb1+Rr'*dp_qb1;
        
        dp_nu   = 1/q(2)^2*(dp_q(1)*q(2)-q(1)*dp_q(2));
        dp_nu11 = 1/q(2)^2*(dp_q1(1)*q(2)-q1(1)*dp_q(2));
        dp_nu12 = 1/q(2)^2*(dp_q1(2)*q(2)-q1(2)*dp_q(2));
        dp_nu21 = 2*dp_nu-dp_nu11;
        dp_nu22 = -dp_nu12;
        
        %local rotations
        
        dp_Re1 = (dp_Rr'*Rg1+Rr'*dp_Rg1)*Ro;
        dp_Re2 = (dp_Rr'*Rg2+Rr'*dp_Rg2)*Ro;
        
        [~,dp_tl1] = logar(Re1,dp_Re1);
        [~,dp_tl2] = logar(Re2,dp_Re2);
        
        % local force vector and tangent stiffness matrix
        
        [dp_fl,dp_rei] = dplocel(dp_u,dp_tl1,dp_tl2,lo,C1,C2);
        dp_re(:,i) = dp_rei;
        
        % transformation to the new local coordinates
        
        [~,dp_De1] = invTs(tl1,dp_tl1);
        [~,dp_De2] = invTs(tl2,dp_tl2);
        
        dp_fe=[dp_fl(1)
            dp_De1'*fl(2:4)+De1'*dp_fl(2:4)
            dp_De2'*fl(5:7)+De2'*dp_fl(5:7)];
        
        dp_H1 = [0   O1     O1
            O1' dp_De1 O3
            O1' O3     dp_De2];
        
        [Dh1_dum,dp_Dh1_dum] = dinvTs(tl1,fl(2:4),dp_tl1,dp_fl(2:4));
        [Dh2_dum,dp_Dh2_dum] = dinvTs(tl2,fl(5:7),dp_tl2,dp_fl(5:7));
        dp_Dh1 = dp_Dh1_dum*De1+Dh1_dum*dp_De1;
        dp_Dh2 = dp_Dh2_dum*De2+Dh2_dum*dp_De2;
        
        dp_Kh = [0   O1     O1
            O1' dp_Dh1 O3
            O1' O3     dp_Dh2];
        
        dp_ke = dp_H1'*kl*H1+H1'*kl*dp_H1+dp_Kh;
        
        % transformation to the global coordinates
        
        dp_r = [-dp_e1;0;0;0;dp_e1;0;0;0];
        
        dp_B = [ dp_r'
            -(1/l^2*(dp_nu*l-nu*dp_l)*e3'+nu/l*dp_e3'), (-dp_nu12/2)*e1'+(1-nu12/2)*dp_e1'+dp_nu11/2*e2'+nu11/2*dp_e2',  1/l^2*(dp_nu*l-nu*dp_l)*e3'+nu/l*dp_e3', 1/2*(-(dp_nu22*e1'+nu22*dp_e1')+dp_nu21*e2'+nu21*dp_e2')
            -1/l^2*(dp_e3'*l-e3'*dp_l), dp_e2', 1/l^2*(dp_e3'*l-e3'*dp_l), 0, 0, 0
            1/l^2*(dp_e2'*l-e2'*dp_l), dp_e3', -1/l^2*(dp_e2'*l-e2'*dp_l), 0, 0, 0
            -(1/l^2*(dp_nu*l-nu*dp_l)*e3'+nu/l*dp_e3'), 1/2*(-dp_nu12*e1'-nu12*dp_e1'+dp_nu11*e2'+nu11*dp_e2') , (1/l^2*(dp_nu*l-nu*dp_l))*e3'+nu/l*dp_e3', (-dp_nu22/2)*e1'+(1-nu22/2)*dp_e1'+dp_nu21/2*e2'+nu21/2*dp_e2'
            -1/l^2*(dp_e3'*l-e3'*dp_l), 0, 0 ,0 ,1/l^2*(dp_e3'*l-e3'*dp_l), dp_e2'
            1/l^2*(dp_e2'*l-e2'*dp_l), 0, 0, 0, -1/l^2*(dp_e2'*l-e2'*dp_l), dp_e3'];
        %         var = B;
        %         dvar(:,i) = reshape(dp_B',[],1);
        
        dp_fg = dp_B'*fe+B'*dp_fe;
        
        dp_A = 1/l^2*((-dp_e1*e1'-e1*dp_e1')*l-(I3-e1*e1')*dp_l);
        
        dp_Dr = [dp_A  O3 -dp_A  O3
            O3    O3  O3    O3
            -dp_A  O3  dp_A  O3
            O3    O3  O3    O3];
        
        dp_G = [0   0           1/l^2*(dp_nu*l-nu*dp_l)  dp_nu12/2  -dp_nu11/2  0  0       0        -1/l^2*(dp_nu*l-nu*dp_l)  dp_nu22/2  -dp_nu21/2  0
            0   0           -1/l^2*dp_l                   0           0     0  0       0         1/l^2*dp_l                        0        0     0
            0  1/l^2*dp_l         0                       0           0     0  0  -1/l^2*dp_l           0                          0        0     0]';
        
        dp_P = -[dp_G';dp_G'];
        
        dp_F = dp_P'*fe(2:7)+P'*dp_fe(2:7);
        
        [~,dp_sF(1:3,:)]   = skew(F(1:3),dp_F(1:3));
        [~,dp_sF(4:6,:)]   = skew(F(4:6),dp_F(4:6));
        [~,dp_sF(7:9,:)]   = skew(F(7:9),dp_F(7:9));
        [~,dp_sF(10:12,:)] = skew(F(10:12),dp_F(10:12));
        
        dp_EE = [dp_Rr O3    O3    O3
            O3    dp_Rr O3    O3
            O3    O3    dp_Rr O3
            O3    O3    O3    dp_Rr];
        
        dp_nab = [0
            1/l^2*((dp_nu*(fe(2)+fe(5))+nu*(dp_fe(2)+dp_fe(5))+dp_fe(3)+dp_fe(6))*l-(nu*(fe(2)+fe(5))+fe(3)+fe(6))*dp_l)
            1/l^2*((dp_fe(4)+dp_fe(7))*l-(fe(4)+fe(7))*dp_l)];
        
        dp_Kg = dp_B'*ke*B+B'*dp_ke*B+B'*ke*dp_B+dp_Dr*fe(1)+Dr*dp_fe(1)-dp_EE*sF*G'*EE'-EE*dp_sF*G'*EE'-EE*sF*dp_G'*EE'-EE*sF*G'*dp_EE'+dp_EE*G*nab*r'+EE*dp_G*nab*r'+EE*G*dp_nab*r'+EE*G*nab*dp_r';
        
        % transformation to the new global coordinates
        
        [~,dp_Dg1] = Ts(tg1,dp_tg1);
        [~,dp_Dg2] = Ts(tg2,dp_tg2);
        
        dp_ft = [dp_fg(1:3)
            dp_Dg1'*fg(4:6)+Dg1'*dp_fg(4:6)
            dp_fg(7:9)
            dp_Dg2'*fg(10:12)+Dg2'*dp_fg(10:12)];
        
        dp_ft_s(:,i) = dp_ft;
        
        [~,dp_Dk1] = dTs(tg1,fg(4:6),dp_tg1,dp_fg(4:6));
        [~,dp_Dk2] = dTs(tg2,fg(10:12),dp_tg2,dp_fg(10:12));
        
        dp_H2 = [O3 O3     O3 O3
            O3 dp_Dg1 O3 O3
            O3 O3     O3 O3
            O3 O3     O3 dp_Dg2];
        
        dp_Kt = dp_H2'*Kg*H2+H2'*dp_Kg*H2+H2'*Kg*dp_H2;
        
        dp_Kt(4:6,4:6) = dp_Kt(4:6,4:6)+dp_Dk1;
        
        dp_Kt(10:12,10:12) = dp_Kt(10:12,10:12)+dp_Dk2;
        
        dp_Kt_s(:,i) = reshape(dp_Kt',[],1);
    end
    
    %% Derivative of ft and Kt wrt fl
    
    dfl_fl = eye(7);
    dfl_fl_s = sparse(dfl_fl);
    
    dfl_Kt_s = sparse(numel(Kt),numel(fl));
    dfl_ft_s = sparse(numel(ft),numel(fl));
    
    for i=1:7
        
        dfl_fe = [dfl_fl_s(1,i)
            De1'*dfl_fl_s(2:4,i)
            De2'*dfl_fl_s(5:7,i)];
        
        [~,dfl_Dh1] = dinvTs(tl1,fl(2:4),zeros(size(tl1)),dfl_fl_s(2:4,i));
        dfl_Dh1       = dfl_Dh1*De1;
        
        [~,dfl_Dh2] = dinvTs(tl2,fl(5:7),zeros(size(tl2)),dfl_fl_s(5:7,i));
        dfl_Dh2       = dfl_Dh2*De2;
        
        dfl_Kh = [0      O1     O1
            O1' dfl_Dh1   O3
            O1'    O3   dfl_Dh2];
        
        dfl_ke = dfl_Kh;
        
        % transformation to the global coordinates
        
        dfl_fg = B'*dfl_fe;
        
        dfl_F  = P'*dfl_fe(2:7);
        
        dfl_sF = [skew(dfl_F(1:3))
            skew(dfl_F(4:6))
            skew(dfl_F(7:9))
            skew(dfl_F(10:12))];
        
        dfl_nab=[0
            (nu*(dfl_fe(2)+dfl_fe(5))+dfl_fe(3)+dfl_fe(6))/l
            (dfl_fe(4)+dfl_fe(7))/l];
        
        dfl_Kg = B'*dfl_ke*B+Dr*dfl_fe(1)-EE*dfl_sF*G'*EE'+EE*G*dfl_nab*r';
        
        % transformation to the new global coordinates
        
        dfl_ft_s(:,i) = [dfl_fg(1:3)
            Dg1'*dfl_fg(4:6)
            dfl_fg(7:9)
            Dg2'*dfl_fg(10:12)];
        
        [~,dfl_Dk1] = dTs(tg1,fg(4:6),zeros(size(tg1)),dfl_fg(4:6));
        [~,dfl_Dk2] = dTs(tg2,fg(10:12),zeros(size(tg2)),dfl_fg(10:12));
        
        dfl_Kt    = H2'*dfl_Kg*H2;
        
        dfl_Kt(4:6,4:6) = dfl_Kt(4:6,4:6)+dfl_Dk1;
        
        dfl_Kt(10:12,10:12)  = dfl_Kt(10:12,10:12)+dfl_Dk2;
        dfl_Kt_s(:,i)        = reshape(dfl_Kt',[],1);
    end
    
    %% Derivative of Kt wrt kl
    
    dkl_Kt_s = sparse(numel(Kt),numel(kl));
    
    for i=1:7
        for j=1:7
            dkl_kl          = sparse(7,7);
            dkl_kl(i,j)     = 1;
            dkl_ke          = H1'*dkl_kl*H1;
            
            % transformation to the global coordinates
            
            dkl_Kg          = B'*dkl_ke*B;
            
            % transformation to the new global coordinates
            
            dkl_Kt          = H2'*dkl_Kg*H2;
            dkl_Kt_s(:,j+(i-1)*7) = reshape(dkl_Kt',[],1);
        end
    end
    
    if tailflag == 1
        %% Derivative of kl and fl wrt C
        
        dC1 = sparse(6,6);
        dC2 = sparse(6,6);
        
        dC1_kl_s = sparse(numel(kl),numel(C1));
        dC2_kl_s = sparse(numel(kl),numel(C2));
        
        dC1_fl_s = sparse(numel(fl),numel(C1));
        dC2_fl_s = sparse(numel(fl),numel(C2));
        
        if morph_flag == 1 && morph.energy
            if isfield(morph,'twist')
                dC1_Mphi = sparse(numel(morph.Mphi),numel(C1));
                dC2_Mphi = sparse(numel(morph.Mphi),numel(C2));
            end
            if isfield(morph,'shear')
                dC1_Mpsi = sparse(numel(morph.Mpsi),numel(C1));
                dC2_Mpsi = sparse(numel(morph.Mpsi),numel(C2));
            end
            if isfield(morph,'fold')
                dC1_Mtheta = sparse(numel(morph.Mtheta),numel(C1));
                dC2_Mtheta = sparse(numel(morph.Mtheta),numel(C2));
            end
        end
        
        for i=1:6
            for j=1:6
                dC1(i,j) = 1;
                if morph_flag == 1
                    if isfield(morph,'twist')
                        [dC1_kl,dC1_fl,dC1_rei,dC1_Ke] = dlocel(u,tl1,tl2,lo,C1,C2,dC1,dC2,morph);
                    else
                        [dC1_kl,dC1_fl,dC1_rei,dC1_Ke] = dlocel(u,tl1,tl2,lo,C1,C2,dC1,dC2);
                    end
                elseif morph_flag == 0
                    [dC1_kl,dC1_fl,dC1_rei,dC1_Ke] = dlocel(u,tl1,tl2,lo,C1,C2,dC1,dC2);
                end
                dC1(i,j) = 0;
                
                dC2(i,j) = 1;
                if morph_flag == 1
                    if isfield(morph,'twist')
                        [dC2_kl,dC2_fl,dC2_rei,dC2_Ke] = dlocel(u,tl1,tl2,lo,C1,C2,dC1,dC2,morph);
                    else
                        [dC2_kl,dC2_fl,dC2_rei,dC2_Ke] = dlocel(u,tl1,tl2,lo,C1,C2,dC1,dC2);
                    end
                elseif morph_flag == 0
                    [dC2_kl,dC2_fl,dC2_rei,dC2_Ke] = dlocel(u,tl1,tl2,lo,C1,C2,dC1,dC2);
                end
                dC2(i,j) = 0;
                
                if morph_flag == 1 && morph.energy
                    if isfield(morph,'twist')
                        dC1_Mphi(:,j+(i-1)*6) = dC1_fl'*morphotw.dqedphi;
                        dC2_Mphi(:,j+(i-1)*6) = dC2_fl'*morphotw.dqedphi;
                    end
                    if isfield(morph,'shear')
                        dC1_Mpsi(:,j+(i-1)*6) = dC1_fl'*morphosh.dqedpsi;
                        dC2_Mpsi(:,j+(i-1)*6) = dC2_fl'*morphosh.dqedpsi;
                    end
                    if isfield(morph,'fold')
                        dC1_Mtheta(:,j+(i-1)*6) = dC1_fl'*morphoth.dqedtheta;
                        dC2_Mtheta(:,j+(i-1)*6) = dC2_fl'*morphoth.dqedtheta;
                    end
                end
                
                dC1_kl_s(:,j+(i-1)*6) = reshape(dC1_kl',[],1);
                dC2_kl_s(:,j+(i-1)*6) = reshape(dC2_kl',[],1);
                dC1_Ke_s(:,j+(i-1)*6) = reshape(dC1_Ke',[],1);
                dC2_Ke_s(:,j+(i-1)*6) = reshape(dC2_Ke',[],1);
                dC1_fl_s(:,j+(i-1)*6) = dC1_fl;
                dC2_fl_s(:,j+(i-1)*6) = dC2_fl;
                pdC1_re(:,j+(i-1)*6)  = dC1_rei;
                pdC2_re(:,j+(i-1)*6)  = dC2_rei;
            end
        end
        
        pdC1_Kt = dkl_Kt_s*dC1_kl_s+dfl_Kt_s*dC1_fl_s;
        pdC2_Kt = dkl_Kt_s*dC2_kl_s+dfl_Kt_s*dC2_fl_s;
        
        pdC1_ft = dfl_ft_s*dC1_fl_s;
        pdC2_ft = dfl_ft_s*dC2_fl_s;
        
        dC1_Ke = dC1_Ke_s;
        dC2_Ke = dC2_Ke_s;
        
        if morph_flag == 1 && morph.energy
            if isfield(morph,'twist')
                morph.dC1_Mphi = dC1_Mphi;
                morph.dC2_Mphi = dC2_Mphi;
            end
            if isfield(morph,'shear')
                morph.dC1_Mpsi = dC1_Mpsi;
                morph.dC2_Mpsi = dC2_Mpsi;
            end
            if isfield(morph,'fold')
                morph.dC1_Mtheta = dC1_Mtheta;
                morph.dC2_Mtheta = dC2_Mtheta;
            end
            varargout{1} = morph;
        end
    else
        pdC1_Kt = 0;
        pdC2_Kt = 0;
        pdC1_ft = 0;
        pdC2_ft = 0;
        pdC1_re = 0;
        pdC2_re = 0;
        dC1_Ke = 0;
        dC2_Ke = 0;
    end
else
    pdC1_Kt = 0;
    pdC2_Kt = 0;
    pdC1_ft = 0;
    pdC2_ft = 0;
    pdC1_re = 0;
    pdC2_re = 0;
    dC1_Ke = 0;
    dC2_Ke = 0;
    dp_Kt_s = 0;
    dp_ft_s = 0;
    dp_re = 0;
end