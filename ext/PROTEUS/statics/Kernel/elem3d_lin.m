 function [re,pdC1_re,pdC2_re,dp_re] = elem3d_lin(x,p,C1,C2,Ro,ders)

tg1 = p(4:6);
tg2 = p(10:12);

ug1 = p(1:3);
ug2 = p(7:9);

ul1 = Ro'*ug1;
ul2 = Ro'*ug2;

tl1 = Ro'*tg1;
tl2 = Ro'*tg2;

x21 = x(4:6)-x(1:3);

lo  = sqrt(x21'*x21);

% local force vector and tangent stiffness matrix

[~,~,re] = locel_lin(ul1,ul2,tl1,tl2,lo,C1,C2);

if ders==1

%% Derivative of Kt and ft wrt p
    dp = eye(12);
    
    for i=1:12
        dp_tg1 = dp(4:6,i);
        dp_tg2 = dp(10:12,i);
        
        dp_ug1 = dp(1:3,i);
        dp_ug2 = dp(7:9,i);

        dp_ul1 = Ro'*dp_ug1;
        dp_ul2 = Ro'*dp_ug2;
        
        dp_tl1 = Ro'*dp_tg1;
        dp_tl2 = Ro'*dp_tg2;
        
        % local force vector and tangent stiffness matrix
        
        [~,dp_rei] = dplocel_lin(dp_ul1,dp_ul2,dp_tl1,dp_tl2,lo,C1,C2);
        dp_re(:,i) = dp_rei;
    end
    
%% Derivative of kl and fl wrt C

    dC1 = sparse(6,6);
    dC2 = sparse(6,6);

    for i=1:6
        for j=1:6
            dC1(i,j) = 1;
            [~,~,dC1_rei] = dlocel_lin(ul1,ul2,tl1,tl2,lo,C1,C2,dC1,dC2);
            dC1(i,j) = 0;

            dC2(i,j) = 1;
            [~,~,dC2_rei] = dlocel_lin(ul1,ul2,tl1,tl2,lo,C1,C2,dC1,dC2);
            dC2(i,j) = 0;

            pdC1_re(:,j+(i-1)*6)  = dC1_rei;
            pdC2_re(:,j+(i-1)*6)  = dC2_rei;
        end
    end
    
else
    pdC1_re = 0;
    pdC2_re = 0;
    dp_re   = 0;
end