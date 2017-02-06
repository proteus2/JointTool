function statics = pgen_str(constant,statics,ders,morphflag)

%% Structural data
EFT      = constant.str.EFT;
p        = statics.str.p;
Ne       = constant.str.Ns;
Ndof     = constant.str.Ndof;

dispind = reshape([EFT(:,1:3);EFT(end,7:9)]',3*(Ne+1),1);
rotind  = reshape([EFT(:,4:6);EFT(end,10:12)]',3*(Ne+1),1);
statics.str.dispind = dispind;
statics.str.rotind  = rotind;

if ders == 1
    dp_p = eye(Ndof);
end
%% Interpolation of the structural deformations

disp     = reshape(p(dispind),3,Ne+1);
rot      = reshape(p(rotind),3,Ne+1);

if ders == 1
    dp_disp2 = dp_p(dispind,:);
    dp_rot2  = dp_p(rotind,:);
    
    for i=1:Ne+1
        dp_disp(1:3,i,:) = dp_disp2((1:3)+3*(i-1),:);
        dp_rot(1:3,i,:)  = dp_rot2((1:3)+3*(i-1),:);
    end
end

if morphflag == 1
    if isfield(statics.morph,'span')
        x = reshape(statics.morph.span.xyz(:,end),3,Ne+1);
        statics.str.x = x;
    else
        x = reshape(constant.str.xyz,3,Ne+1);
        statics.str.x = x;
    end
elseif  morphflag == 0
    x = reshape(constant.str.xyz,3,Ne+1);
    statics.str.x = x;
end

xdef    = x+disp;
if ders == 1
    dp_xdef = dp_disp;
    if morphflag == 1
        if isfield(statics.morph,'span')
            dx_xdef2 = eye(length(statics.morph.span.xyz(:,end)));
            for i=1:Ne+1
                dx_xdef(1:3,i,:) = dx_xdef2((1:3)+3*(i-1),:);
            end
        end
    end
end

%% Beam element data
ds     = xdef(:,2:end)-xdef(:,1:end-1);

if ders == 1
    dp_ds  = dp_xdef(:,2:end,:)-dp_xdef(:,1:end-1,:);
    if morphflag == 1
        if isfield(statics.morph,'span')
            dx_ds = dx_xdef(:,2:end,:)-dx_xdef(:,1:end-1,:);
        end
    end
end

l      = (ds(1,:).^2+ds(2,:).^2+ds(3,:).^2).^.5;

%% Calculation of e1, e2, e3 and cs
if morphflag == 1 && ders == 1
    if isfield(statics.morph,'shear')
        pde2dpsi = zeros(3*Ne,Ne);
        pde3dpsi = zeros(3*Ne,Ne);
    end
    if isfield(statics.morph,'fold')
        pde2dtheta = zeros(3*Ne,Ne);
        pde3dtheta = zeros(3*Ne,Ne);
    end
end
for i=1:Ne
    dsi    = ds(:,i);
    
    if ders == 1
        dp_dsi = squeeze(dp_ds(:,i,:));
        if morphflag == 1
            if isfield(statics.morph,'span')
                dx_dsi = squeeze(dx_ds(:,i,:));
            end
        end
    end
    
    li     = l(i);

    tg1          = rot(:,i);
    if ders == 1
        dp_tg1       = squeeze(dp_rot(:,i,:));
        [Rg1,dp_Rg1] = expon(tg1,dp_tg1);
    else
        [Rg1] = expon(tg1);
    end
    tg2          = rot(:,i+1);
    if ders == 1
        dp_tg2       = squeeze(dp_rot(:,i+1,:));
        [Rg2,dp_Rg2] = expon(tg2,dp_tg2);
    else
        [Rg2] = expon(tg2);
    end

%     %Calculate chord directions
%     cv1         = constant.str.c0((i-1)*3+(1:3));
%     cv2         = constant.str.c0((i-1)*3+(4:6));
%     c1s         = Rg1*cv1;
%     c2s         = Rg2*cv2;
    
    if ders == 1
%         dp_c1s(1,:) = dp_Rg1(1,1,:)*cv1(1)+dp_Rg1(1,2,:)*cv1(2)+dp_Rg1(1,3,:)*cv1(3);
%         dp_c1s(2,:) = dp_Rg1(2,1,:)*cv1(1)+dp_Rg1(2,2,:)*cv1(2)+dp_Rg1(2,3,:)*cv1(3);
%         dp_c1s(3,:) = dp_Rg1(3,1,:)*cv1(1)+dp_Rg1(3,2,:)*cv1(2)+dp_Rg1(3,3,:)*cv1(3);
%         dp_c2s(1,:) = dp_Rg1(1,1,:)*cv2(1)+dp_Rg1(1,2,:)*cv2(2)+dp_Rg1(1,3,:)*cv2(3);
%         dp_c2s(2,:) = dp_Rg1(2,1,:)*cv2(1)+dp_Rg1(2,2,:)*cv2(2)+dp_Rg1(2,3,:)*cv2(3);
%         dp_c2s(3,:) = dp_Rg1(3,1,:)*cv2(1)+dp_Rg1(3,2,:)*cv2(2)+dp_Rg1(3,3,:)*cv2(3);
        % Calculate e1
        dp_li = .5*(2*dsi(1)*dp_dsi(1,:)+2*dsi(2)*dp_dsi(2,:)+2*dsi(3)*dp_dsi(3,:))/li;
    end
    
    e1    = dsi/li;
    if ders == 1
        dp_e1 = 1/li^2*(dp_dsi*li-dsi(:,ones(1,size(dp_li,2))).*[dp_li;dp_li;dp_li]);
    end
    
    if morphflag == 1
        if ders == 1
            if isfield(statics.morph,'span')
                dx_li = .5*(2*dsi(1)*dx_dsi(1,:)+2*dsi(2)*dx_dsi(2,:)+2*dsi(3)*dx_dsi(3,:))/li;
                dx_e1 = 1/li^2*(dx_dsi*li-dsi(:,ones(1,size(dx_li,2))).*[dx_li;dx_li;dx_li]);
            end
        end
        if isfield(statics.morph,'shear')
            if constant.morph.shear.sec(i) == 1
                morph.shear = statics.morph.shear.angle(i,end);
                Rpsi = [cos(morph.shear)   -sin(morph.shear)   0;
                    sin(morph.shear)    cos(morph.shear)   0;
                    0           0        1 ];
                
                Ro = constant.str.R0((i-1)*3+(1:3),1:3)*Rpsi;
            else
                Ro = constant.str.R0((i-1)*3+(1:3),1:3);
            end
        else
            Ro = constant.str.R0((i-1)*3+(1:3),1:3);
        end
        if isfield(statics.morph,'fold')
            if constant.morph.fold.sec(i) == 1
                morph.fold = statics.morph.fold.angle(i,end);
                morph.c0 = constant.str.c0(3*(i-1)+(1:3));
                Rtheta = expon(morph.fold*morph.c0);
                Ro = Rtheta*Ro;
            end
        end
    else
        Ro = constant.str.R0((i-1)*3+(1:3),1:3);
    end
    
    % Calculate e3
    q1    = Rg1*Ro*[0;1;0];
    q2    = Rg2*Ro*[0;1;0];
    q     = (q1+q2)/2;
    
    if ders == 1
        Ro010 = Ro*[0;1;0];
        dp_q1(1,1:Ndof) = dp_Rg1(1,1,1:Ndof)*Ro010(1)+dp_Rg1(1,2,1:Ndof)*Ro010(2)+dp_Rg1(1,3,1:Ndof)*Ro010(3);
        dp_q1(2,1:Ndof) = dp_Rg1(2,1,1:Ndof)*Ro010(1)+dp_Rg1(2,2,1:Ndof)*Ro010(2)+dp_Rg1(2,3,1:Ndof)*Ro010(3);
        dp_q1(3,1:Ndof) = dp_Rg1(3,1,1:Ndof)*Ro010(1)+dp_Rg1(3,2,1:Ndof)*Ro010(2)+dp_Rg1(3,3,1:Ndof)*Ro010(3);
        dp_q2(1,1:Ndof) = dp_Rg2(1,1,1:Ndof)*Ro010(1)+dp_Rg2(1,2,1:Ndof)*Ro010(2)+dp_Rg2(1,3,1:Ndof)*Ro010(3);
        dp_q2(2,1:Ndof) = dp_Rg2(2,1,1:Ndof)*Ro010(1)+dp_Rg2(2,2,1:Ndof)*Ro010(2)+dp_Rg2(2,3,1:Ndof)*Ro010(3);
        dp_q2(3,1:Ndof) = dp_Rg2(3,1,1:Ndof)*Ro010(1)+dp_Rg2(3,2,1:Ndof)*Ro010(2)+dp_Rg2(3,3,1:Ndof)*Ro010(3);
        
        dp_q            = (dp_q1+dp_q2)/2;
    end
    
    e3              = [e1(2)*q(3)-e1(3)*q(2);e1(3)*q(1)-e1(1)*q(3);...
                       e1(1)*q(2)-e1(2)*q(1)];
    if ders == 1
        dp_e3           = [(dp_e1(2,:).*q(3)+e1(2).*dp_q(3,:))-(dp_e1(3,:).*q(2)+e1(3).*dp_q(2,:));...
            (dp_e1(3,:).*q(1)+e1(3).*dp_q(1,:))-(dp_e1(1,:).*q(3)+e1(1).*dp_q(3,:));...
            (dp_e1(1,:).*q(2)+e1(1).*dp_q(2,:))-(dp_e1(2,:).*q(1)+e1(2).*dp_q(1,:))];
        dp_ne3          = 1/2/norm(e3)*(2*e3(1).*dp_e3(1,:)+2*e3(2).*dp_e3(2,:)+2*e3(3).*dp_e3(3,:));
        dp_e3           = 1/norm(e3)^2*[dp_e3(1,:)*norm(e3)-e3(1).*dp_ne3;...
            dp_e3(2,:)*norm(e3)-e3(2).*dp_ne3;...
            dp_e3(3,:)*norm(e3)-e3(3).*dp_ne3];
    end
    e3b             = e3;
    e3              = e3./norm(e3);

    % Calculate e2
    e2              = [e3(2)*e1(3)-e3(3)*e1(2);e3(3)*e1(1)-e3(1)*e1(3);...
                        e3(1)*e1(2)-e3(2)*e1(1)];
    if ders == 1
        dp_e2           = [(dp_e3(2,:)*e1(3)+e3(2)*dp_e1(3,:))-(dp_e3(3,:)*e1(2)+e3(3)*dp_e1(2,:));...
            (dp_e3(3,:)*e1(1)+e3(3)*dp_e1(1,:))-(dp_e3(1,:)*e1(3)+e3(1)*dp_e1(3,:));...
            (dp_e3(1,:)*e1(2)+e3(1)*dp_e1(2,:))-(dp_e3(2,:)*e1(1)+e3(2)*dp_e1(1,:))];
        dp_ne2          = 1/2/norm(e2)*(2*e2(1).*dp_e2(1,:)+2*e2(2).*dp_e2(2,:)+2*e2(3).*dp_e2(3,:));
        dp_e2           = 1/norm(e2)^2*[dp_e2(1,:)*norm(e2)-e2(1).*dp_ne2;...
            dp_e2(2,:)*norm(e2)-e2(2).*dp_ne2;...
            dp_e2(3,:)*norm(e2)-e2(3).*dp_ne2];
    end
    e2b             = e2;
    e2              = e2/norm(e2);
    
    if morphflag == 1 && ders == 1
        if isfield(statics.morph,'span')
            dx_e3           = [(dx_e1(2,:).*q(3))-(dx_e1(3,:).*q(2));...
                (dx_e1(3,:).*q(1))-(dx_e1(1,:).*q(3));...
                (dx_e1(1,:).*q(2))-(dx_e1(2,:).*q(1))];
            dx_ne3          = 1/2/norm(e3b)*(2*e3b(1).*dx_e3(1,:)+2*e3b(2).*dx_e3(2,:)+2*e3b(3).*dx_e3(3,:));
            dx_e3           = 1/norm(e3b)^2*[dx_e3(1,:)*norm(e3b)-e3b(1).*dx_ne3;...
                dx_e3(2,:)*norm(e3b)-e3b(2).*dx_ne3;...
                dx_e3(3,:)*norm(e3b)-e3b(3).*dx_ne3];
            
            % Calculate e2
            dx_e2           = [(dx_e3(2,:)*e1(3)+e3b(2)*dx_e1(3,:))-(dx_e3(3,:)*e1(2)+e3b(3)*dx_e1(2,:));...
                (dx_e3(3,:)*e1(1)+e3b(3)*dx_e1(1,:))-(dx_e3(1,:)*e1(3)+e3b(1)*dx_e1(3,:));...
                (dx_e3(1,:)*e1(2)+e3b(1)*dx_e1(2,:))-(dx_e3(2,:)*e1(1)+e3b(2)*dx_e1(1,:))];
            dx_ne2          = 1/2/norm(e2b)*(2*e2b(1).*dx_e2(1,:)+2*e2b(2).*dx_e2(2,:)+2*e2b(3).*dx_e2(3,:));
            dx_e2           = 1/norm(e2b)^2*[dx_e2(1,:)*norm(e2b)-e2b(1).*dx_ne2;...
                dx_e2(2,:)*norm(e2b)-e2b(2).*dx_ne2;...
                dx_e2(3,:)*norm(e2b)-e2b(3).*dx_ne2];
            
            dx_e1_s((i-1)*3+(1:3),:) = dx_e1;
            dx_e2_s((i-1)*3+(1:3),:) = dx_e2;
            dx_e3_s((i-1)*3+(1:3),:) = dx_e3;
        end
        if isfield(statics.morph,'shear')
            if constant.morph.shear.sec(i) == 1
                
                dpsiRpsi   = [-sin(morph.shear)   -cos(morph.shear)   0;
                    cos(morph.shear)   -sin(morph.shear)   0;
                    0            0       0];
                
                if isfield(statics.morph,'fold') && constant.morph.fold.sec(i) == 1
                    dpsiRo = Rtheta*constant.str.R0((i-1)*3+(1:3),1:3)*dpsiRpsi;
                else
                    dpsiRo = constant.str.R0((i-1)*3+(1:3),1:3)*dpsiRpsi;
                end
                
                dpsiq1 = Rg1*dpsiRo*[0;1;0];
                dpsiq2 = Rg2*dpsiRo*[0;1;0];
                dpsiq = (dpsiq1+dpsiq2)/2;
                
                
                [~,dpsie3b]    = vec(e1,q,zeros(size(e1)),dpsiq);
                dpsine3b = 1/2/norm(e3b)*(2*e3b(1)*dpsie3b(1,:)+2*e3b(2)*dpsie3b(2,:)+2*e3b(3)*dpsie3b(3,:));
                dpsie3 = 1/norm(e3b)^2*(dpsie3b*norm(e3b)-e3b.*dpsine3b);
                
                [~,dpsie2] = vec(e3,e1,dpsie3,zeros(size(dpsie3)));
                pde2dpsi((i-1)*3+(1:3),i) = dpsie2;
                pde3dpsi((i-1)*3+(1:3),i) = dpsie3;
            end
        end
        if isfield(statics.morph,'fold')
            if constant.morph.fold.sec(i) == 1
                
                [~,dthetaRtheta]   = expon(morph.fold*morph.c0,morph.c0);
                
                if isfield(statics.morph,'shear') && constant.morph.shear.sec(i) == 1
                    dthetaRo = dthetaRtheta*constant.str.R0((i-1)*3+(1:3),1:3)*Rpsi;
                else
                    dthetaRo = dthetaRtheta*constant.str.R0((i-1)*3+(1:3),1:3);
                end
                
                dthetaq1 = Rg1*dthetaRo*[0;1;0];
                dthetaq2 = Rg2*dthetaRo*[0;1;0];
                dthetaq = (dthetaq1+dthetaq2)/2;
                
                
                [~,dthetae3b]    = vec(e1,q,zeros(size(e1)),dthetaq);
                dthetane3b = 1/2/norm(e3b)*(2*e3b(1)*dthetae3b(1,:)+2*e3b(2)*dthetae3b(2,:)+2*e3b(3)*dthetae3b(3,:));
                dthetae3 = 1/norm(e3b)^2*(dthetae3b*norm(e3b)-e3b.*dthetane3b);
                
                [~,dthetae2] = vec(e3,e1,dthetae3,zeros(size(dthetae3)));
                pde2dtheta((i-1)*3+(1:3),i) = dthetae2;
                pde3dtheta((i-1)*3+(1:3),i) = dthetae3;
            end
        end
    end
    
    e1_s((i-1)*3+(1:3),1)    = e1;
    e2_s((i-1)*3+(1:3),1)    = e2;
    e3_s((i-1)*3+(1:3),1)    = e3;
    
%     c1s_s((i-1)*3+(1:3),1)    = c1s;
%     c2s_s((i-1)*3+(1:3),1)    = c2s;
    
    if ders == 1
        dp_e1_s((i-1)*3+(1:3),:) = dp_e1;
        dp_e2_s((i-1)*3+(1:3),:) = dp_e2;
        dp_e3_s((i-1)*3+(1:3),:) = dp_e3;
        
%         dp_c1s_s((i-1)*3+(1:3),:) = dp_c1s;
%         dp_c2s_s((i-1)*3+(1:3),:) = dp_c2s;
    end
    
end
% cs_s    = [c1s_s(:,1);c2s_s(end-2:end,1)];

statics.str.e1     = e1_s;
statics.str.e2     = e2_s;
statics.str.e3     = e3_s;
% statics.str.c0     = cs_s;
statics.str.xdef     = reshape(xdef,[],1);

if ders == 1
%     dp_cs_s = [dp_c1s_s;dp_c2s_s(end-2:end,:)];
    statics.sens.dp_e1 = dp_e1_s;
    statics.sens.dp_e2 = dp_e2_s;
    statics.sens.dp_e3 = dp_e3_s;
    
    if morphflag == 1
        if isfield(statics.morph,'shear')
            statics.sens.pdpsi_e2 = pde2dpsi;
            statics.sens.pdpsi_e3 = pde3dpsi;
        end
        if isfield(statics.morph,'fold')
            statics.sens.pdtheta_e2 = pde2dtheta;
            statics.sens.pdtheta_e3 = pde3dtheta;
        end
        if isfield(statics.morph,'span')
            statics.sens.pdx_xdef = dx_xdef2;
            statics.sens.pdx_e1 = dx_e1_s;
            statics.sens.pdx_e2 = dx_e2_s;
            statics.sens.pdx_e3 = dx_e3_s;
        end
    end
    
%     statics.sens.dp_c0 = dp_cs_s;
    statics.sens.dp_xdef = dp_disp2;
end


