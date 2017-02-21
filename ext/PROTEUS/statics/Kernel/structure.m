function statics = structure(constant,statics,ders,tailflag,morphflag,morphen)

% Input to this file is "constant,statics,ders"

Ne   = constant.str.Ns;
Ndof = constant.str.Ndof;

Ks = zeros(Ndof,Ndof);
Fs = zeros(Ndof,1);

if ders == 1
    dp_Fs   = sparse(Ndof,Ndof);
    dp_Ks   = sparse(Ndof*Ndof,Ndof);
    dp_re   = sparse(12*Ne,Ndof);
    if tailflag == 1
        pdC1_Fs = sparse(Ndof,36*Ne);
        pdC2_Fs = sparse(Ndof,36*Ne);
        pdC1_Ks = sparse(Ndof*Ndof,36*Ne);
        pdC2_Ks = sparse(Ndof*Ndof,36*Ne);
        pdC_re  = sparse(12*Ne,36*Ne);
        dC_Kl = sparse(144*Ne,36*Ne);
    end
end
if morphflag == 1
    if isfield(statics.morph,'twist')
        pdphi_Fs   = sparse(Ndof,Ne);
        pdphi_Ks   = sparse(Ndof*Ndof,Ne);
        pdphi_re   = sparse(12*Ne,Ne);
        if morphen == 1
            pdphi_Mphi = sparse(Ne,Ne);
            dp_Mphi    = sparse(Ne,Ndof);
            if ders == 1 && tailflag == 1
                pdC_Mphi  = sparse(Ne,36*Ne);
            end
            if isfield(statics.morph,'shear')
                pdpsi_Mphi = sparse(Ne,Ne);
            end
            if isfield(statics.morph,'fold')
                pdtheta_Mphi = sparse(Ne,Ne);
            end
            if isfield(statics.morph,'span')
                pdx_Mphi = sparse(Ne,length(constant.str.xyz));
            end
        end
    end
    if isfield(statics.morph,'shear')
        pdpsi_Fs   = sparse(Ndof,Ne);
        pdpsi_Ks   = sparse(Ndof*Ndof,Ne);
        pdpsi_re   = sparse(12*Ne,Ne);
        if morphen == 1
            pdpsi_Mpsi = sparse(Ne,Ne);
            dp_Mpsi    = sparse(Ne,Ndof);
            if ders == 1 && tailflag == 1
                pdC_Mpsi  = sparse(Ne,36*Ne);
            end
            if isfield(statics.morph,'twist')
                pdphi_Mpsi = sparse(Ne,Ne);
            end
            if isfield(statics.morph,'fold')
                pdtheta_Mpsi = sparse(Ne,Ne);
            end
            if isfield(statics.morph,'span')
                pdx_Mpsi = sparse(Ne,length(constant.str.xyz));
            end
        end
    end
    if isfield(statics.morph,'fold')
        pdtheta_Fs   = sparse(Ndof,Ne);
        pdtheta_Ks   = sparse(Ndof*Ndof,Ne);
        pdtheta_re   = sparse(12*Ne,Ne);
        if morphen == 1
            pdtheta_Mtheta = sparse(Ne,Ne);
            dp_Mtheta    = sparse(Ne,Ndof);
            if ders == 1 && tailflag == 1
                pdC_Mtheta  = sparse(Ne,36*Ne);
            end
            if isfield(statics.morph,'twist')
                pdphi_Mtheta = sparse(Ne,Ne);
            end
            if isfield(statics.morph,'shear')
                pdpsi_Mtheta = sparse(Ne,Ne);
            end
            if isfield(statics.morph,'span')
                pdx_Mtheta = sparse(Ne,length(constant.str.xyz));
            end
        end
    end
    if isfield(statics.morph,'span')
        pdx_Fs   = sparse(Ndof,Ndof/2);
        pdx_Ks   = sparse(Ndof*Ndof,Ndof/2);
        pdx_re   = sparse(12*Ne,Ndof/2);
        dx_Kl = sparse(144*Ne,Ndof/2);
    end
end

re     = [];
Kl = [];


%% Process beam elements
for i=1:Ne
    dof1 = constant.str.EFT(i,1:6);
    dof2 = constant.str.EFT(i,7:12);
    
    p       = [statics.str.p((i-1)*6+(1:6));statics.str.p(i*6+(1:6))];
       
    C1 = statics.str.C((i-1)*6+(1:6),:);
    C2 = C1;
    
    Ro = constant.str.R0((i-1)*3+(1:3),:);
    
    if morphflag == 1
        morph = [];
        if morphen == 1
            morph.energy = 1;
        else
            morph.energy = 0;
        end
        if isfield(statics.morph,'twist')
            if constant.morph.twist.sec(i) == 1
                morph.twist = statics.morph.twist.angle(i,end);
            end
        end
        if isfield(statics.morph,'shear')
            if constant.morph.shear.sec(i) == 1
                morph.shear = statics.morph.shear.angle(i,end);
            end
        end
        if isfield(statics.morph,'fold')
            if constant.morph.fold.sec(i) == 1
                morph.fold = statics.morph.fold.angle(i,end);
                morph.c0 = constant.str.c0(3*(i-1)+(1:3));
            end
        end
        if isfield(statics.morph,'span')
            x       = [statics.morph.span.xyz((i-1)*3+(1:3),end);statics.morph.span.xyz(i*3+(1:3),end)];
            morph.span = 1;
        else
            x       = [constant.str.xyz((i-1)*3+(1:3));constant.str.xyz(i*3+(1:3))];
        end
        [ft,Kt,rei,Kli,pdC1_Kt,pdC1_ft,pdC2_Kt,pdC2_ft,dp_Kt,dp_ft,pdC1_rei,pdC2_rei,dp_rei,dC1_Kli,dC2_Kli,morph] = elem3d(x,p,C1,C2,Ro,ders,tailflag,morph);
        if morphen == 1
            if isfield(statics.morph,'twist')
                if constant.morph.twist.sec(i) == 1
                    statics.morph.twist.Mphi(i,end) = morph.Mphi;
                end
            end
            if isfield(statics.morph,'shear')
                if constant.morph.shear.sec(i) == 1
                    statics.morph.shear.Mpsi(i,end) = morph.Mpsi;
                end
            end
            if isfield(statics.morph,'fold')
                if constant.morph.fold.sec(i) == 1
                    statics.morph.fold.Mtheta(i,end) = morph.Mtheta;
                end
            end
        end
    else
        x       = [constant.str.xyz((i-1)*3+(1:3));constant.str.xyz(i*3+(1:3))];
        [ft,Kt,rei,Kli,pdC1_Kt,pdC1_ft,pdC2_Kt,pdC2_ft,dp_Kt,dp_ft,pdC1_rei,pdC2_rei,dp_rei,dC1_Kli,dC2_Kli] = elem3d(x,p,C1,C2,Ro,ders,tailflag);
    end
   
    % Store element stiffness matrix
    statics.str.elmStiffMat{i} = Kt;
    
    Fs([dof1';dof2']) = Fs([dof1';dof2'])+ft;
    Ks([dof1';dof2'],[dof1';dof2']) = Ks([dof1';dof2'],[dof1';dof2'])+Kt;
    re = [re;rei];

    Kl = [Kl;Kli];
    
    %% Sensitivities
    if ders==1
        dp_Fs([dof1';dof2'],[dof1';dof2'])  = dp_Fs([dof1';dof2'],[dof1';dof2'])+dp_ft;
        
        dofs = [dof1';dof2'];
        dofsrow = (dofs-1)*Ndof;
        dofmat = reshape([dofsrow(1)+dofs';
                          dofsrow(2)+dofs';
                          dofsrow(3)+dofs';
                          dofsrow(4)+dofs';
                          dofsrow(5)+dofs';
                          dofsrow(6)+dofs';
                          dofsrow(7)+dofs';
                          dofsrow(8)+dofs';
                          dofsrow(9)+dofs';
                          dofsrow(10)+dofs';
                          dofsrow(11)+dofs';
                          dofsrow(12)+dofs']',[],1);
        dp_Ks(dofmat,[dof1';dof2']) = dp_Ks(dofmat,[dof1';dof2'])+dp_Kt;

        dp_re((i-1)*12+(1:12),[dof1';dof2'])    = dp_rei;
        
        if tailflag == 1
            pdC1_Fs([dof1';dof2'],(1:36)+(i-1)*36) = pdC1_Fs([dof1';dof2'],(1:36)+(i-1)*36)+pdC1_ft;
            pdC2_Fs([dof1';dof2'],(1:36)+(i-1)*36) = pdC2_Fs([dof1';dof2'],(1:36)+(i-1)*36)+pdC2_ft;
            pdC1_Ks(dofmat,(1:36)+(i-1)*36) = pdC1_Ks(dofmat,(1:36)+(i-1)*36)+pdC1_Kt;
            pdC2_Ks(dofmat,(1:36)+(i-1)*36) = pdC2_Ks(dofmat,(1:36)+(i-1)*36)+pdC2_Kt;
            dC_Kl((i-1)*144+(1:144),(1:36)+(i-1)*36) = dC1_Kli+dC2_Kli;
            pdC_re((i-1)*12+(1:12),(1:36)+(i-1)*36) = pdC1_rei+pdC2_rei;
        end
    end
    if morphflag == 1
        dofs = [dof1';dof2'];
        dofsrow = (dofs-1)*Ndof;
        dofmat = reshape([dofsrow(1)+dofs';
            dofsrow(2)+dofs';
            dofsrow(3)+dofs';
            dofsrow(4)+dofs';
            dofsrow(5)+dofs';
            dofsrow(6)+dofs';
            dofsrow(7)+dofs';
            dofsrow(8)+dofs';
            dofsrow(9)+dofs';
            dofsrow(10)+dofs';
            dofsrow(11)+dofs';
            dofsrow(12)+dofs']',[],1);
        if isfield(statics.morph,'twist')
            if constant.morph.twist.sec(i) == 1
                pdphi_re((i-1)*12+(1:12),i) = morph.dphire;
                pdphi_Ks(dofmat,i) = pdphi_Ks(dofmat,i)+morph.dphiKt;
                pdphi_Fs([dof1';dof2'],i)  = pdphi_Fs([dof1';dof2'],i)+morph.dphift;
                if morphen == 1
                    pdphi_Mphi(i,i) = morph.dphiMphi;
                    if isfield(statics.morph,'shear') && constant.morph.shear.sec(i) == 1
                        pdpsi_Mphi(i,i) = morph.dpsiMphi;
                    end
                    if isfield(statics.morph,'fold') && constant.morph.fold.sec(i) == 1
                        pdtheta_Mphi(i,i) = morph.dthetaMphi;
                    end
                    if isfield(statics.morph,'span')
                        pdx_Mphi(i,[(i-1)*3+(1:3),i*3+(1:3)]) = morph.dxMphi;
                    end
                    dp_Mphi(i,[dof1';dof2']) = morph.dpMphi;
                    if ders == 1 && tailflag == 1
                        pdC_Mphi(i,(1:36)+(i-1)*36) = morph.dC1_Mphi+morph.dC2_Mphi;
                    end
                end
            end
        end
        if isfield(statics.morph,'shear')
            if constant.morph.shear.sec(i) == 1
                pdpsi_re((i-1)*12+(1:12),i) = morph.dpsire;
                pdpsi_Ks(dofmat,i) = pdpsi_Ks(dofmat,i)+morph.dpsiKt;
                pdpsi_Fs([dof1';dof2'],i)  = pdpsi_Fs([dof1';dof2'],i)+morph.dpsift;
                if morphen == 1
                    pdpsi_Mpsi(i,i) = morph.dpsiMpsi;
                    if isfield(statics.morph,'twist') && constant.morph.twist.sec(i) == 1
                        pdphi_Mpsi(i,i) = morph.dphiMpsi;
                    end
                    if isfield(statics.morph,'fold') && constant.morph.fold.sec(i) == 1
                        pdtheta_Mpsi(i,i) = morph.dthetaMpsi;
                    end
                    if isfield(statics.morph,'span')
                        pdx_Mpsi(i,[(i-1)*3+(1:3),i*3+(1:3)]) = morph.dxMpsi;
                    end
                    dp_Mpsi(i,[dof1';dof2']) = morph.dpMpsi;
                    if ders == 1 && tailflag == 1
                        pdC_Mpsi(i,(1:36)+(i-1)*36) = morph.dC1_Mpsi+morph.dC2_Mpsi;
                    end
                end
            end
        end
        if isfield(statics.morph,'span')
            pdx_re((i-1)*12+(1:12),[(i-1)*3+(1:3),i*3+(1:3)]) = morph.dx_re;
            pdx_Ks(dofmat,[(i-1)*3+(1:3),i*3+(1:3)]) = pdx_Ks(dofmat,[(i-1)*3+(1:3),i*3+(1:3)])+morph.dx_Kt;
            pdx_Fs([dof1';dof2'],[(i-1)*3+(1:3),i*3+(1:3)])  = pdx_Fs([dof1';dof2'],[(i-1)*3+(1:3),i*3+(1:3)])+morph.dx_ft;
            dx_Kl((i-1)*144+(1:144),[(i-1)*3+(1:3),i*3+(1:3)]) = morph.dx_Kl;
        end
        if isfield(statics.morph,'fold')
            if constant.morph.fold.sec(i) == 1
                pdtheta_re((i-1)*12+(1:12),i) = morph.dthetare;
                pdtheta_Ks(dofmat,i) = pdtheta_Ks(dofmat,i)+morph.dthetaKt;
                pdtheta_Fs([dof1';dof2'],i)  = pdtheta_Fs([dof1';dof2'],i)+morph.dthetaft;
                if morphen == 1
                    pdtheta_Mtheta(i,i) = morph.dthetaMtheta;
                    if isfield(statics.morph,'twist') && constant.morph.twist.sec(i) == 1
                        pdphi_Mtheta(i,i) = morph.dphiMtheta;
                    end
                    if isfield(statics.morph,'shear') && constant.morph.shear.sec(i) == 1
                        pdpsi_Mtheta(i,i) = morph.dpsiMtheta;
                    end
                    if isfield(statics.morph,'span')
                        pdx_Mtheta(i,[(i-1)*3+(1:3),i*3+(1:3)]) = morph.dxMtheta;
                    end
                    dp_Mtheta(i,[dof1';dof2']) = morph.dpMtheta;
                    if ders == 1 && tailflag == 1
                        pdC_Mtheta(i,(1:36)+(i-1)*36) = morph.dC1_Mtheta+morph.dC2_Mtheta;
                    end
                end
            end
        end
    end
end

%% Output
statics.str.Fs = Fs;
statics.str.Ks = Ks;
statics.str.Kl = Kl;
statics.str.re = re;

if ders == 1
    statics.sens.dp_Ks = dp_Ks;
    statics.sens.dp_re = dp_re;
    statics.sens.dp_Fs = dp_Fs;
    
    if tailflag == 1
        statics.sens.pdC_Fs = pdC1_Fs+pdC2_Fs;
        statics.sens.pdC_Ks = pdC1_Ks+pdC2_Ks;
        statics.sens.pdC_re = pdC_re;
        statics.sens.dC_Kl = dC_Kl;
    end
end

if morphflag == 1
    if isfield(statics.morph,'twist')
        statics.sens.pdphi_Ks = pdphi_Ks;
        statics.sens.pdphi_re = pdphi_re;
        statics.sens.pdphi_Fs = pdphi_Fs;
        if morphen == 1
            statics.sens.pdphi_Mphi = pdphi_Mphi;
            statics.sens.dp_Mphi = dp_Mphi;
            if ders == 1 && tailflag == 1
                statics.sens.pdC_Mphi = pdC_Mphi;
            end
            if isfield(statics.morph,'shear')
                statics.sens.pdpsi_Mphi = pdpsi_Mphi;
            end
            if isfield(statics.morph,'fold')
                statics.sens.pdtheta_Mphi = pdtheta_Mphi;
            end
            if isfield(statics.morph,'span')
                statics.sens.pdx_Mphi = pdx_Mphi;
            end
        end
    end
    if isfield(statics.morph,'shear')
        statics.sens.pdpsi_Ks = pdpsi_Ks;
        statics.sens.pdpsi_re = pdpsi_re;
        statics.sens.pdpsi_Fs = pdpsi_Fs;
        if morphen == 1
            statics.sens.pdpsi_Mpsi = pdpsi_Mpsi;
            statics.sens.dp_Mpsi = dp_Mpsi;
            if ders == 1 && tailflag == 1
                statics.sens.pdC_Mpsi = pdC_Mpsi;
            end
            if isfield(statics.morph,'twist')
                statics.sens.pdphi_Mpsi = pdphi_Mpsi;
            end
            if isfield(statics.morph,'fold')
                statics.sens.pdtheta_Mpsi = pdtheta_Mpsi;
            end
            if isfield(statics.morph,'span')
                statics.sens.pdx_Mpsi = pdx_Mpsi;
            end
        end
    end
    if isfield(statics.morph,'span')
        statics.sens.pdx_Ks = pdx_Ks;
        statics.sens.pdx_re = pdx_re;
        statics.sens.pdx_Fs = pdx_Fs;
        statics.sens.dx_Kl = dx_Kl;
    end
    if isfield(statics.morph,'fold')
        statics.sens.pdtheta_Ks = pdtheta_Ks;
        statics.sens.pdtheta_re = pdtheta_re;
        statics.sens.pdtheta_Fs = pdtheta_Fs;
        if morphen == 1
            statics.sens.pdtheta_Mtheta = pdtheta_Mtheta;
            statics.sens.dp_Mtheta = dp_Mtheta;
            if ders == 1 && tailflag == 1
                statics.sens.pdC_Mtheta = pdC_Mtheta;
            end
            if isfield(statics.morph,'shear')
                statics.sens.pdpsi_Mtheta = pdpsi_Mtheta;
            end
            if isfield(statics.morph,'twist')
                statics.sens.pdphi_Mtheta = pdphi_Mtheta;
            end
            if isfield(statics.morph,'span')
                statics.sens.pdx_Mtheta = pdx_Mtheta;
            end
        end
    end
end