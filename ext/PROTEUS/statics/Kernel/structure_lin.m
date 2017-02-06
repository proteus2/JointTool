function statics = structure_lin(constant,statics,ders)

% Input to this file is "constant,statics,ders"

Ne   = constant.str.Ns;
Ndof = constant.str.Ndof;

dp_re   = sparse(12*Ne,Ndof);
if ders == 1
    pdC_re  = sparse(12*Ne,36*Ne);
end

re     = [];
if ders == 1
    pdC_re = [];
end

%%
for i=1:Ne
    dof1 = constant.str.EFT(i,1:6);
    dof2 = constant.str.EFT(i,7:12);
    
    x       = [constant.str.xyz((i-1)*3+(1:3));constant.str.xyz(i*3+(1:3))];
    p       = [statics.str.p((i-1)*6+(1:6));statics.str.p(i*6+(1:6))];
    
    C1 = statics.str.C((i-1)*6+(1:6),:);
    C2 = C1;
    Ro = constant.str.R0((i-1)*3+(1:3),:);
    [rei,pdC1_rei,pdC2_rei,dp_rei] = elem3d_lin(x,p,C1,C2,Ro,ders);

    re = [re;rei];
    if ders==1         
        pdC_re((i-1)*12+(1:12),(1:36)+(i-1)*36) = pdC1_rei+pdC2_rei;
        dp_re((i-1)*12+(1:12),[dof1';dof2'])    = dp_rei;
    end
%     keyboard
end

statics.str.re = re;
if ders == 1
    statics.sens.pdC_re = pdC_re;
    statics.sens.dp_re = dp_re;
end