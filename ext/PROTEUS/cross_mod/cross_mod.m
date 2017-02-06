function [crossmod,stringer] = cross_mod(constant,lampar,stringer,sd,ders,tailflag,morphflag)

% Inputs
% constant.cross.yz{}: cell structure containing matrices with the nodal
% coordinates per cross-section
% constant.cross.elmloc{}: cell structure containing the node numbers that
% define an element
% constant.cross.numel{}: number of elements per cross-section
% lampar.A: matrix containing the A matrix of the laminates
% lampar.D: matrix containing the D matrix of the laminates
% lampar.t: vector containing the laminate thicknesses
% sd: 0 for static analysis, 1 for dynamic analysis
% ders: 0 for no derivatives, 1 for derivatives
% tailflag: flag to identify aeroelastic tailoring, 0: tailoring, 1: no
% tailoring


% Conversion between DECAT_V1 and this code to select the correct sensitivies
conv = [1,4,6,2,5,3];

Ns   = constant.str.Ns;
if ders == 1
    Nlam = length(constant.lam.ID);
    crossmod.dCdA = zeros(36*Ns,6*Nlam);
    crossmod.dCdD = zeros(36*Ns,6*Nlam);
    crossmod.dAdt = zeros(Ns,Nlam);
    crossmod.dQdt = zeros(2*Ns,Nlam);
    crossmod.dIdt = zeros(9*Ns,Nlam);
end

crossmod.Elm2D     = cell(Ns,1);
crossmod.yz_discr  = cell(Ns,1);
crossmod.Prop      = cell(Ns,1);
crossmod.Prop_ABD  = cell(Ns,1);
crossmod.laminates = cell(Ns,1);
crossmod.C         = zeros(6*Ns,6);

Numcross = 0;
for j=1:Ns      % Loop over the number of structural elements
    
    for ncross = 1:length(constant.cross.yzlocal{j})
        Numcross = Numcross+1;
        yz_discr = constant.cross.yzlocal{j}{ncross};
        elmloc_discr = constant.cross.elmloc{j}{ncross};
        lamflag = constant.cross.lam{j}{ncross};
        
        NdiscElmt = size(elmloc_discr,1);
        
        LamID   = constant.lam.ID;                % often faster to extract structure rather than open the structure inside the loop
        
        if sd == 1
            MatID   = constant.lam.matID;
            elm_t   = zeros(NdiscElmt,1);
            elm_rho = zeros(NdiscElmt,1);
            for i=1:NdiscElmt
                try
                    elm_t(i)   = lampar.t(LamID==lamflag(i));
                    elm_rho(i) = constant.mat.rho(MatID(LamID==lamflag(i)));
                catch
                    elm_t(i)   = stringer.t(constant.stringer.lamID==lamflag(i));
                    elm_rho(i) = constant.stringer.mA(constant.stringer.lamID==lamflag(i))/(constant.stringer.h(constant.stringer.lamID==lamflag(i))*stringer.t(constant.stringer.lamID==lamflag(i)));
                    stringer.rho(constant.stringer.lamID==lamflag(i)) = elm_rho(i);
                end
            end
        end
        
        laminates        = sort(unique(lamflag));
        NuniqueLam       = length(laminates);
        crossmod.nlam(j,ncross) = NuniqueLam;
        
        lam_update = lamflag;                                                   % renumber laminates in cross-section from 1 to ...
        for i = 1:NuniqueLam
            lam_update(lam_update == laminates(i)) = i;
        end
        Elm2D = [[1:length(lamflag)]',lam_update,elmloc_discr];
        
        Prop_ABD = cell(NuniqueLam,1);
        Prop  =  zeros(21,NuniqueLam);
        
        % Try catch statement: first try lampar, if lampar non-existent,
        % laminate is part of stringer
        laminate_type = [];
        for i = 1:NuniqueLam
            try
                Prop_ABD{i}   = [lampar.A(3*(find(laminates(i)==LamID)-1)+(1:3),:) zeros(3); zeros(3) lampar.D(3*(find(laminates(i)==LamID)-1)+(1:3),:)];
                Prop(:,i) = symmat2vec([lampar.A(3*(find(laminates(i)==LamID)-1)+(1:3),:) zeros(3); zeros(3) lampar.D(3*(find(laminates(i)==LamID)-1)+(1:3),:)],1);
                laminate_type(i) = 1;
            catch
                Prop_ABD{i}   = [stringer.A(3*(find(laminates(i)==constant.stringer.lamID)-1)+(1:3),:) zeros(3); zeros(3) stringer.D(3*(find(laminates(i)==constant.stringer.lamID)-1)+(1:3),:)];
                Prop(:,i) = symmat2vec([stringer.A(3*(find(laminates(i)==constant.stringer.lamID)-1)+(1:3),:) zeros(3); zeros(3) stringer.D(3*(find(laminates(i)==constant.stringer.lamID)-1)+(1:3),:)],1);
                laminate_type(i) = 2;
            end
        end
        
        laminates_skinspar = laminates(laminate_type==1);
               
        if ders == 1 && tailflag == 1
            [S,StrainF,DS,dStraindmat] = DECAT_v1(Elm2D(:,1:4),yz_discr,Prop,'yes','yes','no');
            crossmod.dCcrossdA{j}{ncross} = zeros(36,6*Nlam);
            crossmod.dCcrossdD{j}{ncross} = zeros(36,6*Nlam);
            for i = 1:length(laminates_skinspar)
                for k = 1:6
                    crossmod.dCdA((j-1)*36+(1:36),(find(laminates_skinspar(i)==LamID)-1)*6+k) = crossmod.dCdA((j-1)*36+(1:36),(find(laminates_skinspar(i)==LamID)-1)*6+k)+reshape(symmat2vec(DS(:,18*(i-1)+conv(k)),-1,6)',[],1);
                    crossmod.dCdD((j-1)*36+(1:36),(find(laminates_skinspar(i)==LamID)-1)*6+k) = crossmod.dCdD((j-1)*36+(1:36),(find(laminates_skinspar(i)==LamID)-1)*6+k)+reshape(symmat2vec(DS(:,18*(i-1)+12+conv(k)),-1,6)',[],1);
                    crossmod.dCcrossdA{j}{ncross}(:,(find(laminates_skinspar(i)==LamID)-1)*6+k) = reshape(symmat2vec(DS(:,18*(i-1)+conv(k)),-1,6)',[],1);
                    crossmod.dCcrossdD{j}{ncross}(:,(find(laminates_skinspar(i)==LamID)-1)*6+k) = reshape(symmat2vec(DS(:,18*(i-1)+12+conv(k)),-1,6)',[],1);
                end
            end
        else
            [S,StrainF] = DECAT_v1(Elm2D(:,1:4),yz_discr,Prop,'yes','no','no');
        end
        crossmod.Elm2D{j}{ncross}    = Elm2D;
        crossmod.Prop{j}{ncross}      = Prop;
        crossmod.Prop_ABD{j}{ncross}  = Prop_ABD;
        crossmod.laminates{j}{ncross} = laminates;
        crossmod.laminate_type{j}{ncross} = laminate_type; % 1: Skin+spars, 2: Stringers
        crossmod.C((j-1)*6+(1:6),1:6) = crossmod.C((j-1)*6+(1:6),1:6)+S;        
        crossmod.Ccross{j}{ncross} = S;
        crossmod.Numcross{j}(ncross,1) = Numcross;
        
        crossmod.StrainF{j}{ncross}=StrainF;
        if ders == 1 && tailflag == 1
            crossmod.dStraindmat{j}{ncross} = dStraindmat;
        end
        
        if sd == 1
            if ders == 1 && tailflag == 1
                [mA,mQ2,mQ3,mI11,mI22,mI23,mI33,dmAdt,dmQ2dt,dmQ3dt,dmI11dt,dmI22dt,dmI23dt,dmI33dt] = dcrossprop(yz_discr,elm_t,elmloc_discr,elm_rho,1);
                if isfield(constant,'stringer')
                    Nlam_string = length(constant.stringer.lamID);
                    dtconv = sparse(lamflag,(1:NdiscElmt), ones(length(lamflag),1),Nlam+Nlam_string,NdiscElmt);
                    dtconv(Nlam+1:end,:) = [];
                else
                    dtconv = sparse(lamflag,(1:NdiscElmt), ones(length(lamflag),1),Nlam,NdiscElmt);
                end
                
                if ncross>1
                    try
                        crossmod.dmAdt(j,:)         = crossmod.dmAdt(j,:)+(dtconv*dmAdt)';
                        crossmod.dmQdt(2*(j-1)+1,:) = crossmod.dmQdt(2*(j-1)+1,:)+(dtconv*dmQ2dt)';
                        crossmod.dmQdt(2*(j-1)+2,:) = crossmod.dmQdt(2*(j-1)+2,:)+(dtconv*dmQ3dt)';
                        crossmod.dmIdt(9*(j-1)+1,:) = crossmod.dmIdt(9*(j-1)+1,:)+(dtconv*dmI11dt)';
                        crossmod.dmIdt(9*(j-1)+5,:) = crossmod.dmIdt(9*(j-1)+5,:)+(dtconv*dmI22dt)';
                        crossmod.dmIdt(9*(j-1)+6,:) = crossmod.dmIdt(9*(j-1)+6,:)+(dtconv*dmI23dt)';
                        crossmod.dmIdt(9*(j-1)+8,:) = crossmod.dmIdt(9*(j-1)+8,:)+(dtconv*dmI23dt)';
                        crossmod.dmIdt(9*(j-1)+9,:) = crossmod.dmIdt(9*(j-1)+9,:)+(dtconv*dmI33dt)';
                    catch err
                        keyboard
                    end
                else
                    crossmod.dmAdt(j,:)         = dtconv*dmAdt;
                    crossmod.dmQdt(2*(j-1)+1,:) = dtconv*dmQ2dt;
                    crossmod.dmQdt(2*(j-1)+2,:) = dtconv*dmQ3dt;
                    crossmod.dmIdt(9*(j-1)+1,:) = dtconv*dmI11dt;
                    crossmod.dmIdt(9*(j-1)+5,:) = dtconv*dmI22dt;
                    crossmod.dmIdt(9*(j-1)+6,:) = dtconv*dmI23dt;
                    crossmod.dmIdt(9*(j-1)+8,:) = dtconv*dmI23dt;
                    crossmod.dmIdt(9*(j-1)+9,:) = dtconv*dmI33dt;
                end
            else
                [mA,mQ2,mQ3,mI11,mI22,mI23,mI33,~,~,~,~,~,~,~] = dcrossprop(yz_discr,elm_t,elmloc_discr,elm_rho,0);
            end
            if ncross>1
                crossmod.mA(j,1)               = crossmod.mA(j,1)+mA;
                crossmod.mI((j-1)*3+(1:3),1:3) = crossmod.mI((j-1)*3+(1:3),1:3)+[mI11,0,0;0,mI22,mI23;0,mI23,mI33];
                crossmod.mQ(j,:)               = crossmod.mQ(j,:)+[mQ2,mQ3];
            else
                crossmod.mA(j,1)               = mA;
                crossmod.mI((j-1)*3+(1:3),1:3) = [mI11,0,0;0,mI22,mI23;0,mI23,mI33];
                crossmod.mQ(j,:)               = [mQ2,mQ3];
            end
        end
        
        if 0
            figure(3+j)
            hold on
            for iplot=1:size(crossmod.Elm2D{j},1)
                plot(crossmod.yz_discr{j}(crossmod.Elm2D{j}(iplot,3:4),1),crossmod.yz_discr{j}(crossmod.Elm2D{j}(iplot,3:4),2),'-xb')
            end
            xtext = (crossmod.yz_discr{j}(crossmod.Elm2D{j}(:,3),1)+crossmod.yz_discr{j}(crossmod.Elm2D{j}(:,4),1))/2;
            ztext = (crossmod.yz_discr{j}(crossmod.Elm2D{j}(:,3),2)+crossmod.yz_discr{j}(crossmod.Elm2D{j}(:,4),2))/2;
            for itext = 1:length(xtext)
                text(xtext(itext),ztext(itext),num2str(crossmod.laminates{j}(crossmod.Elm2D{j}(itext,2))))
            end
        end
    end
end
