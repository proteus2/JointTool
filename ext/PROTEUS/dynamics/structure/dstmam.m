function [str] = dstmam(str,constant,sens,tailflag)

% Local element matrices
for i=1:str.Nel
    C1 = str.elm.C((i-1)*6+(1:6),1:6);
    C2 = str.elm.C((i-1)*6+(1:6),1:6);
    Imat = str.elm.I((i-1)*3+(1:3),1:3);
    Q = str.elm.Q(i,:);
    
    [Ml,Ml_lumped] = locel(str.elm.ell(i),C1,C2,str.elm.A(i),Q,Imat,str.rho,constant,constant.str.R0((i-1)*3+(1:3),:));
    Mlg_structural(i,:,:) = squeeze(str.elm.R(i,:,:))*Ml*squeeze(str.elm.R(i,:,:))';
    
    mlumped=zeros(12,12);
    for j=1:length(constant.lumped.mass);
        for k=1:length(constant.lumped.mass{j});
            
            if (i==constant.lumped.element{j}(k));
                add = squeeze(Ml_lumped{j}(k,:,:));
            else
                add=zeros(12,12);
            end
            %
            mlumped = mlumped + add;
        end
    end
    
    Mlg_lumped(i,:,:)=squeeze(str.elm.R(i,:,:))*mlumped*squeeze(str.elm.R(i,:,:))';
    
end

% Assemble both structural and lumped contributions to mass matrix
Mlg=Mlg_structural+Mlg_lumped;

% Global element matrices
% NOTE: it should work with constant.str.EFT_full
[str.M] = globel(Mlg,str.eft); 

% Clamped boundary condition at the root
frdof = constant.str.frdof;

% Store full matrices
str.Kfull = str.K;
str.Mfull = str.M;

% Store fxdof 
str.Krs = str.K(constant.str.fxdof,frdof); 
str.Mrs = str.M(constant.str.fxdof,frdof);

% Store frdof
str.K = str.K(frdof,frdof); 
str.M = str.M(frdof,frdof);

[V,D] = eig(str.K,str.M);

[lambda,idx] = sort(diag(D));
eigvec       = V(:,idx);

if length(lambda)>=10
    str.eig=sqrt(lambda(1:10));
    str.cycles = str.eig/2/pi;
    str.eigvec = eigvec(:,1:10);
else
    str.eig=sqrt(lambda(1:end));
    str.cycles = str.eig/2/pi;
    str.eigvec = eigvec;
end

if sens == 1
    if tailflag == 1
        str.dMdC = sparse(numel(str.M),numel(str.elm.C));
        for i=1:str.Nel
            C1 = str.elm.C((i-1)*6+(1:6),1:6);
            C2 = str.elm.C((i-1)*6+(1:6),1:6);
            Imat = str.elm.I((i-1)*3+(1:3),1:3);
            Q = str.elm.Q(i,:);
            for j=1:6
                for k=1:6
                    % Only analyse elements
                    % Local element matrices
                    dMlg_str =  zeros(str.Nel,12,12);
                    
                    dC1 = zeros(6);
                    dC2 = zeros(6);
                    dC1(j,k) = 1;
                    [dmlC1,~] = dlocel(str.elm.ell(i),C1,C2,dC1,dC2,str.elm.A(i),Q,Imat,str.rho,constant,constant.str.R0((i-1)*3+(1:3),:),0,0);
                    dMlg_str(i,:,:) = squeeze(str.elm.R(i,:,:))*dmlC1*squeeze(str.elm.R(i,:,:))';
                    
                    dC1 = zeros(6);
                    dC2 = zeros(6);
                    dC2(j,k) = 1;
                    [~,dmlC2] = dlocel(str.elm.ell(i),C1,C2,dC1,dC2,str.elm.A(i),Q,Imat,str.rho,constant,constant.str.R0((i-1)*3+(1:3),:),0,0);
                    dMlg_str(i,:,:) = squeeze(dMlg_str(i,:,:)) + squeeze(str.elm.R(i,:,:))*dmlC2*squeeze(str.elm.R(i,:,:))';
                    
                    %%%% addition of (multiple) external mass
                    dMlg_nonstr =  zeros(str.Nel,12,12);
                    
                    dmlC1_nonstr=zeros(12,12);
                    dmlC2_nonstr=zeros(12,12);
                    
                    for index1=1:length(constant.lumped.mass)
                        for index2=1:length(constant.lumped.mass{index1})
                            
                            if(i==constant.lumped.element{index1}(index2))
                                
                                dC1 = zeros(6);
                                dC2 = zeros(6);
                                dC1(j,k) = 1;
                                [~,~,dmc1,~] = dlocel(str.elm.ell(i),C1,C2,dC1,dC2,str.elm.A(i),Q,Imat,str.rho,constant,constant.str.R0((i-1)*3+(1:3),:),index1,index2);
                                
                                dC1 = zeros(6);
                                dC2 = zeros(6);
                                dC2(j,k) = 1;
                                [~,~,~,dmc2] = dlocel(str.elm.ell(i),C1,C2,dC1,dC2,str.elm.A(i),Q,Imat,str.rho,constant,constant.str.R0((i-1)*3+(1:3),:),index1,index2);
                                
                                add1=squeeze(dmc1);
                                add2=squeeze(dmc2);
                                
                            else
                                add1=zeros(12,12);
                                add2=zeros(12,12);
                            end
                            
                            dmlC1_nonstr=dmlC1_nonstr+add1;
                            dmlC2_nonstr=dmlC2_nonstr+add2;
                            
                        end
                    end
                    dMlg_nonstr(i,:,:) = squeeze(str.elm.R(i,:,:))*dmlC1_nonstr*squeeze(str.elm.R(i,:,:))'+ squeeze(str.elm.R(i,:,:))*dmlC2_nonstr*squeeze(str.elm.R(i,:,:))';
                    
                    % Add contribution of both structural and non-structural
                    % masses
                    dMlg=dMlg_str+dMlg_nonstr;
                    
                    % Global element matrices
                    [dMg] = globel(dMlg,str.eft);
                    
                    % Clamped boundary condition at the root
                    dMgrs = dMg(constant.str.fxdof,constant.str.frdof);
                    
                    dMg=dMg(constant.str.frdof,constant.str.frdof);
                    
                    str.dMrsdC(:,k+(j-1)*6+(i-1)*6*6)=sparse(reshape(dMgrs',[],1));
                    str.dMdC(:,k+(j-1)*6+(i-1)*6*6)=sparse(reshape(dMg',[],1));
                end
            end
            % Only analyse elements
            % Local element matrices
            dmlgdA = zeros(str.Nel,12,12);
            for ii=1:2
                dmlgdQ{ii,1} = zeros(str.Nel,12,12);
            end
            for ii=1:3
                for jj=1:3
                    dmlgdI{ii,jj} = zeros(str.Nel,12,12);
                end
            end
            
            [dmldA,dmldQ,dmldI] = dlocel_m(str.elm.ell(i),C1,C2,str.rho);
            dmlgdA(i,:,:) = squeeze(str.elm.R(i,:,:))*dmldA*squeeze(str.elm.R(i,:,:))';
            for ii=1:2
                dmlgdQ{ii,1}(i,:,:) = squeeze(str.elm.R(i,:,:))*dmldQ{ii,1}*squeeze(str.elm.R(i,:,:))';
            end
            for ii=1:3
                for jj=1:3
                    dmlgdI{ii,jj}(i,:,:) = squeeze(str.elm.R(i,:,:))*dmldI{ii,jj}*squeeze(str.elm.R(i,:,:))';
                end
            end
            
            % Global element matrices
            [dMgdA] = globel(dmlgdA,str.eft);
            for ii=1:2
                [dMgdQ{ii,1}] = globel(dmlgdQ{ii,1},str.eft);
            end
            for ii=1:3
                for jj=1:3
                    [dMgdI{ii,jj}] = globel(dmlgdI{ii,jj},str.eft);
                end
            end
            
            % Clamped boundary condition at the root
            dMgrsdA=dMgdA(constant.str.fxdof,constant.str.frdof);
            dMgdA=dMgdA(constant.str.frdof,constant.str.frdof);
            
            str.dMrsdA(:,i)=sparse(reshape(dMgrsdA',[],1));
            str.dMdA(:,i)=sparse(reshape(dMgdA',[],1));
            
            for ii=1:2
                dMgrsdQ2=dMgdQ{ii,1}(constant.str.fxdof,constant.str.frdof);
                dMgdQ2=dMgdQ{ii,1}(constant.str.frdof,constant.str.frdof);
                str.dMrsdQ(:,2*(i-1)+ii) = sparse(reshape(dMgrsdQ2',[],1));
                str.dMdQ(:,2*(i-1)+ii) = sparse(reshape(dMgdQ2',[],1));
            end
            
            
            for ii=1:3
                for jj=1:3
                    dMgrsdI2=dMgdI{ii,jj}(constant.str.fxdof,constant.str.frdof);
                    dMgdI2=dMgdI{ii,jj}(constant.str.frdof,constant.str.frdof);
                    str.dMrsdI(:,3*3*(i-1)+3*(ii-1)+jj) = sparse(reshape(dMgrsdI2',[],1));
                    str.dMdI(:,3*3*(i-1)+3*(ii-1)+jj) = sparse(reshape(dMgdI2',[],1));
                end
            end
        end
    end
    for i=1:str.Nel
        C1 = str.elm.C((i-1)*6+(1:6),1:6);
        C2 = str.elm.C((i-1)*6+(1:6),1:6);
        Imat = str.elm.I((i-1)*3+(1:3),1:3);
        Q = str.elm.Q(i,:);
        [dm1dL,dmlumpeddL] = dlocel_l(str.elm.ell(i),C1,C2,str.elm.A(i),Q,Imat,str.rho,constant,constant.str.R0((i-1)*3+(1:3),:));
        dMlg = zeros(str.Nel,12,12);
        dMlg_lumped = zeros(str.Nel,12,12);
        dMlg(i,:,:) = squeeze(str.elm.R(i,:,:))*dm1dL*squeeze(str.elm.R(i,:,:))';
        
        dmlumped=zeros(12,12);
        for j=1:length(constant.lumped.mass);
            for k=1:length(constant.lumped.mass{j});
                
                if (i==constant.lumped.element{j}(k));
                    add = squeeze(dmlumpeddL{j}(k,:,:));
                else
                    add=zeros(12,12);
                end
                %
                dmlumped = dmlumped + add;
            end
        end
        
        dMlg_lumped(i,:,:)=squeeze(str.elm.R(i,:,:))*dmlumped*squeeze(str.elm.R(i,:,:))';
        
        dMlg = dMlg+dMlg_lumped;
        
        % Global element matrices
        [dM] = globel(dMlg,str.eft);
        
        % Clamped boundary condition at the root
        dMrs =dM(constant.str.fxdof,constant.str.frdof);
        dM =dM(constant.str.frdof,constant.str.frdof);
        
        str.dMrsdL(:,i) = sparse(reshape(dMrs',[],1));
        str.dMdL(:,i) = sparse(reshape(dM',[],1));
        
        [Ml,Ml_lumped] = locel(str.elm.ell(i),C1,C2,str.elm.A(i),Q,Imat,str.rho,constant,constant.str.R0((i-1)*3+(1:3),:));
        
        for j=1:9
            dMlg = zeros(str.Nel,12,12);
            dR = zeros(12);
            dR(ceil(j/3),j-(ceil(j/3)-1)*3) = 1;
            dR(ceil(j/3)+3,j-(ceil(j/3)-1)*3+3) = 1;
            dR(ceil(j/3)+6,j-(ceil(j/3)-1)*3+6) = 1;
            dR(ceil(j/3)+9,j-(ceil(j/3)-1)*3+9) = 1;
            
            mlumped=zeros(12,12);
            for index1=1:length(constant.lumped.mass);
                for index2=1:length(constant.lumped.mass{index1});
                    
                    if (i==constant.lumped.element{index1}(index2));
                        add = squeeze(Ml_lumped{index1}(index2,:,:));
                    else
                        add=zeros(12,12);
                    end
                    %
                    mlumped = mlumped + add;
                end
            end
            
            dMlg(i,:,:) = dR*(Ml+mlumped)*squeeze(str.elm.R(i,:,:))'+squeeze(str.elm.R(i,:,:))*(Ml+mlumped)*dR';
            
            % Global element matrices
            [dM] = globel(dMlg,str.eft);
            
            % Clamped boundary condition at the root
            dMrs =dM(constant.str.fxdof,constant.str.frdof);
            dM =dM(constant.str.frdof,constant.str.frdof);
            
            str.dMrsdR(:,j+(i-1)*9) = sparse(reshape(dMrs',[],1));
            str.dMdR(:,j+(i-1)*9) = sparse(reshape(dM',[],1));
        end
    end
    
    str.dMrsdxyz0 = str.dMrsdL*str.elm.delldxyz0;
    str.dMdxyz0 = str.dMdL*str.elm.delldxyz0;
    
end
