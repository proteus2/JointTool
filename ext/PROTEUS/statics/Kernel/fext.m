function statics = fext(constant,statics,sc,ders,trimflag,tailflag,alfa,varargin)

if isfield(constant.general,'engine')
    ENGINE = constant.general.engine;
else
    ENGINE = 0;
end

if length(varargin) == 1
    spanflag = 1;
else
    spanflag = 0;
end
    
dt = 0;

if isfield(constant.fext,'dmagnitudedt') == 1
    dt = 1;
end

% Check whether span extension and gravity are present
dextflag = 0;
if spanflag == 1 && isfield(constant.fext,'dmagnitudedext') == 1
    dextflag = 1;
end

%%% analogous with structure routine
%%% non constant input : p
%%% output : Fext , Kfext

if ENGINE
    Ndof = length(statics.str.p_full);
else
    Ndof = length(statics.str.p); 
end

T=[cos(alfa),0,-sin(alfa);
            0,1,0;
            sin(alfa),0,cos(alfa)];
if trimflag == 1
    dTda = [-sin(alfa),0,-cos(alfa);
            0,0,0;
            cos(alfa),0,-sin(alfa)];
end

Kfext_total=zeros(Ndof);
Fext_total=zeros(Ndof,1);
dFext_totaldsc = sparse(numel(Fext_total),1);

if ders == 1
    dKfext_totaldp = sparse(numel(Kfext_total),Ndof);
    dFext_totaldp = sparse(numel(Fext_total),Ndof);
    
    if tailflag == 1
        dKfext_totaldt = sparse(numel(Kfext_total),length(constant.lam.ID));
        dFext_totaldt = sparse(numel(Fext_total),length(constant.lam.ID));
    end    
end

if spanflag == 1
    dKfext_totaldext = sparse(numel(Kfext_total),constant.str.Ns);
    dFext_totaldext = sparse(numel(Fext_total),constant.str.Ns);
end

if trimflag == 1
    dKfext_totalda = sparse(numel(Kfext_total),1);
    dFext_totalda = sparse(numel(Fext_total),1);
    if ders == 1 && tailflag == 1
        dFext_totaldadt = sparse(numel(Fext_total),length(constant.lam.ID));
    end
end

for i=1:constant.fext.nr_fextfields;
    
    if dextflag == 1 && isempty(constant.fext.dmagnitudedext{i})==0
        dext = 1;
    else
        dext = 0;
    end
    
    for j=1:constant.fext.nr_fext{i};
        
        if constant.fext.alphaflag{i}(j) == 1
            trim = trimflag;
        else
            trim = 0;
        end
        
        dof1=constant.fext.dof1{i}(j,:);
        dof2=constant.fext.dof2{i}(j,:);
        xi=constant.fext.xi{i}(j,:);
        v02=constant.fext.v0{i}(j,:);
        
        Nex=[sc*T*constant.fext.magnitude{i}(j,1:3)']';
        Mex=[sc*T*constant.fext.magnitude{i}(j,4:6)']';
        v0 = (T*v02')';
        
        dNexdsc = [T*constant.fext.magnitude{i}(j,1:3)']';
        dMexdsc = [T*constant.fext.magnitude{i}(j,4:6)']';
        
        if trim == 1
            dNexda = [sc*dTda*constant.fext.magnitude{i}(j,1:3)']';
            dMexda = [sc*dTda*constant.fext.magnitude{i}(j,4:6)']';
            dv0da = (dTda*v02');
        end
        
        if ders == 1 && tailflag == 1
            if dt == 1 && isempty(constant.fext.dmagnitudedt{i})==0
                if trim == 1
                    dNexdadt = sc*dTda*constant.fext.dmagnitudedt{i}(6*(j-1)+(1:3),:);
                    dMexdadt = sc*dTda*constant.fext.dmagnitudedt{i}(6*(j-1)+(4:6),:);
                    dv0dadt = (dTda*constant.fext.dv0dt{i}(3*(j-1)+(1:3),:));
                end
                dNexdt = sc*T*constant.fext.dmagnitudedt{i}(6*(j-1)+(1:3),:);
                dMexdt = sc*T*constant.fext.dmagnitudedt{i}(6*(j-1)+(4:6),:);
                dv0dt = T*constant.fext.dv0dt{i}(3*(j-1)+(1:3),:);
            end
        end
        
        if dext == 1
            dNexdext = sc*T*constant.fext.dmagnitudedext{i}(6*(j-1)+(1:3),:);
            dMexdext = sc*T*constant.fext.dmagnitudedext{i}(6*(j-1)+(4:6),:);
        end
        
        %displacements adjacent nodes
        if ENGINE
            u1     = statics.str.p_full(dof1(1:3)');
            theta1 = statics.str.p_full(dof1(4:6)');
            u2     = statics.str.p_full(dof2(1:3)');
            theta2 = statics.str.p_full(dof2(4:6)');
        else
            u1     =statics.str.p(dof1(1:3)');
            theta1 =statics.str.p(dof1(4:6)');
            u2     =statics.str.p(dof2(1:3)');
            theta2 =statics.str.p(dof2(4:6)');
        end
        
        if ders == 1
           du1dp = sparse(3,Ndof);
           du1dp(:,dof1(1:3))=eye(3);
           dtheta1dp = sparse(3,Ndof);
           dtheta1dp(:,dof1(4:6))=eye(3);
           du2dp = sparse(3,Ndof);
           du2dp(:,dof2(1:3))=eye(3);
           dtheta2dp = sparse(3,Ndof);
           dtheta2dp(:,dof2(4:6))=eye(3);
        end
        
        %displacements artificial node
        theta_a=(1-xi)*theta1+xi*theta2;
        ua=(1-xi)*u1+xi*u2;
        
        if ders == 1
           dtheta_adp =  (1-xi)*dtheta1dp+xi*dtheta2dp;
           duadp =  (1-xi)*du1dp+xi*du2dp;
        end
        
        if ders == 1
            [Ra,dRadp] = expon(theta_a,full(dtheta_adp));
        else
            Ra=expon(theta_a);
        end
        
        O3=zeros(3);
        I3=eye(3);
        
        
        B=[I3*(1-xi),-skewmatrix(Ra*v0')*(1-xi),xi*I3,-skewmatrix(Ra*v0')*xi;
            O3,(1-xi)*I3,O3,xi*I3];
        
        if trim == 1
            dBda = [O3,-skewmatrix(Ra*dv0da)*(1-xi),O3,-skewmatrix(Ra*dv0da)*xi;
                O3,O3,O3,O3];
        end
        
        if ders == 1
            if dt == 1 && isempty(constant.fext.dmagnitudedt{i})==0 && tailflag == 1
                for k=1:size(dv0dt,2)
                    if trim == 1
                        dBdadt(:,k) = reshape([O3,-skewmatrix(Ra*dv0dadt(:,k))*(1-xi),O3,-skewmatrix(Ra*dv0dadt(:,k))*xi;
                            O3,O3,O3,O3]',[],1);
                    end
                    dBdt(:,k) = reshape([O3,-skewmatrix(Ra*dv0dt(:,k))*(1-xi),O3,-skewmatrix(Ra*dv0dt(:,k))*xi;
                        O3,O3,O3,O3]',[],1);
                end
            end
            
            for k=1:Ndof
                dBdp(:,k) = reshape([O3,-skewmatrix(squeeze(dRadp(:,:,k))*v0')*(1-xi),O3,-skewmatrix(squeeze(dRadp(:,:,k))*v0')*xi;
                    O3,O3,O3,O3]',[],1);
            end
        end
        
        if (constant.fext.follower{i}(j)==0)
            
            R=eye(3);
            Km=zeros(12,12);
            if ders == 1
               dRdp = zeros(size(dRadp));
               dKmdp = sparse(numel(Km),Ndof);
            end
            
            if trim == 1
                dKmda = sparse(12,12);
            end
        else
            %%% Note gravity is never a follower force, so Km is
            %%% independent of the thickness and span extension
            R=Ra;
            material=[O3,-skewmatrix(R*Nex')*(1-xi),O3,-skewmatrix(R*Nex')*xi;
                      O3,-skewmatrix(R*Mex')*(1-xi),O3,-skewmatrix(R*Mex')*xi];
            Km=B'*material;
            if ders == 1
                dRdp = dRadp;
                for k=1:Ndof
                    dmaterialdp(:,k) = reshape([O3,-skewmatrix(squeeze(dRdp(:,:,k))*Nex')*(1-xi),O3,-skewmatrix(squeeze(dRdp(:,:,k))*Nex')*xi;
                        O3,-skewmatrix(squeeze(dRdp(:,:,k))*Mex')*(1-xi),O3,-skewmatrix(squeeze(dRdp(:,:,k))*Mex')*xi]',[],1);
                end

                for k=1:Ndof
                        dKmdp(:,k) = reshape((reshape(dBdp(:,k),size(B,2),size(B,1))*material+B'*(reshape(dmaterialdp(:,k),size(material,2),size(material,1)))')',[],1);
                end
                dmaterialdsc=[O3,-skewmatrix(R*dNexdsc')*(1-xi),O3,-skewmatrix(R*dNexdsc')*xi;
                O3,-skewmatrix(R*dMexdsc')*(1-xi),O3,-skewmatrix(R*dMexdsc')*xi];
                dKmdsc = B'*dmaterialdsc;
            end
            
            if trim == 1
                dmaterialda = [O3,-skewmatrix(R*dNexda')*(1-xi),O3,-skewmatrix(R*dNexda')*xi;
                    O3,-skewmatrix(R*dMexda')*(1-xi),O3,-skewmatrix(R*dMexda')*xi];
                dKmda = dBda'*material+B'*dmaterialda;
            end
        end
        
        v=skewmatrix(R*Nex')*skewmatrix(Ra*v0');
            
        dvdsc = skewmatrix(R*dNexdsc')*skewmatrix(Ra*v0');
        
        if ders == 1
            for k=1:Ndof
                dvdp(:,k) = reshape((skewmatrix(squeeze(dRdp(:,:,k))*Nex')*skewmatrix(Ra*v0')+skewmatrix(R*Nex')*skewmatrix(squeeze(dRadp(:,:,k))*v0'))',[],1);
            end
            
            if dt == 1 && isempty(constant.fext.dmagnitudedt{i})==0 && tailflag == 1
                for k=1:size(dv0dt,2)
                    dvdt(:,k) = reshape((skewmatrix(R*dNexdt(:,k))*skewmatrix(Ra*v0')+skewmatrix(R*Nex')*skewmatrix(Ra*dv0dt(:,k)))',[],1);
                end
            end
        end
        
        if dext == 1
            for k=1:size(dNexdext,2)
                dvdext(:,k) = reshape((skewmatrix(R*dNexdext(:,k))*skewmatrix(Ra*v0'))',[],1);
            end
        end
        
        if trim == 1
           dvda = skewmatrix(R*dNexda')*skewmatrix(Ra*v0')+skewmatrix(R*Nex')*skewmatrix(Ra*dv0da); 
        end
        
        Kg=[O3,O3,O3,O3;
            O3,(1-xi)^2*v,O3,(1-xi)*xi*v;
            O3,O3,O3,O3;
            O3,(xi)*(1-xi)*v,O3,(xi)^2*v];
        
        if ders == 1
            for k=1:Ndof
                dKgdp(:,k) = reshape(([O3,O3,O3,O3;
                    O3,(1-xi)^2*reshape(dvdp(:,k),3,3)',O3,(1-xi)*xi*reshape(dvdp(:,k),3,3)';
                    O3,O3,O3,O3;
                    O3,(xi)*(1-xi)*reshape(dvdp(:,k),3,3)',O3,(xi)^2*reshape(dvdp(:,k),3,3)'])',[],1);
            end
            
            if dt == 1 && isempty(constant.fext.dmagnitudedt{i})==0 && tailflag == 1
                for k=1:size(dv0dt,2)
                    dKgdt(:,k) = reshape(([O3,O3,O3,O3;
                        O3,(1-xi)^2*reshape(dvdt(:,k),3,3)',O3,(1-xi)*xi*reshape(dvdt(:,k),3,3)';
                        O3,O3,O3,O3;
                        O3,(xi)*(1-xi)*reshape(dvdt(:,k),3,3)',O3,(xi)^2*reshape(dvdt(:,k),3,3)'])',[],1);
                end
            end
        end
        
        if dext == 1
            for k=1:size(dNexdext,2)
                dKgdext(:,k) = reshape(([O3,O3,O3,O3;
                    O3,(1-xi)^2*reshape(dvdext(:,k),3,3)',O3,(1-xi)*xi*reshape(dvdext(:,k),3,3)';
                    O3,O3,O3,O3;
                    O3,(xi)*(1-xi)*reshape(dvdext(:,k),3,3)',O3,(xi)^2*reshape(dvdext(:,k),3,3)'])',[],1);
            end
        end
        
        if trim == 1
            dKgda = [O3,O3,O3,O3;
                O3,(1-xi)^2*dvda,O3,(1-xi)*xi*dvda;
                O3,O3,O3,O3;
                O3,(xi)*(1-xi)*dvda,O3,(xi)^2*dvda];
        end
        
        fg=B'*[R*Nex';R*Mex'];

        dfgdsc = B'*[R*dNexdsc';R*dMexdsc'];
        if ders == 1
            for k=1:Ndof
                dfgdp(:,k) = reshape(dBdp(:,k),size(B,2),size(B,1))*[R*Nex';R*Mex']+B'*[squeeze(dRdp(:,:,k))*Nex';squeeze(dRdp(:,:,k))*Mex'];
            end
            
            if dt == 1 && isempty(constant.fext.dmagnitudedt{i})==0 && tailflag == 1
                for k=1:size(dv0dt,2)
                    if trim == 1
                        dfgdadt(:,k) = reshape(dBdadt(:,k),size(B,2),size(B,1))*[R*Nex';R*Mex']+dBda'*[R*dNexdt(:,k);R*dMexdt(:,k)]+...
                            B'*[R*dNexdadt(:,k);R*dMexdadt(:,k)]+reshape(dBdt(:,k),size(B,2),size(B,1))*[R*dNexda';R*dMexda'];
                    end
                    dfgdt(:,k) = reshape(dBdt(:,k),size(B,2),size(B,1))*[R*Nex';R*Mex']+B'*[R*dNexdt(:,k);R*dMexdt(:,k)];
                end
            end
        end
        
        if dext == 1
            for k=1:size(dNexdext,2)
                dfgdext(:,k) = B'*[R*dNexdext(:,k);R*dMexdext(:,k)];
            end
        end
            
        if trim == 1
           dfgda = dBda'*[R*Nex';R*Mex']+B'*[R*dNexda';R*dMexda']; 
        end
        %transformation to the new global coordinates
        
        
        if ders == 1
            for k=1:Ndof
                [Dg1,dDg1dpi] = Ts(theta1,full(dtheta1dp(:,k)));
                [Dg2,dDg2dpi] = Ts(theta2,full(dtheta2dp(:,k)));
                dDg1dp(:,:,k) = dDg1dpi;
                dDg2dp(:,:,k) = dDg2dpi;
            end
        else
            Dg1 = Ts(theta1);
            Dg2 = Ts(theta2);
        end
        
        H=[I3 O3  O3 O3
            O3 Dg1 O3 O3
            O3 O3  I3 O3
            O3 O3  O3 Dg2];
        
        if ders == 1
            for k=1:Ndof
                dHdp(:,k) = reshape(([O3,O3,O3,O3;
                    O3,squeeze(dDg1dp(:,:,k)),O3,O3;
                    O3,O3,O3,O3;
                    O3,O3,O3,squeeze(dDg2dp(:,:,k))])',[],1);
            end
        end
        
        ft=H'*fg;
        
        dftdsc = H'*dfgdsc;
        
        if ders == 1
            for k=1:Ndof
                dHdpmat = reshape(dHdp(:,k),12,12)';
                dftdp(:,k) = dHdpmat'*fg + H'*dfgdp(:,k);
            end
            
            if dt == 1 && isempty(constant.fext.dmagnitudedt{i})==0 && tailflag == 1
                dftdt = H'*dfgdt;
                if trim == 1
                    dftdadt = H'*dfgdadt;
                end
            end
        end
        
        if dext == 1 
            dftdext = H'*dfgdext;
        end
        
        if trim == 1
            dftda = H'*dfgda;
        end
        
        F=zeros(Ndof,1);
        F([dof1';dof2'])=F([dof1';dof2'])+ft;
        
        dFdsc = sparse(Ndof,1);
        dFdsc([dof1';dof2']) = dftdsc;
        
        if ders == 1
            dFdp = sparse(Ndof,Ndof);
            dFdp([dof1';dof2'],:) = dftdp;
            
            if dt == 1 && isempty(constant.fext.dmagnitudedt{i})==0 && tailflag == 1
                dFdt = sparse(Ndof,size(dv0dt,2));
                dFdt([dof1';dof2'],:) = dftdt;
                if trim == 1
                    dFdadt = sparse(Ndof,size(dv0dt,2));
                    dFdadt([dof1';dof2'],:) = dftdadt;
                end
            end
        end
        
        if dext == 1
            dFdext = sparse(Ndof,size(dNexdext,2));
            dFdext([dof1';dof2'],:) = dftdext;
        end
        
        if trim == 1
            dFda = sparse(Ndof,1);
            dFda([dof1';dof2']) = dftda;
        end
        
        %Kext
        if ders == 1
            for k=1:Ndof
                [Dk1,dDk1dpi] = dTs(theta1,fg(4:6),full(dtheta1dp(:,k)),dfgdp(4:6,k));
                [Dk2,dDk2dpi] = dTs(theta2,fg(10:12),full(dtheta2dp(:,k)),dfgdp(10:12,k));
                dDk1dp(:,:,k) = dDk1dpi;
                dDk2dp(:,:,k) = dDk2dpi;
            end
            
            if dt == 1 && isempty(constant.fext.dmagnitudedt{i})==0 && tailflag == 1
                for k=1:size(dv0dt,2)
                    [~,dDk1dti] = dTs(theta1,fg(4:6),zeros(3,1),dfgdt(4:6,k));
                    [~,dDk2dti] = dTs(theta2,fg(10:12),zeros(3,1),dfgdt(10:12,k));
                    dDk1dt(:,:,k) = dDk1dti;
                    dDk2dt(:,:,k) = dDk2dti;
                end
            end
        else
            %momenten zijn input voor dTS
            Dk1 = dTs(theta1,fg(4:6));
            Dk2 = dTs(theta2,fg(10:12));
        end
        
        if dext == 1
            for k=1:size(dNexdext,2)
                [~,dDk1dexti] = dTs(theta1,fg(4:6),zeros(3,1),dfgdext(4:6,k));
                [~,dDk2dexti] = dTs(theta2,fg(10:12),zeros(3,1),dfgdext(10:12,k));
                dDk1dext(:,:,k) = dDk1dexti;
                dDk2dext(:,:,k) = dDk2dexti;
            end
        end
        
        if trim == 1
            [~,dDk1da] = dTs(theta1,fg(4:6),zeros(3,1),dfgda(4:6,:));
            [~,dDk2da] = dTs(theta2,fg(10:12),zeros(3,1),dfgda(10:12,:));
        end
        
        Kt = H'*(Kg+Km)*H;
        Kt(4:6,4:6) = Kt(4:6,4:6)+Dk1;
        Kt(10:12,10:12) = Kt(10:12,10:12)+Dk2;
        
        if ders == 1
            for k=1:Ndof
                dHdpmat = reshape(dHdp(:,k),12,12)';
                dKgdpmat = reshape(dKgdp(:,k),12,12)';
                dKmdpmat = reshape(dKmdp(:,k),12,12)';
                
                dKtdpmat = dHdpmat'*(Kg+Km)*H+H'*(dKgdpmat+dKmdpmat)*H+H'*(Kg+Km)*dHdpmat;
                dKtdpmat(4:6,4:6) = dKtdpmat(4:6,4:6)+squeeze(dDk1dp(:,:,k));
                dKtdpmat(10:12,10:12) = dKtdpmat(10:12,10:12)+squeeze(dDk2dp(:,:,k));
                
                dKtdp(:,k) = reshape(dKtdpmat',[],1);
            end
                
            if dt == 1 && isempty(constant.fext.dmagnitudedt{i})==0 && tailflag == 1
                for k=1:size(dv0dt,2)
                    dKgdtmat = reshape(dKgdt(:,k),12,12)';
                    
                    dKtdtmat = H'*(dKgdtmat)*H;
                    dKtdtmat(4:6,4:6) = dKtdtmat(4:6,4:6)+squeeze(dDk1dt(:,:,k));
                    dKtdtmat(10:12,10:12) = dKtdtmat(10:12,10:12)+squeeze(dDk2dt(:,:,k));
                    
                    dKtdt(:,k) = reshape(dKtdtmat',[],1);
                end
            end
        end
        
        if dext == 1
            for k=1:size(dNexdext,2)
                dKgdextmat = reshape(dKgdext(:,k),12,12)';
                
                dKtdextmat = H'*(dKgdextmat)*H;
                dKtdextmat(4:6,4:6) = dKtdextmat(4:6,4:6)+squeeze(dDk1dext(:,:,k));
                dKtdextmat(10:12,10:12) = dKtdextmat(10:12,10:12)+squeeze(dDk2dext(:,:,k));
                
                dKtdext(:,k) = reshape(dKtdextmat',[],1);
            end
        end
            
        if trim == 1
            dKtda = H'*(dKgda+dKmda)*H;
            dKtda(4:6,4:6) = dKtda(4:6,4:6)+dDk1da;
            dKtda(10:12,10:12) = dKtda(10:12,10:12)+dDk2da;
        end
        
        K=zeros(Ndof);
        K([dof1';dof2'],[dof1';dof2'])=K([dof1';dof2'],[dof1';dof2'])+Kt;
        
        if ders == 1
            for k=1:Ndof
                dKdpmat = sparse(Ndof,Ndof);
                dKdpmat([dof1';dof2'],[dof1';dof2']) = reshape(dKtdp(:,k),12,12)';
                dKdp(:,k) = reshape(dKdpmat',[],1);
            end
            
            if dt == 1 && isempty(constant.fext.dmagnitudedt{i})==0 && tailflag == 1
                for k=1:size(dv0dt,2)
                    dKdtmat = sparse(Ndof,Ndof);
                    dKdtmat([dof1';dof2'],[dof1';dof2']) = reshape(dKtdt(:,k),12,12)';
                    dKdt(:,k) = reshape(dKdtmat',[],1);
                end
            end
        end
        
        if dext == 1
            for k=1:size(dNexdext,2)
                dKdextmat = sparse(Ndof,Ndof);
                dKdextmat([dof1';dof2'],[dof1';dof2']) = reshape(dKtdext(:,k),12,12)';
                dKdext(:,k) = reshape(dKdextmat',[],1);
            end
        end
        
        if trim == 1
            dKdamat = sparse(Ndof,Ndof);
            dKdamat([dof1';dof2'],[dof1';dof2']) = dKtda;
            dKda = reshape(dKdamat',[],1);
        end
        
        Kfext{i}(j,:,:)=K;
        Fext{i}(j,:)=F;
        
        Kfext_total=Kfext_total+squeeze(K);
        Fext_total=Fext_total+F;
                    
        dFext_totaldsc = dFext_totaldsc+dFdsc;
            
        if ders == 1
            dFext_totaldp = dFext_totaldp+dFdp;
            if dt == 1 && isempty(constant.fext.dmagnitudedt{i})==0 && tailflag == 1
                dFext_totaldt = dFext_totaldt+dFdt;
            end
            
            dKfext_totaldp = dKfext_totaldp+dKdp;
            if dt == 1 && isempty(constant.fext.dmagnitudedt{i})==0 && tailflag == 1
                dKfext_totaldt = dKfext_totaldt+dKdt;
                if trim == 1
                    dFext_totaldadt = dFext_totaldadt+dFdadt;
                end
            end
        end
        
        if dext == 1
            dFext_totaldext = dFext_totaldext+dFdext;
            dKfext_totaldext = dKfext_totaldext+dKdext;
        end
        
        if trim == 1
            dFext_totalda = dFext_totalda+dFda;
            dKfext_totalda = dKfext_totalda+dKda;
        end
    end
end

statics.str.Kfext=Kfext_total;
statics.str.Fext=Fext_total;
statics.sens.dFextdsc = dFext_totaldsc;

if ders == 1
    statics.sens.dFextdp = dFext_totaldp;
    statics.sens.dKfextdp = dKfext_totaldp;
    
    if tailflag == 1
        statics.sens.dFextdt = dFext_totaldt;
        statics.sens.dKfextdt = dKfext_totaldt;
        if trim == 1
            statics.sens.dFextdadt = dFext_totaldadt;
        end
    end
end

if spanflag == 1
    statics.sens.dFextdext = dFext_totaldext;
    statics.sens.dKfextdext = dKfext_totaldext;
end

if trimflag == 1
    statics.sens.dFextda = dFext_totalda;
    statics.sens.dKfextda = dKfext_totalda;
end

function [skewmatrix]=skewmatrix(vector)
skewmatrix=[0,-vector(3),vector(2);vector(3),0,-vector(1);-vector(2),vector(1),0];





