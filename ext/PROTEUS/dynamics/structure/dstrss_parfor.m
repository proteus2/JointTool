function [str] = dstrss_parfor(str,sens,tailflag,varargin)

if length(varargin)==1
    morphflag = 1;
    morph = varargin{1};
else
    morphflag = 0;
end

Minv = str.M\eye(size(str.M));
sizeK = size(str.K,1);

str.Ass = [zeros(sizeK),-Minv*str.K; eye(sizeK),zeros(sizeK)];
str.Bss = [Minv;zeros(sizeK)];

str.nstr = size(str.Ass,1); % Number of structural DOF

if sens == 1
    
    if tailflag == 1
        
        dAssddv = cell(1,size(str.dKddv,2));
        dBssddv = cell(1,size(str.dKddv,2));
        
        dKddv = str.dKddv';
        dMddv = str.dMddv';
        ndv = size(str.dKddv,2);
        Kmat = str.K;
        
        parfor i=1:ndv
            dKij=reshape(dKddv(i,:),sizeK,sizeK)';
            strdAssK=[sparse(sizeK,sizeK),-Minv*dKij;sparse(sizeK,2*sizeK)];
            strdBssK=sparse(2*sizeK,sizeK);
            dMij=reshape(dMddv(i,:),sizeK,sizeK)';
            strdAssM=[zeros(sizeK),Minv*dMij*Minv*Kmat;zeros(sizeK,2*sizeK)];
            strdBssM=[-Minv*dMij*Minv;zeros(sizeK,sizeK)];
            
            dAssddv{i} = reshape((strdAssK+strdAssM)',[],1);
            dBssddv{i} = reshape((strdBssK+strdBssM)',[],1);
        end
        str.dAssddv = cell2mat(dAssddv);
        str.dBssddv = cell2mat(dBssddv);
        clear dKddv
        clear dMddv
        clear dAssddv
        clear dBssddv
    end
    if morphflag == 1
        if morph.twist
            str.dAssdphi = sparse(numel(str.Ass),size(str.dKdphi,2));
            str.dBssdphi = sparse(numel(str.Bss),size(str.dKdphi,2));
            
            for i=1:size(str.dKdphi,2)
                dKij=reshape(str.dKdphi(:,i),size(str.K,2),size(str.K,1))';
                strdAssK=[zeros(sizeK),-Minv*dKij;zeros(sizeK,2*sizeK)];
                strdBssK=[zeros(2*sizeK,sizeK)];
                dMij=reshape(str.dMdphi(:,i),size(str.M,2),size(str.M,1))';
                strdAssM=[zeros(sizeK),Minv*dMij*Minv*str.K;zeros(sizeK,2*sizeK)];
                strdBssM=[-Minv*dMij*Minv;zeros(sizeK,sizeK)];
                
                str.dAssdphi(:,i) = reshape((strdAssK+strdAssM)',[],1);
                str.dBssdphi(:,i) = reshape((strdBssK+strdBssM)',[],1);
            end
        end
        if morph.shear
            str.dAssdpsi = sparse(numel(str.Ass),size(str.dKdpsi,2));
            str.dBssdpsi = sparse(numel(str.Bss),size(str.dKdpsi,2));
            
            for i=1:size(str.dKdpsi,2)
                dKij=reshape(str.dKdpsi(:,i),size(str.K,2),size(str.K,1))';
                strdAssK=[zeros(sizeK),-Minv*dKij;zeros(sizeK,2*sizeK)];
                strdBssK=[zeros(2*sizeK,sizeK)];
                dMij=reshape(str.dMdpsi(:,i),size(str.M,2),size(str.M,1))';
                strdAssM=[zeros(sizeK),Minv*dMij*Minv*str.K;zeros(sizeK,2*sizeK)];
                strdBssM=[-Minv*dMij*Minv;zeros(sizeK,sizeK)];
                
                str.dAssdpsi(:,i) = reshape((strdAssK+strdAssM)',[],1);
                str.dBssdpsi(:,i) = reshape((strdBssK+strdBssM)',[],1);
            end
        end
        if morph.fold
            str.dAssdtheta = sparse(numel(str.Ass),size(str.dKdtheta,2));
            str.dBssdtheta = sparse(numel(str.Bss),size(str.dKdtheta,2));
            
            for i=1:size(str.dKdtheta,2)
                dKij=reshape(str.dKdtheta(:,i),size(str.K,2),size(str.K,1))';
                strdAssK=[zeros(sizeK),-Minv*dKij;zeros(sizeK,2*sizeK)];
                strdBssK=[zeros(2*sizeK,sizeK)];
                dMij=reshape(str.dMdtheta(:,i),size(str.M,2),size(str.M,1))';
                strdAssM=[zeros(sizeK),Minv*dMij*Minv*str.K;zeros(sizeK,2*sizeK)];
                strdBssM=[-Minv*dMij*Minv;zeros(sizeK,sizeK)];
                
                str.dAssdtheta(:,i) = reshape((strdAssK+strdAssM)',[],1);
                str.dBssdtheta(:,i) = reshape((strdBssK+strdBssM)',[],1);
            end
        end
        if morph.span
            str.dAssdext = sparse(numel(str.Ass),size(str.dKdext,2));
            str.dBssdext = sparse(numel(str.Bss),size(str.dKdext,2));
            
            for i=1:size(str.dKdext,2)
                dKij=reshape(str.dKdext(:,i),size(str.K,2),size(str.K,1))';
                strdAssK=[zeros(sizeK),-Minv*dKij;zeros(sizeK,2*sizeK)];
                strdBssK=[zeros(2*sizeK,sizeK)];
                dMij=reshape(str.dMdext(:,i),size(str.M,2),size(str.M,1))';
                strdAssM=[zeros(sizeK),Minv*dMij*Minv*str.K;zeros(sizeK,2*sizeK)];
                strdBssM=[-Minv*dMij*Minv;zeros(sizeK,sizeK)];
                
                str.dAssdext(:,i) = reshape((strdAssK+strdAssM)',[],1);
                str.dBssdext(:,i) = reshape((strdBssK+strdBssM)',[],1);
            end
        end
        if morph.camber
            str.dAssdparam = sparse(numel(str.Ass),size(str.dKdparam,2));
            str.dBssdparam = sparse(numel(str.Bss),size(str.dKdparam,2));
            
            for i=1:size(str.dKdparam,2)
                dKij=reshape(str.dKdparam(:,i),size(str.K,2),size(str.K,1))';
                strdAssK=[zeros(sizeK),-Minv*dKij;zeros(sizeK,2*sizeK)];
                strdBssK=[zeros(2*sizeK,sizeK)];
                dMij=reshape(str.dMdparam(:,i),size(str.M,2),size(str.M,1))';
                strdAssM=[zeros(sizeK),Minv*dMij*Minv*str.K;zeros(sizeK,2*sizeK)];
                strdBssM=[-Minv*dMij*Minv;zeros(sizeK,sizeK)];
                
                str.dAssdparam(:,i) = reshape((strdAssK+strdAssM)',[],1);
                str.dBssdparam(:,i) = reshape((strdBssK+strdBssM)',[],1);
            end
        end
    end
end