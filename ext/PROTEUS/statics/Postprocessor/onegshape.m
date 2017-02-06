function [statics] = onegshape(constant,statics,ders,tailflag,morphflag)

% Input: 1g twist at structural node locations wrt y-axis
%           structural deformations

% Rotate 1,0,0 with initial twist and then structural twist
% Project on xz-plane
% Compute angle from arctan

% Interpolate twist distribution
theta_ini = interp1(constant.inp.xyz(2:3:end),constant.inp.theta,constant.str.xyz(2:3:end));

% if sens == 1
%     [dthetady,~,dthetady0] = dinterpolation(str.xyz{k}(2:3:end),str.theta_ini{k},xyz0(2:3:end),loc1,loc2);
%     
%     dthetainidxyz  = sparse(numel(theta_ini),numel(str.xyz{k}));
%     
%     dthetainidxyz(:,1:3:end) = dthetady0*dxyz0dxyz(2:3:end,1:3:end);
%     dthetainidxyz(:,2:3:end) = dthetady + dthetady0*dxyz0dxyz(2:3:end,2:3:end);
%     dthetainidxyz(:,3:3:end) = dthetady0*dxyz0dxyz(2:3:end,3:3:end);
% end
cd ../Kernel
for i=1:constant.str.Ns+1
    Rtheta_ini = [cos(theta_ini(i)),0,sin(theta_ini(i));
        0,1,0
        -sin(theta_ini(i)),0,cos(theta_ini(i))];
    
    if ders == 1
        [Rtheta,dRthetadploc] = expon(statics.str.p(6*(i-1)+(4:6)),eye(3));
    else
        Rtheta = expon(statics.str.p(6*(i-1)+(4:6)));
    end

    c0in = [1;0;0];
    
    c0 = Rtheta*Rtheta_ini*c0in;
    
    statics.str.twist(i,1) = atand(-c0(3)/c0(1));
    
    if ders == 1
        dc0dp = zeros(3,length(statics.str.p));
        for k = 1:3
            dc0dp(:,6*(i-1)+3+k) =squeeze(dRthetadploc(:,:,k))*Rtheta_ini*c0in;
        end
        
        dtwistdp = rad2deg(-1/(1+(c0(3)/c0(1))^2)*((c0(1)*dc0dp(3,:)-c0(3)*dc0dp(1,:))/c0(1)^2));
        
        if tailflag == 1
            statics.sens.dC_twist(i,:) = dtwistdp*statics.sens.dC_p;
            statics.sens.dtwistdt(i,:) = dtwistdp*statics.sens.dpdt;
        end
        if morphflag == 1
            if isfield(statics.morph,'twist')
                statics.sens.dtwistdphi(i,:) = dtwistdp*statics.sens.dpdphi;
            end
            if isfield(statics.morph,'camber')
                statics.sens.dtwistdparam(i,:) = dtwistdp*statics.sens.dpdparam;
            end
            if isfield(statics.morph,'fold')
                statics.sens.dtwistdtheta(i,:) = dtwistdp*statics.sens.dpdtheta;
            end
            if isfield(statics.morph,'shear')
                statics.sens.dtwistdpsi(i,:) = dtwistdp*statics.sens.dpdpsi;
            end
            if isfield(statics.morph,'span')
                statics.sens.dtwistdext(i,:) = dtwistdp*statics.sens.dpdext;
            end
        end
    end
end
cd ../Postprocessor
