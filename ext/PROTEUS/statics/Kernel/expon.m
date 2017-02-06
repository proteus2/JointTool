function [R,dR,d2R]=expon(t,dt)

if nargin==1
    flag = 1;
    dR = 0;
else
    flag = 2;
end

I  = eye(3,3);

al  = sqrt(t(1)^2+t(2)^2+t(3)^2);

Rsk = skew(t);

if al<eps
    R  = I;
else
    Rsk2        = Rsk^2;
    R  = I+sin(al)/al*Rsk+0.5*(sin(al/2)/(al/2))^2*Rsk2;
end

if flag==2
    % dRdx
    [Rsk,dRsk]  = skew(t,dt);
    if al<eps
        dR = dRsk;
    else
        dal = 1/2/sqrt(t(1)^2+t(2)^2+t(3)^2)*(2*t(1).*dt(1,:)+2*t(2).*dt(2,:)+2*t(3).*dt(3,:));
        c1 = 1/al^2*(cos(al)*al-sin(al))*dal;
        c2 = 4*sin(al/2)/al^3*(.5*cos(al/2)*al-sin(al/2))*dal;
        for i=1:size(Rsk,1)
            for j=1:size(Rsk,2)
                Rskc1(i,j,:) = c1.*Rsk(i,j);
                Rskc2(i,j,:) = c2.*Rsk2(i,j);
            end    
        end
        for i=1:size(dt,2)
            RdRdRR(:,:,i) = Rsk*dRsk(:,:,i)+dRsk(:,:,i)*Rsk;
        end
        dR = Rskc1+sin(al)/al.*dRsk+Rskc2+2*(sin(al/2)/al)^2.*RdRdRR;
    end
    
    %d2Rdx2
    if al<eps
        d2R = zeros(3,3,size(dt,2),size(dt,2));
    else
        for i=1:size(dt,2)
            for j=1:size(dt,2)
                d2al = -1/al^3*(t(1).*dt(1,j)+t(2).*dt(2,j)+t(3).*dt(3,j))*(t(1).*dt(1,i)+t(2).*dt(2,i)+t(3).*dt(3,i))+...
                    1/al*(dt(1,j)*dt(1,i)+dt(2,j)*dt(2,i)+dt(3,j)*dt(3,i));
                
                c1 = 1/al^2*(cos(al)*al-sin(al));
                c2 = (al^2*(-al*sin(al))-(cos(al)*al-sin(al))*2*al)/al^4;
                c3 = 4*(al^3*1/2*cos(al/2)-sin(al/2)*3*al^2)/al^6*(al/2*cos(al/2)-sin(al/2));
                c4 = 4*sin(al/2)/al^3*(-1/2*al/2*sin(al/2));
                c5 = 4*sin(al/2)/al^3*(.5*cos(al/2)*al-sin(al/2));
                c6 = 0.5*(sin(al/2)/(al/2))^2;
                
                d2R(:,:,i,j) = c1*Rsk*d2al+...
                    c1*dRsk(:,:,j)*dal(:,i)+...
                    c2*Rsk*dal(:,i)*dal(:,j)+...
                    c1*dal(:,j)*dRsk(:,:,i)+...
                    c3*Rsk2*dal(:,i)*dal(:,j)+...
                    c4*Rsk2*dal(:,i)*dal(:,j)+...
                    c5*dal(:,i)*(Rsk*dRsk(:,:,j)+dRsk(:,:,j)*Rsk)+...
                    c5*Rsk2*d2al+...
                    c5*dal(:,j)*(Rsk*dRsk(:,:,i)+dRsk(:,:,i)*Rsk)+...
                    c6*(dRsk(:,:,i)*dRsk(:,:,j)+dRsk(:,:,j)*dRsk(:,:,i));
            end
        end
    end
end