function [t,dt,d2t,ddbt]=logar(R,dR,d2R,dbR,ddbR)

if nargin==1
    flag = 1;
    dt = 0;
    d2t = 0;
elseif nargin==2
    flag = 2;
    d2t = 0;
elseif nargin==3
    flag = 3;
    ddbt = 0;
elseif nargin==5
    flag = 4;
    d2t = 0;
end

u=[R(3,2)-R(2,3)
   R(1,3)-R(3,1)
   R(2,1)-R(1,2)];

nu=norm(u);

if nu==0
  t=[0;0;0];
else
  t=asin(nu/2)/nu*u;
end

if or(or(flag==2,flag==3),flag==4)
    for i=1:size(dR,3)
        du(:,i) = [dR(3,2,i)-dR(2,3,i)
                   dR(1,3,i)-dR(3,1,i)
                   dR(2,1,i)-dR(1,2,i)];
    end
    dnu = 1/2/norm(u)*(2*u(1)*du(1,:)+2*u(2)*du(2,:)+2*u(3)*du(3,:));
    if nu==0
      dt = .5*du;
    else
      dt = [(1/nu^2*(dnu/sqrt(4-nu^2)*nu-asin(nu/2)*dnu))*u(1);...
          (1/nu^2*(dnu/sqrt(4-nu^2)*nu-asin(nu/2)*dnu))*u(2);...
          (1/nu^2*(dnu/sqrt(4-nu^2)*nu-asin(nu/2)*dnu))*u(3)]+asin(nu/2)/nu.*du;
    end
end
if flag==3
    d2u = [d2R(3,2)-d2R(2,3)
          d2R(1,3)-d2R(3,1)
          d2R(2,1)-d2R(1,2)];
    d2nu = 1/2/norm(u)^2*((2*(du(1)*du(1)+u(1)*d2u(1))+2*(du(2)*du(2)+u(2)*d2u(2))+...
        2*(du(3)*du(3)+u(3)*d2u(3)))*norm(u)-(2*u(1)*du(1)+2*u(2)*du(2)+2*u(3)*du(3))*dnu);
    if nu==0
      d2t = .5*d2u;
    else
      d2t = d2nu*u/sqrt(4-nu^2)/nu+dnu^2*u/(4-nu^2)^(3/2)-2*dnu^2*u/sqrt(4-nu^2)/nu^2+...
          2*dnu*du/sqrt(4-nu^2)/nu+2*asin(nu/2)*u*dnu^2/nu^3-2*asin(nu/2)*du*dnu/nu^2-...
          asin(nu/2)*u*d2nu/nu^2+asin(nu/2)*d2u/nu;
    end
end
if flag==4
    for i=1:size(dR,3)
        ddbu(:,i) = [ddbR(3,2,i)-ddbR(2,3,i)
                     ddbR(1,3,i)-ddbR(3,1,i)
                     ddbR(2,1,i)-ddbR(1,2,i)];
    end
    dbu = [dbR(3,2)-dbR(2,3)
           dbR(1,3)-dbR(3,1)
           dbR(2,1)-dbR(1,2)];
    dbnu = 1/2/norm(u)*(2*u(1)*dbu(1,:)+2*u(2)*dbu(2,:)+2*u(3)*dbu(3,:));
    ddbnu = 1/2/norm(u)^2*((2*(dbu(1)*du(1,:)+u(1)*ddbu(1,:))+2*(dbu(2)*du(2,:)+u(2)*ddbu(2,:))+...
        2*(dbu(3)*du(3,:)+u(3)*ddbu(3,:)))*norm(u)-(2*u(1)*du(1,:)+2*u(2)*du(2,:)+2*u(3)*du(3,:))*dbnu);
    if nu==0
      ddbt = .5*ddbu;
    else
      ddbt = ([ddbnu*u(1);ddbnu*u(2);ddbnu*u(3)])./(sqrt(4-nu^2)*nu)+([dnu.*u(1).*dbnu;dnu.*u(2).*dbnu;dnu.*u(3).*dbnu])./(4-nu^2)^(3/2)-...
          2*([dnu.*u(1).*dbnu;dnu.*u(2).*dbnu;dnu.*u(3).*dbnu])./(sqrt(4-nu^2)*nu^2)+...
          ([dnu*dbu(1);dnu*dbu(2);dnu*dbu(3)])./(sqrt(4-nu^2)*nu)+1/nu^3*([2*asin(nu/2).*u(1).*dnu.*dbnu;...
          2*asin(nu/2).*u(2).*dnu.*dbnu;2*asin(nu/2).*u(3).*dnu.*dbnu])-([asin(.5*nu).*dbu(1).*dnu;...
          asin(.5*nu).*dbu(2).*dnu;asin(.5*nu).*dbu(3).*dnu])./(nu^2)-...
          ([asin(.5*nu).*u(1).*ddbnu;asin(.5*nu).*u(2).*ddbnu;asin(.5*nu).*u(3).*ddbnu])./(nu^2)+...
          (dbnu.*du)./(sqrt(4-nu^2)*nu)-(asin(.5*nu).*du.*dbnu)./(nu^2)+(asin(.5*nu).*ddbu)./(nu);
    end
end