function varargout = dinvTs(varargin)

if nargin==2
    flag = 1;
    t = varargin{1};
    v = varargin{2};
elseif nargin==4
    flag = 2;
    t  = varargin{1};
    v  = varargin{2};
    dt = varargin{3};
    dv = varargin{4};
end

nt = norm(t);

if nt==0
  vsk = skew(v);
  Dh  = -1/2*vsk;
else
  a   = nt/2;
  eta =(sin(a)-a*cos(a))/(nt^2*sin(a));
  miu = (nt*(nt+sin(nt))-8*sin(a)^2)/(4*nt^4*sin(a)^2);
  I3  = eye(3);
  M   = skew(t);
  M1  = skew(v);
  M2  = t*v'-2*v*t'+(t'*v)*I3;
  M3  = M*M*v*t';
  Dh  = eta*M2+miu*M3-1/2*M1;
end

varargout{1} = Dh;

if flag==2
    if nt==0
      [~,dvsk] = skew(v,dv);
      dDh        = -1/2*dvsk;
    else
      dnt       = 1/2/norm(t)*2*(t(1)*dt(1)+t(2)*dt(2)+t(3)*dt(3));
      da        = dnt/2;
      deta      = 1/(nt^2*sin(a))^2*((cos(a)*da-da*cos(a)+a*sin(a)*da)*...
                  (nt^2*sin(a))-(sin(a)-a*cos(a))*(2*nt*dnt*sin(a)+nt^2*cos(a)*da));
      dmiu      = 1/(4*nt^4*sin(a)^2)^2*( (dnt*(nt+sin(nt))+nt*(dnt+cos(nt)*dnt)-16*sin(a)*cos(a)*da)*...
                  (4*nt^4*sin(a)^2)-(nt*(nt+sin(nt))-8*sin(a)^2)*(16*nt^3*dnt*sin(a)^2+8*nt^4*sin(a)*cos(a)*da) );
      [~,dM]  = skew(t,dt);
      [~,dM1] = skew(v,dv);
      dM2       = dt*v'+t*dv'-2*(dv*t'+v*dt')+(dt'*v+t'*dv)*I3;
      dM3       = dM*M*v*t'+M*dM*v*t'+M*M*dv*t'+M*M*v*dt';
      dDh       = deta*M2+eta*dM2+dmiu*M3+miu*dM3-1/2*dM1;
    end
    varargout{2} = dDh;
end