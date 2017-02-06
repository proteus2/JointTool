function varargout = dTs(varargin)

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
  Dk = 1/2*skew(v);
else
  e  = t/nt;
  ev = vec(e,v);

  a1 = (cos(nt)-sin(nt)/nt)/nt;
  M1 = v*e'-(e'*v)*(e*e');

  a2 = (1-sin(nt)/nt)/nt;
  M2 = e*v'-2*(e'*v)*(e*e')+(e'*v)*eye(3);

  a3 = sin(nt)/nt-(2*sin(nt/2)/nt)^2;
  M3 = ev*e';

  a4 = 2*(sin(nt/2)/nt)^2;
  M4 = skew(v);

  Dk = a1*M1+a2*M2-a3*M3+a4*M4;
end

varargout{1} = Dk;

if flag==2
    if nt==0
      [~,dDk] = skew(v,dv);
      dDk     = 1/2*dDk;
    else
      dnt       = 1/2/nt*(2*t(1)*dt(1)+2*t(2)*dt(2)+2*t(3)*dt(3));
      de        = 1/nt^2*(dt*nt-t*dnt);
      [~,dev] = vec(e,v,de,dv);

      da1       = 1/nt^2*((-sin(nt)*dnt-1/nt^2*(cos(nt)*dnt*nt-sin(nt)*dnt))*nt-(cos(nt)-sin(nt)/nt)*dnt);
      dM1       = dv*e'+v*de'-(de'*v)*(e*e')-(e'*dv)*(e*e')-(e'*v)*(de*e')-(e'*v)*(e*de');

      da2       = 1/nt^2*((-1/nt^2*(cos(nt)*dnt*nt-sin(nt)*dnt))*nt-(1-sin(nt)/nt)*dnt);
      dM2       = de*v'+e*dv'-2*(de'*v)*(e*e')-2*(e'*dv)*(e*e')-2*(e'*v)*(de*e')-2*(e'*v)*(e*de')+(de'*v+e'*dv)*eye(3);

      da3       = 1/nt^2*(cos(nt)*dnt*nt-sin(nt)*dnt)-(2*(2*sin(nt/2)/nt)*1/nt^2*(cos(nt/2)*dnt*nt-2*sin(nt/2)*dnt));
      dM3       = dev*e'+ev*de';

      da4       = 4*(sin(nt/2)/nt)*1/nt^2*(1/2*cos(nt/2)*dnt*nt-sin(nt/2)*dnt);
      [~,dM4] = skew(v,dv);

      dDk = da1*M1+a1*dM1+da2*M2+a2*dM2-da3*M3-a3*dM3+da4*M4+a4*dM4;
    end
    varargout{2} = dDk;
end