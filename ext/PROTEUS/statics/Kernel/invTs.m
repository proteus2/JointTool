function varargout = invTs(varargin)

if nargin==1
    flag = 1;
    t = varargin{1};
elseif nargin==2
    flag = 2;
    t  = varargin{1};
    dt = varargin{2};
end

nt = norm(t);

I = eye(3,3);

if nt==0
  De = I;
else
  b  = nt/2;
  a  = (sin(b)-b*cos(b))/(nt^2*sin(b));
  M  = skew(t);
  De = I-1/2*M+a*M*M;
end

varargout{1} = De;

if flag==2
    if nt==0
        [~,dM] = skew(t,dt);
        dDe = -1/2*dM;
        varargout{2} = dDe;
    else
        dnt    = 1/2/nt*(2*t(1)*dt(1)+2*t(2)*dt(2)+2*t(3)*dt(3));
        db     = dnt/2;
        da     = 1/(nt^2*sin(b))^2*((cos(b)*db-db*cos(b)+b*sin(b)*db)*(nt^2*sin(b))-(sin(b)-b*cos(b))*(2*nt*dnt*sin(b)+nt^2*cos(b)*db));
        [~,dM] = skew(t,dt);
        dDe    = -1/2*dM+da*M*M+a*dM*M+a*M*dM;
        varargout{2} = dDe;
    end
end