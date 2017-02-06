function varargout = Ts(varargin)

if nargin==1
    flag = 1;
    t = varargin{1};
elseif nargin==2
    flag = 2;
    t  = varargin{1};
    dt = varargin{2};
end

nt = norm(t);
I  = eye(3,3);

if nt==0
    Dg = I;
else
    a  = 2*(sin(nt/2)/nt)^2;
    b  = (1-sin(nt)/nt)/nt^2;
    M  = skew(t);
    Dg = I+a*M+b*M*M;
end

varargout{1} = Dg;

if flag==2
    if nt==0
        [~,dM] = skew(t,dt);
        dDg = 1/2*dM;
    else
        dnt      = 1/2/nt*(2*t(1)*dt(1)+2*t(2)*dt(2)+2*t(3)*dt(3));
        da       = 4*(sin(nt/2)/nt)*1/nt^2*(1/2*cos(nt/2)*dnt*nt-sin(nt/2)*dnt);
        db       = 1/nt^4*((-1/nt^2*(cos(nt)*dnt*nt-sin(nt)*dnt))*nt^2-(1-sin(nt)/nt)*2*nt*dnt);
        [~,dM] = skew(t,dt);
        try
            dDg      = da*M+a*dM+db*M*M+b*dM*M+b*M*dM;
        catch err
            keyboard
        end
    end
    varargout{2} = dDg;
end