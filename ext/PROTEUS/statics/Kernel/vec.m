function varargout = vec(varargin)

if nargin==2
    flag = 1;
    u = varargin{1};
    v = varargin{2};
elseif nargin==4
    flag = 2;
    u  = varargin{1};
    v  = varargin{2};
    du = varargin{3};
    dv = varargin{4};
elseif nargin == 6
    flag = 3;
    u  = varargin{1};
    v  = varargin{2};
    du = varargin{3};
    dv = varargin{4};
    d2u = varargin{5};
    d2v = varargin{6};
elseif nargin == 10
    flag = 4;
    u  = varargin{1};
    v  = varargin{2};
    du = varargin{3};
    dv = varargin{4};
    d2u = varargin{5};
    d2v = varargin{6};
    dbu = varargin{7}; % Derivative wrt second variable
    dbv = varargin{8}; % Derivative wrt second variable
    ddbu = varargin{9}; % Cross derivative wrt first and second variable
    ddbv = varargin{10}; % Cross derivative wrt first and second variable
end

w    = zeros(3,1);
w(1) = u(2)*v(3)-u(3)*v(2);
w(2) = u(3)*v(1)-u(1)*v(3);
w(3) = u(1)*v(2)-u(2)*v(1);

varargout{1} = w;

if flag==2 || flag==3 || flag==4
    dw(1,:) = du(2,:)*v(3)+u(2)*dv(3,:)-du(3,:)*v(2)-u(3)*dv(2,:);
    dw(2,:) = du(3,:)*v(1)+u(3)*dv(1,:)-du(1,:)*v(3)-u(1)*dv(3,:);
    dw(3,:) = du(1,:)*v(2)+u(1)*dv(2,:)-du(2,:)*v(1)-u(2)*dv(1,:);
    varargout{2} = dw;
end

if flag==3
    d2w(1,:) = d2u(2,:)*v(3)+du(2,:)*dv(3,:)+du(2,:)*dv(3,:)+u(2)*d2v(3,:)-...
        d2u(3,:)*v(2)-du(3,:)*dv(2,:)-du(3,:)*dv(2,:)-u(3)*d2v(2,:);
    d2w(2,:) = d2u(3,:)*v(1)+du(3,:)*dv(1,:)+du(3,:)*dv(1,:)+u(3)*d2v(1,:)-...
        d2u(1,:)*v(3)-du(1,:)*dv(3,:)-du(1,:)*dv(3,:)-u(1)*d2v(3,:);
    d2w(3,:) = d2u(1,:)*v(2)+du(1,:)*dv(2,:)+du(1)*dv(2,:)+u(1)*d2v(2,:)-...
        d2u(2,:)*v(1)-du(2,:)*dv(1,:)-du(2,:)*dv(1,:)-u(2)*d2v(1,:);
    varargout{3} = d2w;
end
if flag==4
    ddbw(1,:) = ddbu(2,:)*v(3)+du(2,:)*dbv(3)+dbu(2)*dv(3,:)+u(2)*ddbv(3,:)-...
        ddbu(3,:)*v(2)-du(3,:)*dbv(2)-dbu(3)*dv(2,:)-u(3)*ddbv(2,:);
    ddbw(2,:) = ddbu(3,:)*v(1)+du(3,:)*dbv(1)+dbu(3)*dv(1,:)+u(3)*ddbv(1,:)-...
        ddbu(1,:)*v(3)-du(1,:)*dbv(3)-dbu(1)*dv(3,:)-u(1)*ddbv(3,:);
    ddbw(3,:) = ddbu(1,:)*v(2)+du(1,:)*dbv(2)+dbu(1)*dv(2,:)+u(1)*ddbv(2,:)-...
        ddbu(2,:)*v(1)-du(2,:)*dbv(1)-dbu(2)*dv(1,:)-u(2)*ddbv(1,:);
    varargout{4} = ddbw;
end