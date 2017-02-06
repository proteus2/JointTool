% function [sk,dsk] = skew(x,dx)
function varargout = skew(varargin)

if nargin==1
    flag = 1;
    x = varargin{1};
elseif nargin==2
    flag = 2;
    x  = varargin{1};
    dx = varargin{2};
end

sk = [  0     -x(3)  x(2);
        x(3)   0    -x(1);
       -x(2)   x(1)  0  ];
   
varargout{1} = sk;

if flag==2
    if size(x,2)==size(dx,2)
        dsk(1,2) =-dx(3,:);
        dsk(1,3) = dx(2,:);
        dsk(2,1) = dx(3,:);
        dsk(2,3) =-dx(1,:);
        dsk(3,1) =-dx(2,:);
        dsk(3,2) = dx(1,:);
    else
        dsk(1,2,:) =-dx(3,:);
        dsk(1,3,:) = dx(2,:);
        dsk(2,1,:) = dx(3,:);
        dsk(2,3,:) =-dx(1,:);
        dsk(3,1,:) =-dx(2,:);
        dsk(3,2,:) = dx(1,:);        
    end
    varargout{2} = dsk;
end