function plotBeamElements(constant,varargin)

if isempty(varargin)
    delta = zeros(constant.str.Ns+1,6);
    faceAlpha = 1.0;
    edgeAlpha = 0.1;
elseif length(varargin)==1
    delta = varargin{1};
    faceAlpha = 1.0;
    edgeAlpha = 0.1;
elseif length(varargin)==2
    delta = zeros(constant.str.Ns+1);
    faceAlpha = varargin{1};
    edgeAlpha = varargin{2};
elseif length(varargin)==3
    delta = varargin{1};
    faceAlpha = varargin{2};
    edgeAlpha = varargin{3};
else
    error('Too many input arguments')
end

GRAPHFLAG = 1;

ux    = delta(:,1);
uy    = delta(:,2);
uz    = delta(:,3);
rx    = delta(:,4);
ry    = delta(:,5);
rz    = delta(:,6);

global FLAGSPAR

Zoffset = 0; % Offset the top and bottom skin for visualisation

lamvec = zeros(length(constant.lam.ID),1);
for i=1:length(constant.lam.TopID)
    lamvec(constant.lam.TopID{i}) = 1;
    lamvec(constant.lam.BotID{i}) = 2;
end

if FLAGSPAR
    for i=1:length(constant.lam.SparID)
        lamvec(constant.lam.SparID{i}) = 3;
    end
end

for j=1:constant.str.Ns
    
    xyz1 = constant.str.xyz(3*(j-1)+(1:3)) + [ux(j),uy(j),uz(j)]';
    xyz2 = constant.str.xyz(3*(j)+(1:3))   + [ux(j+1),uy(j+1),uz(j+1)]';
    
    % Calculate rotation matrix
    R = expon([rx(j),ry(j),rz(j)]);
    
    crossyz = constant.cross.yzlocal{j}{1};
    crossxyz = constant.str.R0(3*(j-1)+(1:3),:)*[zeros(1,size(crossyz,1));crossyz'];
    
    crossloc = constant.cross.elmloc{j}{1}(1:size(constant.cross.lam{j}{1}(constant.cross.lam{j}{1}~=constant.stringer.lamID(j)),1),:);
    crosslam = constant.cross.lam{j}{1}(constant.cross.lam{j}{1}~=constant.stringer.lamID(j));
    
    for i=1:size(crossloc,1)
        if lamvec(crosslam(i))==1
            Zlam = Zoffset;
        elseif lamvec(crosslam(i))==2
            Zlam = -Zoffset;
        else
            Zlam = 0;
        end

        if GRAPHFLAG 
            hold on
            
            X0  =        [crossxyz(1,crossloc(i,1))+xyz1(1), crossxyz(1,crossloc(i,2))+xyz1(1), crossxyz(1,crossloc(i,2))+xyz2(1), crossxyz(1,crossloc(i,1))+xyz2(1)] - [xyz1(1),xyz1(1),xyz2(1),xyz2(1)];
            Y0  =        [crossxyz(2,crossloc(i,1))+xyz1(2), crossxyz(2,crossloc(i,2))+xyz1(2), crossxyz(2,crossloc(i,2))+xyz2(2), crossxyz(2,crossloc(i,1))+xyz2(2)] - [xyz1(2),xyz1(2),xyz2(2),xyz2(2)];
            Z0  = Zlam + [crossxyz(3,crossloc(i,1))+xyz1(3), crossxyz(3,crossloc(i,2))+xyz1(3), crossxyz(3,crossloc(i,2))+xyz2(3), crossxyz(3,crossloc(i,1))+xyz2(3)] - [xyz1(3),xyz1(3),xyz2(3),xyz2(3)];
            
            XYZ = R*[X0;Y0;Z0];
            
            X   = XYZ(1,:) + [xyz1(1),xyz1(1),xyz2(1),xyz2(1)];
            Y   = XYZ(2,:) + [xyz1(2),xyz1(2),xyz2(2),xyz2(2)];
            Z   = XYZ(3,:) + [xyz1(3),xyz1(3),xyz2(3),xyz2(3)];
            
            fill3(X,Y,Z,1,'facealpha',faceAlpha,'edgealpha',edgeAlpha); 
        
        end
    end
end

colormap('jet')
colormap(flipud(colormap))

% Labels
% axis equal
xlabel('Chord [m]')
ylabel('Span [m]')
zlabel('Height [m]')
view(-50,25)
axis equal

set(gcf, 'Position', [100, 100, 1049, 600]);
movegui(gcf,'east')
end

% Utilities
function [R,dR]=expon(t,dt)

if nargin==1
    flag = 1;
    dR = 0;
else flag = 2;
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
end

end

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
end