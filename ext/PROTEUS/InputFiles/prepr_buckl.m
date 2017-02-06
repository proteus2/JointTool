function [buckl] = prepr_buckl(constant,str,lumped,panels)

% Locate the ribs
if isfield(lumped,'type')
    NRibs = find(strcmp(lumped.type,'Ribs'));
else
    NRibs=[];
end

if isempty(NRibs)
    buckl_loc = constant.lam.xyzb(:,2);
else
    buckl_loc = lumped.location{NRibs}(:,2);
end

buckl.yloc = buckl_loc;

for i=1:length(buckl_loc)-1
    % Find element that corresponds to the start of a buckling panel
    nsec1 = find(str.xyz(2:3:end)<=buckl_loc(i),1,'last');
    buckl.nsec1{i} = nsec1;
    
    % Find element that corresponds to the end of a buckling panel
    nsec2 = find(str.xyz(2:3:end)>=buckl_loc(i+1),1,'first')-1;
    buckl.nsec2{i} = nsec2;
    
    
    xyz1 = [interp1(str.xyz(2:3:end),str.xyz(1:3:end),buckl_loc(i)),buckl_loc(i),interp1(str.xyz(2:3:end),str.xyz(3:3:end),buckl_loc(i))]';
    xyz2 = [interp1(str.xyz(2:3:end),str.xyz(1:3:end),buckl_loc(i+1)),buckl_loc(i+1),interp1(str.xyz(2:3:end),str.xyz(3:3:end),buckl_loc(i+1))]';
    
    % Store start and end location of buckling panel on the beam
    buckl.xyz{i} = [xyz1;xyz2];
    
    tic
    if nsec1 == nsec2 % Buckling panel covers a single element
        for ncross=1:size(panels.elmloc{nsec1},2)
            buckl.Npan{i}(1,ncross) = size(panels.elmloc{nsec1}{ncross},1);
            for k=1:size(panels.elmloc{nsec1}{ncross},1)
                % Store local coordinates of a panel
                buckl.yzpanel{i}{ncross}{k} = [panels.yzlocal{nsec1}{ncross}(panels.elmloc{nsec1}{ncross}(k,1),:);panels.yzlocal{nsec1}{ncross}(panels.elmloc{nsec1}{ncross}(k,2),:)];
                % Define panel geometry
                w = norm(buckl.yzpanel{i}{ncross}{k}(2,:)-buckl.yzpanel{i}{ncross}{k}(1,:));
                l = norm(xyz2-xyz1);
                
                x = [0;l;l;0];
                y= [0;0;w;w];
                
                % Store laminate corresponding to each panel
                buckl.lam{i}{ncross}(k,1) = panels.lam{nsec1}{ncross}(k);
                buckl.type{i}{ncross}(k,1) = panels.type{nsec1}{ncross}(k);
                
                % Store cross-sectional elements that belong to the buckling
                % panel
                crosselm = panels.elmnum{nsec1}{ncross}{k};
                buckl.crosselm{i}{ncross}{k} = find(ismember(panels.elmconv{nsec1}{ncross},crosselm));
                
                % Store the transformation for each of the panels
                buckl.stiff{i}{ncross}{k} = RR(x,y,10);
            end
        end
    else % Buckling panel covers multiple elements
       for j=nsec1:nsec2
           nsec = j;
           for ncross=1:size(panels.elmloc{nsec},2)
               if j == nsec1 || length(Npan)<ncross
                   Npan(ncross) = 0;
               end
               buckl.Npan{i}(j-nsec1+1,ncross) = size(panels.elmloc{nsec}{ncross},1);
               for k=1:size(panels.elmloc{nsec}{ncross},1)
                   % Store local coordinates of a panel
                   buckl.yzpanel{i}{ncross}{k+Npan(ncross)} = [panels.yzlocal{nsec}{ncross}(panels.elmloc{nsec}{ncross}(k,1),:);panels.yzlocal{nsec}{ncross}(panels.elmloc{nsec}{ncross}(k,2),:)];
                   % Define panel geometry
                   w = norm(buckl.yzpanel{i}{ncross}{k+Npan(ncross)}(2,:)-buckl.yzpanel{i}{ncross}{k+Npan(ncross)}(1,:));
                   l = norm(xyz2-xyz1);
                   
                   x = [0;l;l;0];
                   y= [0;0;w;w];
                   
                   % Store laminate corresponding to each panel
                   buckl.lam{i}{ncross}(k+Npan(ncross),1) = panels.lam{nsec}{ncross}(k);
                   buckl.type{i}{ncross}(k+Npan(ncross),1) = panels.type{nsec}{ncross}(k);
                   
                   % Store cross-sectional elements that belong to the buckling
                   % panel
                   crosselm = panels.elmnum{nsec}{ncross}{k};
                   buckl.crosselm{i}{ncross}{k+Npan(ncross)} = find(ismember(panels.elmconv{nsec}{ncross},crosselm));
                   
                   % Store the transformation for each of the panels
                   buckl.stiff{i}{ncross}{k+Npan(ncross)} = RR(x,y,10);
               end
               Npan(ncross) = Npan(ncross) + size(panels.elmloc{nsec}{ncross},1);
           end
       end
    end    
end
end

function [stiff] = RR(x,y,n)

 [L,L1,L2,wg,sg] = shapefun(n);
 ng = length(sg);
 
 K11 = zeros(n*(n+1)/2); K12 = K11; K16 = K11; K22 = K11; K26 = K11; K66 = K11;
 Kxx = K11; Kxy = K11; Kyy = K11;
 
 for iy = 1:ng
     for ix = 1:ng
         [xg,Ag,Jig,Hgx,Hgy] = mapping(x,y,sg(ix),sg(iy));
         w = zeros(n*(n+1)/2,1);
         thx = w; thy = w; kxx = w; kxy = w; kyy = w; 
         k = 0;
         for r = 1:n
             for t = 1:r
                k=k+1;
                Pval = L(t,ix)*L(r-t+1,iy);
                Psx = L1(t,ix)*L(r-t+1,iy); Psy = L(t,ix)*L1(r-t+1,iy);
                Psxx = L2(t,ix)*L(r-t+1,iy); Psxy = L1(t,ix)*L1(r-t+1,iy); Psyy = L(t,ix)*L2(r-t+1,iy);
                [w(k),thx(k),thy(k),kxx(k),kxy(k),kyy(k)] = tansformvalues(Pval,Psx,Psy,Psxx,Psxy,Psyy,Jig,Hgx,Hgy);
             end
         end
         f = wg(ix)*wg(iy)*Ag;
         K11 = K11 + kxx*kxx'*f;
         K12 = K12 + (kxx*kyy'+kyy*kxx')*f;
         K16 = K16 + (kxx*kxy'+kxy*kxx')*f;
         K22 = K22 + kyy*kyy'*f;
         K26 = K26 + (kyy*kxy'+kxy*kyy')*f;
         K66 = K66 + kxy*kxy'*f;
         Kxx = Kxx - thx*thx'*f;
         Kxy = Kxy - (thx*thy'+thy*thx')*f;
         Kyy = Kyy - thy*thy'*f;
     end
 end
 
%  K11(abs(K11)<=1e-12)=0;
%  K12(abs(K12)<=1e-12)=0;
%  K16(abs(K16)<=1e-12)=0;
%  K22(abs(K22)<=1e-12)=0;
%  K26(abs(K26)<=1e-12)=0;
%  K66(abs(K66)<=1e-12)=0;
%  Kxx(abs(Kxx)<=1e-12)=0;
%  Kxy(abs(Kxy)<=1e-12)=0;
%  Kyy(abs(Kyy)<=1e-12)=0;
 
 stiff.K11 = K11;
 stiff.K12 = K12;
 stiff.K16 = K16;
 stiff.K22 = K22;
 stiff.K26 = K26;
 stiff.K66 = K66;
 stiff.Kxx = Kxx;
 stiff.Kxy = Kxy;
 stiff.Kyy = Kyy;
 
end
function [xg,Ag,Jig,Hgx,Hgy] = mapping(x,y,sx,sy)

 Ng = [(1-sx)*(1-sy)/4,(1+sx)*(1-sy)/4,(1+sx)*(1+sy)/4,(1-sx)*(1+sy)/4];
 N1g = [-(1-sy)/4,(1-sy)/4,(1+sy)/4,-(1+sy)/4;
        -(1-sx)/4,-(1+sx)/4,(1+sx)/4,(1-sx)/4];            
 xg = Ng*x; yg = Ng*y;
 Jg = [N1g*x N1g*y];
 Hgx = ([1/4,-1/4,1/4,-1/4]*x)*[0 1;1 0];
 Hgy = ([1/4,-1/4,1/4,-1/4]*y)*[0 1;1 0];         
 Ag = Jg(1,1)*Jg(2,2)-Jg(2,1)*Jg(1,2);
 Jig = inv(Jg);

end

function [w,thx,thy,kxx,kxy,kyy] = tansformvalues(Pval,Psx,Psy,Psxx,Psxy,Psyy,Jig,Hx,Hy)

w = Pval;
thx = Jig(1,1)*Psx+Jig(1,2)*Psy;
thy = Jig(2,1)*Psx+Jig(2,2)*Psy;
H = Jig*([Psxx Psxy;Psxy Psyy]-thx*Hx-thy*Hy)*Jig';
kxx = H(1,1);
kyy = H(2,2);
kxy = H(1,2)+H(2,1);

end

function [L,L1,L2,w,x] = shapefun(n)

if n < 1
    L = []; L1 = []; L2 = []; w = []; x = [];
    return
end


% Gauss weights and abscissas
T = diag([n:-1:1]./sqrt(2*[n+1:-1:2]-1)./sqrt(2*[n:-1:1]-1),-1); T = T + T';
[V,xg] = eig(T);
w = V(end,:).^2;
x = diag(xg)';
w = 2*w/sum(w);
% Lobatto shape functions and first and second derivatives
L = zeros(n,length(x)); L1 = L; L2 = L;

L2(1,:) = 1;
L1(1,:) = x;
L(1,:) = (x.*L1(1,:)-1)/2;

if n == 1, return; end

L2(2,:) = 3*x;
L1(2,:) = 3*L(1,:)+1;
L(2,:) = (x.*L1(2,:)-L1(1,:))/3;

for k = 3:n
    L1(k,:) = ((2*k-1)*x.*L1(k-1,:)-(k-1)*L1(k-2,:))/k;
    L(k,:) = (x.*L1(k,:)-L1(k-1,:))/(k+1);
    L2(k,:) = k*L1(k-1,:)+x.*L2(k-1,:);
end

end