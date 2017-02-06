function [str,aero,Wing,Nlam] = prepr_stat(xyz,xyzqc,theta,c,xref,Ne,Na,fixnodes,Nlam,lumped,model)

aero = [];
Wing = [];
str  = [];
% xyz  = reshape(xyz',[],1);          % structural node xyz coordinates
c0   = [cos(theta),zeros(length(theta),1),-sin(theta)];            % Initial chord direction
c0   = reshape(c0',[],1);           % c0 is the local vector direction due to twist angle

aero.Na  = Na;

% xyz_str  = reshape(xyz,3,Ne+1)';
wingl = 0;

for i=1:size(xyz,1)-1
    wingl = wingl + norm(xyz(i+1,:)-xyz(i,:));
end

elml = wingl/Ne;

index = find(fixnodes);
xyz_fix = [];
sec_fix = [];
for i=1:length(index)-1
    yfix = linspace(xyz(index(i),2),xyz(index(i+1),2),Nlam.Span(i)+1)';
%-------------------------------------------------------------------------%
% Shift yfix to closed rib location
%-------------------------------------------------------------------------%
%     if ~strcmp(model,'EMB')
%         for j=1:length(yfix)-2
%             if isfield(lumped,'location')
%                 [delta,ind] = min(abs(lumped.location{1}(:,2)-yfix(j+1)));
%             else
%                 delta=[];
%             end
%             if isempty(delta) || delta>=elml/2
%                 yfix(j+1) = yfix(j+1);
%             else
%                 yfix(j+1) = lumped.location{1}(ind,2);
%             end
%         end
%     end
%-------------------------------------------------------------------------%
    xfix = interp1(xyz(:,2),xyz(:,1),yfix);
    zfix = interp1(xyz(:,2),xyz(:,3),yfix);
    xyz_fix = [xyz_fix(1:end-1,:);
               xfix,yfix,zfix];
           
    % Define structural sections for morphing       
    sec_fix(end+(1:length(yfix)-1),1) = i;
end
Nlam.xyz_lam = xyz_fix;

sec_str = [];
for i=1:size(xyz_fix,1)-1
    numelloc = ceil(norm(xyz_fix(i+1,:)-xyz_fix(i,:))/elml);
    if i==1
        xyz_str(1:numelloc+1,:) = [linspace(xyz_fix(i,1),xyz_fix(i+1,1),numelloc+1)',...
                                   linspace(xyz_fix(i,2),xyz_fix(i+1,2),numelloc+1)',...
                                   linspace(xyz_fix(i,3),xyz_fix(i+1,3),numelloc+1)'];
        
        sec_str(end+(1:numelloc),1) = sec_fix(i);
    else
        xyz_str(end:end+numelloc,:) = [linspace(xyz_fix(i,1),xyz_fix(i+1,1),numelloc+1)',...
                                   linspace(xyz_fix(i,2),xyz_fix(i+1,2),numelloc+1)',...
                                   linspace(xyz_fix(i,3),xyz_fix(i+1,3),numelloc+1)'];
        
        sec_str(end+(1:numelloc),1) = sec_fix(i);
    end
end

c0_str   = [interp1(xyz(:,2),c0(1:3:end),xyz_str(:,2)),...
            interp1(xyz(:,2),c0(2:3:end),xyz_str(:,2)),...
            interp1(xyz(:,2),c0(3:3:end),xyz_str(:,2))];        
c_str    = interp1(xyz(:,2),c,xyz_str(:,2));                            % chord (m)
xref_str = interp1(xyz(:,2),xref,xyz_str(:,2)).*c_str;
str.xref = xref_str;
str.c    = c_str;
str.sec  = sec_str;

if 1 % Beam coordinate system (e1,e2,e3)
    e1_str = (xyz_str(2:end,:)-xyz_str(1:end-1,:))';                      % Y local Coord. directions (1st coord of beam axis)
    c0_sec = (c0_str(1:end-1,:)+c0_str(2:end,:))'/2;                      % normalised avg direction of chord axis between 2 nodes (affected by theta only)
    e3_str = cross(e1_str,c0_sec,1);                                      % Z local Coord. directions perpend. to local y and x axis(3rd coord of beam axis)
    e2_str = cross(e3_str,e1_str,1);                                         % X local coord. beam axis
    
    % apply @norm to each col of ei_str, then normalised by dividing, then reshape in column format
    e1_str   = reshape( e1_str ./ repmat( arrayfun(@(idx) norm(e1_str(:,idx)), 1:size(e1_str,2)) ,3,1) ,[],1);   % Normalised, reshaped e1
    e2_str   = reshape( e2_str ./ repmat( arrayfun(@(idx) norm(e2_str(:,idx)), 1:size(e2_str,2)) ,3,1) ,[],1);   % Normalised, reshaped e2
    e3_str   = reshape( e3_str ./ repmat( arrayfun(@(idx) norm(e3_str(:,idx)), 1:size(e3_str,2)) ,3,1) ,[],1);   % Normalised, reshaped e3
    R0 = [e1_str e2_str e3_str]; 
end

xyzqc_str(:,1) = interp1(xyz(:,2),xyzqc(:,1),xyz_str(:,2));
xyzqc_str(:,2) = interp1(xyz(:,2),xyzqc(:,2),xyz_str(:,2));
xyzqc_str(:,3) = interp1(xyz(:,2),xyzqc(:,3),xyz_str(:,2));

xyz_str = reshape(xyz_str',[],1);
c0_str  = reshape(c0_str',[],1);
xyzqc_str = reshape(xyzqc_str',[],1);

str.xyzqc     = xyzqc_str;
str.xyz       = xyz_str;    
str.c0        = c0_str;
str.e1        = e1_str;
str.e2        = e2_str;
str.e3        = e3_str;
str.R0  = R0; % replaces R0

str.Ns = size(xyz_str,1)/3-1;
clear Ne;

% Element freedom table
EFT       = bsxfun(@plus,[0:6:6*(str.Ns-1)]',[1:12]);                       % Elmt Numbering
clampel   = find(str.xyz(2:3:end)==0);                                      % Clamped elements
fxdof     = 6*(clampel-1)+(1:6);
Ndof      = EFT(end,end);                                                   % Number of Dofs
frdof     = setdiff(1:EFT(end,end),fxdof);                                  % Number of Free Dofs ?

str.EFT   = EFT;
str.fxdof = fxdof;
str.Ndof  = Ndof;
str.frdof = frdof;

str.dof.wing.elm  = EFT(:,6)/6;
str.dof.wing.conn = [EFT(:,6)/6, EFT(:,12)/6];
str.dof.wing.nod  = [1; EFT(:,12)/6];
str.dof.wing.coo  = (1:3*size(str.dof.wing.nod))'; 
str.dof.wing.xn   = (1:3:3*size(str.dof.wing.nod))'; 
str.dof.wing.yn   = (2:3:3*size(str.dof.wing.nod))'; 
str.dof.wing.zn   = (3:3:3*size(str.dof.wing.nod))'; 
str.dof.wing.all  = unique(reshape(EFT',[],1));
str.dof.wing.dx   = str.dof.wing.all(1:6:end);
str.dof.wing.dy   = str.dof.wing.all(2:6:end);
str.dof.wing.dz   = str.dof.wing.all(3:6:end);
str.dof.wing.tx   = str.dof.wing.all(4:6:end);
str.dof.wing.ty   = str.dof.wing.all(5:6:end);
str.dof.wing.tz   = str.dof.wing.all(6:6:end);

% Wing.span         = 2*xyz_str(end-1);
% Wing.planformArea = 2*trapz(xyz(:,2),c);
% Wing.AR           = Wing.span^2/Wing.planformArea;
% Wing.CLa          = 2*pi/(1+2/Wing.AR);

Wing.span         = 2*xyz_str(end-1);
Wing.planformArea = trapz(xyz(:,2),c);
Wing.AR           = Wing.span^2/Wing.planformArea;
Wing.CLa          = 2*pi/(1+2/Wing.AR);

% ----
% Display
fprintf('The wing configuration has the following parameters: \n')
fprintf('Wing span [m]      Wing Planform area [m^2]      Aspect ratio [-]     CLa[-]\n')
fprintf('   %3.1f                %3.1f                       %3.1f            %3.1f\n',Wing.span,Wing.planformArea,Wing.AR,Wing.CLa)

