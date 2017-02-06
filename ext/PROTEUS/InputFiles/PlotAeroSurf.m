function [constant] = PlotAeroSurf(constant,wing_data,xyz,theta,c,xref,xsref,FLAGSPAR)
% =====                                                              ==== %
%                       AerodynSurfCoord [17/07/2015]                     %
%                                                                         %
%  This file calculate the twisted wing surface coord. points
%  it also calculates the wingbox coordinates (similar to WingBoxFileGen) 
%  These coord. are only used for plotting at the moment.
%
%  The Coordinates are stored in
%
% =====                                                              ==== %

global GRAPHFLAG

AeroSurfXYZ    = [];

ninp  = length(constant.inp.Aerofoil.NodeProfilesUpper);
nairf = 20;     % number of point used to represent the WB

wb_le_right = wing_data(:,9);    % in % of chord
wb_le_left  = wb_le_right(2:end);
wb_le_vec   = [flipud(wb_le_left); wb_le_right];

wb_te_right = wing_data(:,10);   % in % of chord
wb_te_left  = wb_te_right(2:end);
wb_te_vec   = [flipud(wb_te_left); wb_te_right];

sparloc_right = wing_data(:,12:end);
sparloc_left  = sparloc_right(2:end,:);
sparloc_vec   = [flipud(sparloc_left); sparloc_right];

for j = 1 : ninp
   
    % load and rotate normalised aerofoil node profile
    R          = [cos(theta(j)) sin(theta(j)); -sin(theta(j)) cos(theta(j))];   % rotation matrix
    x_local    = [constant.inp.Aerofoil.NodeProfilesUpper{j}(:,1);constant.inp.Aerofoil.NodeProfilesLower{j}(:,1)]; 
    z_local    = [constant.inp.Aerofoil.NodeProfilesUpper{j}(:,2);constant.inp.Aerofoil.NodeProfilesLower{j}(:,2)];
    x_local    = x_local - xref(j)+xsref(j);                                % twisted about the quarter chord 
    twisted_xz = R*[x_local' ; z_local'];                                   % Apply rotation
    x_twisted  = twisted_xz(1,:)';                               
    z_twisted  = twisted_xz(2,:)';
         
    x_twisted_denorm = x_twisted*c(j) +  xyz(j,1);      %(m)  de-normalised data
    z_twisted_denorm = z_twisted*c(j) +  xyz(j,3);      %(m) de-normalised data
    
    
    AeroSurfXYZ = [AeroSurfXYZ,[x_twisted_denorm xyz(j,2)*ones(size(x_twisted)) z_twisted_denorm]]; %#ok<AGROW>
    

    % ---
    if 1 && GRAPHFLAG == 1    % Plot aerofoil contour at node location on the 3D wing planform
        figure(1); hold on
        plot3(x_twisted_denorm,xyz(j,2)*ones(size(x_twisted)),z_twisted_denorm,'blue')
    end
    % ---
    
%     wb_le  = wing_data(j,9);    % in % of chord
%     wb_te  = wing_data(j,10);   % in % of chord
    
    wb_le = wb_le_vec(j);
    wb_te = wb_te_vec(j);

    WB_X = linspace(wb_le,wb_te,nairf)' - xref(j)+xsref(j);

    Nupper     = size(constant.inp.Aerofoil.NodeProfilesUpper{j},1);
    WB_Z_Upper = interp1(x_twisted(1:Nupper),z_twisted(1:Nupper),WB_X,'PCHIP');                     % normalised WB upper Z coord.
    WB_Z_Lower = interp1(x_twisted(Nupper+1:end),z_twisted(Nupper+1:end),flipud(WB_X),'PCHIP');     % normalised WB lower Z coord.
    
    wb_xlocal = c(j)*[WB_X;flipud(WB_X)]+xyz(j,1);                          % de-normalised
    wb_zlocal = c(j)*[WB_Z_Upper;WB_Z_Lower]+xyz(j,3);                      % de-normalised
    
    TopSkin_xyz{j} = [wb_xlocal(1:nairf)      xyz(j,2)*ones(nairf,1) wb_zlocal(1:nairf)];       %#ok<AGROW>
    BotSkin_xyz{j} = [wb_xlocal(nairf+1:end)  xyz(j,2)*ones(nairf,1) wb_zlocal(nairf+1:end)];   %#ok<AGROW>
    
    if FLAGSPAR
        
        sparloc = sparloc_vec(j,:);
        nspars  = length(sparloc)-sum(isnan(sparloc));
        for i=1:length(sparloc)
            if ~isnan(sparloc(i))
                Spars_Xloc = c(j)*(sparloc(i)- xref(j)+xsref(j))+xyz(j,1);
                if ~exist('Spars_xyzTop') || length(Spars_xyzTop)<i
                    Spars_xyzTop{i}(1,1:3)     = [Spars_Xloc xyz(j,2) interp1(wb_xlocal(1:nairf),wb_zlocal(1:nairf),Spars_Xloc)];
                    Spars_xyzBot{i}(1,1:3)     = [Spars_Xloc xyz(j,2) interp1(wb_xlocal(nairf+1:end),wb_zlocal(nairf+1:end),Spars_Xloc)];
                else
                    Spars_xyzTop{i}(end+1,1:3) = [Spars_Xloc xyz(j,2) interp1(wb_xlocal(1:nairf),wb_zlocal(1:nairf),Spars_Xloc)];
                    Spars_xyzBot{i}(end+1,1:3) = [Spars_Xloc xyz(j,2) interp1(wb_xlocal(nairf+1:end),wb_zlocal(nairf+1:end),Spars_Xloc)];
                end
            end
        end
    end
    
    if GRAPHFLAG == 1 && j>1     % plots the spars (in fill color)
        figure(1); hold on
        for ii = 2 : length(TopSkin_xyz{j})
            Skin_xyz_Local1 = TopSkin_xyz{j-1};
            Skin_xyz_Local2 = TopSkin_xyz{j};
            
            fill3([Skin_xyz_Local1(ii-1,1) Skin_xyz_Local2(ii-1,1) Skin_xyz_Local2(ii,1) Skin_xyz_Local1(ii,1) Skin_xyz_Local1(ii-1,1)],...
                [Skin_xyz_Local1(ii-1,2) Skin_xyz_Local2(ii-1,2) Skin_xyz_Local2(ii,2) Skin_xyz_Local1(ii,2) Skin_xyz_Local1(ii-1,2)],...
                [Skin_xyz_Local1(ii-1,3) Skin_xyz_Local2(ii-1,3) Skin_xyz_Local2(ii,3) Skin_xyz_Local1(ii,3) Skin_xyz_Local1(ii-1,3)], [1 1 0.5],'facealpha',.2,'EdgeAlpha',.1)
            
            Skin_xyz_Local1 = BotSkin_xyz{j-1};
            Skin_xyz_Local2 = BotSkin_xyz{j};
            fill3([Skin_xyz_Local1(ii-1,1) Skin_xyz_Local2(ii-1,1) Skin_xyz_Local2(ii,1) Skin_xyz_Local1(ii,1) Skin_xyz_Local1(ii-1,1)],...
                [Skin_xyz_Local1(ii-1,2) Skin_xyz_Local2(ii-1,2) Skin_xyz_Local2(ii,2) Skin_xyz_Local1(ii,2) Skin_xyz_Local1(ii-1,2)],...
                [Skin_xyz_Local1(ii-1,3) Skin_xyz_Local2(ii-1,3) Skin_xyz_Local2(ii,3) Skin_xyz_Local1(ii,3) Skin_xyz_Local1(ii-1,3)], [1 1 0.5],'facealpha',.2,'EdgeAlpha',.1)
            
        end        
    end
    % ---
end

AerodySurfaceX  = AeroSurfXYZ(:,1:3:end);
AerodySurfaceY  = AeroSurfXYZ(:,2:3:end);
AerodySurfaceZ  = AeroSurfXYZ(:,3:3:end);
Npoints         = length(AerodySurfaceZ);
Wing_TE_xyz     =  [AerodySurfaceX(end,:)' AerodySurfaceY(end,:)' AerodySurfaceZ(end,:)'];
Wing_LE_xyz     =  [AerodySurfaceX(round(Npoints/2),:)' AerodySurfaceY(round(Npoints/2),:)' AerodySurfaceZ(round(Npoints/2),:)'];


if GRAPHFLAG == 1 % plot the aerodynamic surface
    if FLAGSPAR
        for i=1:length(Spars_xyzTop)
            for j=2:size(Spars_xyzTop{i},1)
                fill3([Spars_xyzTop{i}(j-1,1) Spars_xyzTop{i}(j,1) Spars_xyzBot{i}(j,1) Spars_xyzBot{i}(j-1,1) Spars_xyzTop{i}(j-1,1)],...
                    [Spars_xyzTop{i}(j-1,2) Spars_xyzTop{i}(j,2) Spars_xyzBot{i}(j,2) Spars_xyzBot{i}(j-1,2) Spars_xyzTop{i}(j-1,2)],...
                    [Spars_xyzTop{i}(j-1,3) Spars_xyzTop{i}(j,3) Spars_xyzBot{i}(j,3) Spars_xyzBot{i}(j-1,3) Spars_xyzTop{i}(j-1,3)], [1 0.5 1])
            end
        end
    end
    handlesurf = surf(AeroSurfXYZ(:,1:3:end),AeroSurfXYZ(:,2:3:end),AeroSurfXYZ(:,3:3:end),'EdgeColor','none','FaceColor',[0 0 1]);
    set(handlesurf,'facealpha',(0.1))
    grid on
    axis equal
    view([45,30])
    
    plot3(Wing_LE_xyz(:,1),Wing_LE_xyz(:,2),Wing_LE_xyz(:,3),'black','LineWidth',2)
    plot3(Wing_TE_xyz(:,1),Wing_TE_xyz(:,2),Wing_TE_xyz(:,3),'green','LineWidth',2)

end

constant.Coord3D.AeroSurfXYZ = AeroSurfXYZ;
constant.Coord3D.Wing_LE_xyz = Wing_LE_xyz;
constant.Coord3D.Wing_TE_xyz = Wing_TE_xyz;
constant.Coord3D.TopSkin_xyz = TopSkin_xyz;
constant.Coord3D.BotSkin_xyz = BotSkin_xyz;
if FLAGSPAR
    constant.Coord3D.Spars_xyzTop   = Spars_xyzTop;
    constant.Coord3D.Spars_xyzBot   = Spars_xyzBot;
end

