function [lumped] = AddRibs(wing_data,constant,xyz,theta,xref,xsref,c)

global GRAPHFLAG xlsFileName

% =====                                                              ==== %
%                        Ribs Input file [13/04/2015]                     %
%                                                                         %
%  Estimates the Ribs mass based on thickness, material and chord         %
%  Mass is then accounted for as a lumped mass and external forces
% =====                                                              ==== %


rib_data      = xlsread(xlsFileName,'Ribs');
NRibs         = size(rib_data,1);

if NRibs~=0
    
    RibDensity    = rib_data(:,3);
    rib_data(:,3) = [];
    RibThickness  = rib_data(:,2);
    y_ribs        = rib_data(:,1);
    
    
    % ---
    if 1 % linear interpolation of profiles at Ribs locations
        
        ProfileRibs  = cell(NRibs,1);                                                    % at section
        
        for i = 1 : NRibs
            y_sec      = y_ribs(i);
            
            IndexProfile1 = find(y_sec>=xyz(:,2),1,'last');
            IndexProfile2 = find(y_sec<=xyz(:,2),1);
            
            ProfileRibs{i} = ([constant.inp.Aerofoil.NodeProfilesUpper{IndexProfile2};constant.inp.Aerofoil.NodeProfilesLower{IndexProfile2}] - [constant.inp.Aerofoil.NodeProfilesUpper{IndexProfile1};constant.inp.Aerofoil.NodeProfilesLower{IndexProfile1}])...
                /(xyz(IndexProfile2,2)-xyz(IndexProfile1,2)+eps)  * (y_sec-xyz(IndexProfile1,2)) + [constant.inp.Aerofoil.NodeProfilesUpper{IndexProfile1};constant.inp.Aerofoil.NodeProfilesLower{IndexProfile1}];
        end
        
        if 0   % plot Aerofoil profile interpolated (visual check)
            for iaerofoil = 1:size(rib_data,1)
                figure (2)
                hold all
                xlabel('y')
                ylabel('z')
                plot(ProfileRibs{iaerofoil}(:,1),ProfileRibs{iaerofoil}(:,2),'--')
            end
        end
    end
    % ---
    
    % ---
    if 1    % interpolate with xyz
        
        rib_Xcentroid = zeros(NRibs,1);
        rib_Ycentroid = zeros(NRibs,1);
        rib_Zcentroid = zeros(NRibs,1);
        RibMassDist   = zeros(NRibs,1);
        
        wb_le_right = wing_data(:,9);
        wb_le_left  = wb_le_right(2:end);
        wb_le0      = [flipud(wb_le_left); wb_le_right];

        wb_te_right = wing_data(:,10);
        wb_te_left  = wb_te_right(2:end);
        wb_te0      = [flipud(wb_te_left); wb_te_right];

        WB_le          = interp1(xyz(:,2),wb_le0,y_ribs);
        WB_te          = interp1(xyz(:,2),wb_te0,y_ribs);
        
        WB_RibsXCoord0 = cell(NRibs,1);                                                   % wing box X points at node location
        WB_RibsXCoord  = cell(NRibs,1);
        c_ribs         = interp1(xyz(:,2),c,y_ribs); % chord at ribs location
        theta_ribs     = interp1(xyz(:,2),theta,y_ribs);
        xref_ribs   = interp1(xyz(:,2),xref,y_ribs);
        xsref_ribs  = interp1(xyz(:,2),xsref,y_ribs);
        
        
        for ii = 1:NRibs
            N = size(ProfileRibs{ii},1)/2;
            
            WB_RibsXCoord0{ii} = (linspace(WB_le(ii,1),WB_te(ii,1),20)- xref_ribs(ii)+xsref_ribs(ii))*c_ribs(ii)+ interp1(xyz(:,2),xyz(:,1),y_ribs(ii));
            
            % Rotate the airfoil
            R         = eye(2);%[cos(theta_ribs(ii)) sin(theta_ribs(ii)); -sin(theta_ribs(ii)) cos(theta_ribs(ii))];
            
            x_local    =  ProfileRibs{ii}(:,1);                                    % normalised profile at ribs location
            z_local    =  ProfileRibs{ii}(:,2);
            x_local    =  x_local - xref_ribs(ii)+xsref_ribs(ii);                                           % twisted about the quarter chord
            
            twisted_xz = R*[x_local' ; z_local'];
            x_twisted  = twisted_xz(1,:)' * c_ribs(ii) + interp1(xyz(:,2),xyz(:,1),y_ribs(ii)); % Real twisted profile at ribs location (define from x=0 to x=c)
            
            z_twisted  = twisted_xz(2,:)' * c_ribs(ii) + interp1(xyz(:,2),xyz(:,3),y_ribs(ii)) ;
            
            RibsBox_ZCoord_Upper{ii}  = (interp1(x_twisted(1:N),z_twisted(1:N),WB_RibsXCoord0{ii},'PCHIP'));
            RibsBox_ZCoord_Lower{ii}  = (interp1(x_twisted(N+1:end),z_twisted(N+1:end),WB_RibsXCoord0{ii},'PCHIP'));
            
            if 0 && GRAPHFLAG==1
                figure(1); hold on;
                %             plot3(x_twisted,y_ribs(ii)*ones(size(z_twisted)),z_twisted,'cyan')
                plot3(WB_RibsXCoord0{ii},y_ribs(ii)*ones(size(WB_RibsXCoord0{ii})),RibsBox_ZCoord_Upper{ii},'bo','MarkerSize',2)
                plot3(WB_RibsXCoord0{ii},y_ribs(ii)*ones(size(WB_RibsXCoord0{ii})),RibsBox_ZCoord_Lower{ii},'bo','MarkerSize',2)
            end
            
            if 1
                
                x_offset = -min(WB_RibsXCoord0{ii});
                z_offset = -min(RibsBox_ZCoord_Lower{ii});
                x        = WB_RibsXCoord0{ii} + x_offset;
                z_up     = RibsBox_ZCoord_Upper{ii} + z_offset;
                z_low    = RibsBox_ZCoord_Lower{ii} + z_offset;
                
                xz_up    = [x' z_up'];
                xz_low   = [x' z_low'];
                
                RibArea        = trapz(xz_up(:,1),xz_up(:,2)) -  trapz(xz_low(:,1),xz_low(:,2));
                rib_Xcentroid(ii,1) = 1/RibArea * trapz(xz_up(:,1),xz_up(:,1).*(xz_up(:,2)-xz_low(:,2))) - x_offset;
                rib_Ycentroid(ii,1) = y_ribs(ii);
                rib_Zcentroid(ii,1) = 1/RibArea * trapz(xz_up(:,1),0.5*(xz_up(:,2).^2 - xz_low(:,2).^2)) - z_offset;
                
                RibVolume         = RibArea   * RibThickness(ii);
                RibMassDist(ii,1) = RibVolume * RibDensity(ii);
                if 0
                    figure(); hold on;
                    plot(x_twisted,z_twisted,'cyan')
                    plot(WB_RibsXCoord0{ii},RibsBox_ZCoord_Upper{ii},'bo')
                    plot(WB_RibsXCoord0{ii},RibsBox_ZCoord_Lower{ii},'bo')
                    %                 plot(xz_up(:,1),xz_up(:,2),'go')
                    %                 plot(xz_low(:,1),xz_low(:,2),'ro')
                    plot(rib_Xcentroid(ii),rib_Zcentroid(ii),'X','MarkerSize',6)
                end
                
                % Twist rib centroid with wing twist
                R         = [cos(theta_ribs(ii)) sin(theta_ribs(ii)); -sin(theta_ribs(ii)) cos(theta_ribs(ii))];
                
                twisted_xz = R*[rib_Xcentroid(ii,1)- interp1(xyz(:,2),xyz(:,1),y_ribs(ii)); rib_Zcentroid(ii,1)-interp1(xyz(:,2),xyz(:,3),y_ribs(ii))]+[interp1(xyz(:,2),xyz(:,1),y_ribs(ii)); interp1(xyz(:,2),xyz(:,3),y_ribs(ii))];
                rib_Xcentroid(ii) = twisted_xz(1);
                rib_Zcentroid(ii) = twisted_xz(2);
            end
            
            % Store Mass Input (Lumped Formulation)
            lumped.type{1}         = 'Ribs';
            lumped.mass{1}(ii,1)   = RibMassDist(ii,1);                             % [kg]
            lumped.location{1}(ii,:) = [rib_Xcentroid(ii,1), rib_Ycentroid(ii,1), rib_Zcentroid(ii,1)];
            
            IG = [0, 0, 0 ;
                0, 0, 0 ;
                0, 0, 0];
            lumped.IG{1}(:,3*(ii-1)+(1:3))=IG;
        end
    end
    % ---
    
    % ---
    if 1 && GRAPHFLAG == 1    % plot lumped mass
        figure(1); hold on
        plot3(lumped.location{1}(:,1),lumped.location{1}(:,2),lumped.location{1}(:,3),'redd','MarkerFaceColor','none','MarkerSize',6)
    end
    % ---
else
    lumped=struct;
end
