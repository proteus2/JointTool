function [constant,varargout] = LoadAerofoilProfiles(constant,xyz,varargin)

if length(varargin) == 2
    fixnodes = varargin{1};
    morph = varargin{2};
    camberflag = 1;
else
    camberflag = 0;
end

global GRAPHFLAG xlsFileName

% =====                                                              ==== %
%                       LoadAerofoilProfiles [17/07/2015]                 %
%                                                                         %
%  This file loads aerofoils profiles (from an excel file)
%  Interpolate the corresponding profiles at nodes and section locations
%  The profiles and their camber are stored in 
%       constant.inp.Aerofoil.NodeProfiles     
%       constant.inp.Aerofoil.SecProfiles 
%       constant.inp.Aerofoil.SecCamber
%
% Note
%       Keep checking interpolation with new data is added
% =====                                                              ==== %

% Define the y-location of the fixed nodes, in case camber morphing is
% present

if camberflag == 1
    yfix = xyz(fixnodes==1,2);
end

% ---
if 1    % Load Profiles for Excel file called Aerofoil
    % aerofoil have been extracted for IGES file, small errors are likely)

    [~,sheets] = xlsfinfo(xlsFileName);
    index = find(strcmp(sheets,'ProfileYLocation'));
    
    ProfilesLocation = xlsread(xlsFileName,index); 
    
    if camberflag == 1
        [~,~,ProfilesFixInd] = intersect(yfix,ProfilesLocation);
%         ProfilesFix = zeros(numel(ProfilesLocation),1);
%         ProfilesFix(ProfilesFixInd) = 1;
        CamberAirfoil = zeros(numel(ProfilesLocation),1);
%         CamberIni = interp1(yfix,morph.inp.camber.ini,ProfilesLocation);
        CamberIni = zeros(numel(ProfilesLocation),1);
        
        for i=1:length(ProfilesFixInd)
           if morph.inp.camber.loc(i)==1
               if i == 1
                   CamberAirfoil(ProfilesFixInd(i):ProfilesFixInd(i+1)-1) = 1;
                   if morph.inp.camber.loc(i+1) == 1
                       CamberIni(ProfilesFixInd(i):ProfilesFixInd(i+1)-1) = interp1(yfix(i:i+1),morph.inp.camber.ini(i:i+1),ProfilesLocation((ProfilesFixInd(i):ProfilesFixInd(i+1)-1)));
                   else
                       CamberIni(ProfilesFixInd(i):ProfilesFixInd(i+1)-1) = morph.inp.camber.ini(i);
                   end
               elseif i == length(ProfilesFixInd)
                   CamberAirfoil(ProfilesFixInd(i-1)+1:ProfilesFixInd(i)) = 1;
                   if morph.inp.camber.loc(i-1) == 1
                       CamberIni(ProfilesFixInd(i-1)+1:ProfilesFixInd(i)) = interp1(yfix(i-1:i),morph.inp.camber.ini(i-1:i),ProfilesLocation((ProfilesFixInd(i-1)+1:ProfilesFixInd(i))));
                   else
                       CamberIni(ProfilesFixInd(i-1)+1:ProfilesFixInd(i)) = morph.inp.camber.ini(i);
                   end
               else
                   CamberAirfoil(ProfilesFixInd(i):ProfilesFixInd(i+1)-1) = 1;
                   if morph.inp.camber.loc(i+1) == 1
                       CamberIni(ProfilesFixInd(i):ProfilesFixInd(i+1)-1) = interp1(yfix(i:i+1),morph.inp.camber.ini(i:i+1),ProfilesLocation((ProfilesFixInd(i):ProfilesFixInd(i+1)-1)));
                   else
                       CamberIni(ProfilesFixInd(i):ProfilesFixInd(i+1)-1) = morph.inp.camber.ini(i);
                   end
                   CamberAirfoil(ProfilesFixInd(i-1)+1:ProfilesFixInd(i)) = 1;
                   if morph.inp.camber.loc(i-1) == 1
                       CamberIni(ProfilesFixInd(i-1)+1:ProfilesFixInd(i)) = interp1(yfix(i-1:i),morph.inp.camber.ini(i-1:i),ProfilesLocation((ProfilesFixInd(i-1)+1:ProfilesFixInd(i))));
                   else
                       CamberIni(ProfilesFixInd(i-1)+1:ProfilesFixInd(i)) = morph.inp.camber.ini(i);
                   end
               end
           end
        end
    end
    
    NumAerofoil      = numel(ProfilesLocation);
    
    for iSheet = 1 : NumAerofoil
        if camberflag == 1 && CamberAirfoil(iSheet) == 1
            Aerofoil = xlsread(xlsFileName,index+iSheet);
            camber_level{iSheet} = Aerofoil(1,3:4:end)';
            clear XZ_UpperCoord XZ_LowerCoord
            for i=1:size(Aerofoil,2)/4
                XZ_UpperCoord(:,2*(i-1)+(1:2)) = Aerofoil(3:end,4*(i-1)+(1:2));
                XZ_LowerCoord(:,2*(i-1)+(1:2)) = Aerofoil(3:end,4*(i-1)+(3:4));
            end
%             XZ_UpperCoord(isnan(XZ_UpperCoord(:,1)),:)=[];                      % remove lines in case diff. number of points are used for lower and upper coord.
%             XZ_LowerCoord(isnan(XZ_LowerCoord(:,1)),:)=[];                      % remove lines in case diff. number of points are used for lower and upper coord.
            
            Profile_Upper{iSheet} = XZ_UpperCoord;
            Profile_Lower{iSheet} = XZ_LowerCoord;
        else
            Aerofoil = xlsread(xlsFileName,index+iSheet);
            
            XZ_UpperCoord = Aerofoil(3:end,1:2);
            XZ_LowerCoord = Aerofoil(3:end,3:4);
            
%             XZ_UpperCoord(isnan(XZ_UpperCoord(:,1)),:)=[];                      % remove lines in case diff. number of points are used for lower and upper coord.
%             XZ_LowerCoord(isnan(XZ_LowerCoord(:,1)),:)=[];                      % remove lines in case diff. number of points are used for lower and upper coord.
            
            Profile_Upper{iSheet} = XZ_UpperCoord;
            Profile_Lower{iSheet} = XZ_LowerCoord;
        end
    end
end
% ---

% ---
if 1    % Refined profiles with interpolation - All airfoil input will have the same number of points
    
    N1 = 100;                                                               % # of points used to model the leading edge
    N2 = 50;                                                                % # of points used to model the mid profile
    N3 = 75;                                                                % # of points used to model the trailing edge
    N  = N1 + N2 + N3;                                                      % Number of airfoil data points for the upper and lower surface

    x_new = unique([linspace(0,0.15,N1),linspace(0.15,0.65,N2+1),linspace(0.65,1,N3+1)]');      % new x-axis discretization
    
    for i = 1 : NumAerofoil
        if camberflag == 1 && CamberAirfoil(i) == 1
            for k=1:size(Profile_Upper{i},2)/2
                airf_coord_upper = [flipud(x_new), interp1(Profile_Upper{i}(~isnan(Profile_Upper{i}(:,2*(k-1)+1)),2*(k-1)+1) ,Profile_Upper{i}(~isnan(Profile_Upper{i}(:,2*(k-1)+1)),2*(k-1)+2),flipud(x_new) ,'PCHIP')];     % upper surface xy coord.
                airf_coord_lower = [x_new, interp1(Profile_Lower{i}(~isnan(Profile_Lower{i}(:,2*(k-1)+1)),2*(k-1)+1) ,Profile_Lower{i}(~isnan(Profile_Lower{i}(:,2*(k-1)+1)),2*(k-1)+2) ,x_new,'PCHIP')];                     % lower surface xy coord.
                airf_coord{i}(:,2*(k-1)+(1:2))    = [airf_coord_upper; airf_coord_lower];
            end
        else
            airf_coord_upper = [flipud(x_new), interp1(Profile_Upper{i}(~isnan(Profile_Upper{i}(:,1)),1) ,Profile_Upper{i}(~isnan(Profile_Upper{i}(:,1)),2),flipud(x_new) ,'PCHIP')];     % upper surface xy coord.
            airf_coord_lower = [x_new, interp1(Profile_Lower{i}(~isnan(Profile_Lower{i}(:,1)),1) ,Profile_Lower{i}(~isnan(Profile_Lower{i}(:,1)),2) ,x_new,'PCHIP')];                     % lower surface xy coord.
            airf_coord{i}    = [airf_coord_upper; airf_coord_lower];
        end
    end
end
% ---

% ---
if 0 && GRAPHFLAG == 1   % plot Aerofoil profile
    for iaerofoil = 1 : length(airf_coord)
        figure (2)
        hold all
        xlabel('y')
        ylabel('z')
        plot3(airf_coord{iaerofoil}(:,1),ProfilesLocation(iaerofoil)*ones(size(airf_coord{iaerofoil},1),1),airf_coord{iaerofoil}(:,2))
    end
end
% ---
% --- Define initial airfoils
% for iairf = 1:NumAerofoil
%     if camberflag == 1 && CamberAirfoil(iairf) == 1
%         airf_coord_ini{iairf}(:,1)    = interp1(camber_level{iairf}',airf_coord{iairf}(:,1:2:end)',CamberIni(iairf))';
%         airf_coord_ini{iairf}(:,2)    = interp1(camber_level{iairf}',airf_coord{iairf}(:,2:2:end)',CamberIni(iairf))';
%     else
%         airf_coord_ini{iairf}    = airf_coord{iairf};
%     end
% end

% ---
% ---
if 1 % linear interpolation of profiles at node locations and sections
    ProfileNode = cell(size(xyz,1),1);     % at node
    ProfileNodeUpper = cell(size(xyz,1),1);     % at node
    ProfileNodeLower = cell(size(xyz,1),1);     % at node
    if camberflag == 1
        morph.camber.loc = zeros(size(xyz,1),1);
        camberlocs = find(morph.inp.camber.loc==1);
%         morph.camber.param = zeros(size(xyz,1),1);
        morph.camber.dparamddv = zeros(size(xyz,1),sum(morph.inp.camber.loc));
        morph.camber.ini  = zeros(size(xyz,1),1);
%         morph.camber.low = zeros(size(xyz,1),1);
%         morph.camber.high = zeros(size(xyz,1),1);
        morph.camber.data = cell(1,size(xyz,1));
        morph.camber.NodesUpper = cell(1,size(xyz,1));
        morph.camber.NodesLower = cell(1,size(xyz,1));
        morph.camber.level = cell(1,size(xyz,1));
        morph.camber.axis = interp1(yfix,morph.inp.camber.axis,xyz(:,2));
    end
    for iNode = 1 : size(xyz,1)
        IndexProfile1 = find(xyz(iNode,2)>=ProfilesLocation,1,'last');
        IndexProfile2 = find(xyz(iNode,2)<=ProfilesLocation,1);  
        
        if camberflag == 1
            IndexProfileFix1 = find(xyz(iNode,2)>=ProfilesLocation(ProfilesFixInd),1,'last');
            IndexProfileFix2 = find(xyz(iNode,2)<=ProfilesLocation(ProfilesFixInd),1);
        end
        
        if isempty(IndexProfile1)                                                % towards root
            if camberflag == 1
                if CamberAirfoil(1) == 1
                    for i=1:size(airf_coord{1},2)/2
                        airfdat(:,2*(i-1)+(1:2)) = [airf_coord{1}(:,2*(i-1)+1), airf_coord{1}(:,2*(i-1)+2)];
                    end
                    morph.camber.loc(iNode) = 1;
                    morph.camber.level{iNode} = camber_level{1};
%                     morph.camber.param(iNode) = morph.inp.camber.param(1);
                    morph.camber.dparamddv(iNode,1) = 1;
                    
                    morph.camber.ini(iNode) = CamberIni(1);
%                     morph.camber.low(iNode) = morph.inp.camber.low(1);
%                     morph.camber.high(iNode) = morph.inp.camber.high(1);
                    ProfileNode{iNode}(:,1) = interp1(morph.camber.level{iNode},airfdat(:,1:2:end)',morph.camber.ini(iNode))';
                    ProfileNode{iNode}(:,2) = interp1(morph.camber.level{iNode},airfdat(:,2:2:end)',morph.camber.ini(iNode))';
                    airfin = [ProfileNode{iNode}(N:-1:1,1), 0.5*(ProfileNode{iNode}(N:-1:1,2) + ProfileNode{iNode}(N+1:end,2))];
                    shiftin = interp1(airfin(:,1),airfin(:,2),morph.camber.axis(iNode));
                    for i=1:size(airf_coord{1},2)/2
                        airf = [airf_coord{1}(:,2*(i-1)+1), airf_coord{1}(:,2*(i-1)+2)];
                        camberdat = [airf(N:-1:1,1), 0.5*(airf(N:-1:1,2) + airf(N+1:end,2))];
                        shift = interp1(camberdat(:,1),camberdat(:,2),morph.camber.axis(iNode))-shiftin;
                        morph.camber.NodesUpper{iNode}(:,i) = reshape([airf_coord{1}(N:-1:1,2*(i-1)+1), airf_coord{1}(N:-1:1,2*(i-1)+2)-shift]',[],1);
                        morph.camber.NodesLower{iNode}(:,i) = reshape([airf_coord{1}(N:-1:1,2*(i-1)+1), airf_coord{1}(N+1:end,2*(i-1)+2)-shift]',[],1);
                        morph.camber.data{iNode}(:,i) = reshape([airf_coord{1}(N:-1:1,2*(i-1)+1), 0.5*(airf_coord{1}(N:-1:1,2*(i-1)+2) + airf_coord{1}(N+1:end,2*(i-1)+2))-shift]',[],1);
                    end
                else
                    ProfileNode{iNode} = airf_coord{1};
                end
            else
                ProfileNode{iNode} = airf_coord{1};
            end
        elseif isempty(IndexProfile2)                                            % towards tip
            if camberflag == 1
                if CamberAirfoil(end) == 1
                    for i=1:size(airf_coord{1},2)/2
                        airfdat(:,2*(i-1)+(1:2)) = [airf_coord{end}(:,2*(i-1)+1), airf_coord{end}(:,2*(i-1)+2)];
                    end
                    morph.camber.loc(iNode) = 1;
                    morph.camber.level{iNode} = camber_level{end};
%                     morph.camber.param(iNode) = morph.inp.camber.param(end);
                    morph.camber.dparamddv(iNode,end) = 1;
                    
                    morph.camber.ini(iNode) = CamberIni(end);
%                     morph.camber.low(iNode) = morph.inp.camber.low(end);
%                     morph.camber.high(iNode) = morph.inp.camber.high(end);
                    ProfileNode{iNode}(:,1) = interp1(morph.camber.level{iNode},airfdat(:,1:2:end)',morph.camber.ini(iNode))';
                    ProfileNode{iNode}(:,2) = interp1(morph.camber.level{iNode},airfdat(:,2:2:end)',morph.camber.ini(iNode))';
                    airfin = [ProfileNode{iNode}(N:-1:1,1), 0.5*(ProfileNode{iNode}(N:-1:1,2) + ProfileNode{iNode}(N+1:end,2))];
                    shiftin = interp1(airfin(:,1),airfin(:,2),morph.camber.axis(iNode));
                    for i=1:size(airf_coord{end},2)/2
                        airf = [airf_coord{end}(:,2*(i-1)+1), airf_coord{end}(:,2*(i-1)+2)];
                        camberdat = [airf(N:-1:1,1), 0.5*(airf(N:-1:1,2) + airf(N+1:end,2))];
                        shift = interp1(camberdat(:,1),camberdat(:,2),morph.camber.axis(iNode))-shiftin;
                        morph.camber.NodesUpper{iNode}(:,i) = reshape([airf_coord{end}(N:-1:1,2*(i-1)+1), airf_coord{end}(N:-1:1,2*(i-1)+2)-shift]',[],1);
                        morph.camber.NodesLower{iNode}(:,i) = reshape([airf_coord{end}(N:-1:1,2*(i-1)+1), airf_coord{end}(N+1:end,2*(i-1)+2)-shift]',[],1);
                        morph.camber.data{iNode}(:,i) = reshape([airf_coord{end}(N:-1:1,2*(i-1)+1), 0.5*(airf_coord{end}(N:-1:1,2*(i-1)+2) + airf_coord{end}(N+1:end,2*(i-1)+2))-shift]',[],1);
                    end
                else
                    ProfileNode{iNode} = airf_coord{end};
                end
            else
                ProfileNode{iNode} = airf_coord{end};
            end
        elseif IndexProfile1 == IndexProfile2
            if camberflag == 1
                if CamberAirfoil(IndexProfile1) == 1
                    for i=1:size(airf_coord{IndexProfile1},2)/2
                        airfdat(:,2*(i-1)+(1:2)) = [airf_coord{IndexProfile1}(:,2*(i-1)+1), airf_coord{IndexProfile1}(:,2*(i-1)+2)];
                    end
                    morph.camber.loc(iNode) = 1;
                    morph.camber.level{iNode} = camber_level{IndexProfile1};
%                     morph.camber.param(iNode) = morph.inp.camber.param(IndexProfile1);
                    morph.camber.dparamddv(iNode,camberlocs==IndexProfileFix1) = 1;
                    
                    morph.camber.ini(iNode) = CamberIni(IndexProfile1);
%                     morph.camber.low(iNode) = morph.inp.camber.low(IndexProfile1);
%                     morph.camber.high(iNode) = morph.inp.camber.high(IndexProfile1);
                    ProfileNode{iNode}(:,1) = interp1(morph.camber.level{iNode},airfdat(:,1:2:end)',morph.camber.ini(iNode))';
                    ProfileNode{iNode}(:,2) = interp1(morph.camber.level{iNode},airfdat(:,2:2:end)',morph.camber.ini(iNode))';
                    airfin = [ProfileNode{iNode}(N:-1:1,1), 0.5*(ProfileNode{iNode}(N:-1:1,2) + ProfileNode{iNode}(N+1:end,2))];
                    shiftin = interp1(airfin(:,1),airfin(:,2),morph.camber.axis(iNode));
                    for i=1:size(airf_coord{IndexProfile1},2)/2
                        camberdat = [airf_coord{IndexProfile1}(N:-1:1,2*(i-1)+1), 0.5*(airf_coord{IndexProfile1}(N:-1:1,2*(i-1)+2) + airf_coord{IndexProfile1}(N+1:end,2*(i-1)+2))];
                        shift = interp1(camberdat(:,1),camberdat(:,2),morph.camber.axis(iNode))-shiftin;
                        morph.camber.NodesUpper{iNode}(:,i) = reshape([airf_coord{IndexProfile1}(N:-1:1,2*(i-1)+1), airf_coord{IndexProfile1}(N:-1:1,2*(i-1)+2)-shift]',[],1);
                        morph.camber.NodesLower{iNode}(:,i) = reshape([airf_coord{IndexProfile1}(N:-1:1,2*(i-1)+1), airf_coord{IndexProfile1}(N+1:end,2*(i-1)+2)-shift]',[],1);
                        morph.camber.data{iNode}(:,i) = reshape([airf_coord{IndexProfile1}(N:-1:1,2*(i-1)+1), 0.5*(airf_coord{IndexProfile1}(N:-1:1,2*(i-1)+2) + airf_coord{IndexProfile1}(N+1:end,2*(i-1)+2))-shift]',[],1);
                    end
                else
                    ProfileNode{iNode} = airf_coord{IndexProfile1};
                end
            else
                ProfileNode{iNode} = airf_coord{IndexProfile1};
            end
        else                                                                % in between two profiles
            if camberflag == 1
                if CamberAirfoil(IndexProfile1) == 1 && CamberAirfoil(IndexProfile2) == 1
                    for i=1:size(airf_coord{IndexProfile1},2)/2
                        airfdat(:,2*(i-1)+(1:2)) = (airf_coord{IndexProfile2}(:,2*(i-1)+(1:2)) - airf_coord{IndexProfile1}(:,2*(i-1)+(1:2)))...
                            /(ProfilesLocation(IndexProfile2)-ProfilesLocation(IndexProfile1))  * (xyz(iNode,2)-ProfilesLocation(IndexProfile1)) + airf_coord{IndexProfile1}(:,2*(i-1)+(1:2));
                    end
                    morph.camber.level{iNode} = camber_level{IndexProfile1};
                    morph.camber.loc(iNode) = 1;
%                     morph.camber.param(iNode) = interp1([ProfilesLocation(IndexProfile1),ProfilesLocation(IndexProfile2)],[morph.inp.camber.param(IndexProfile1),morph.inp.camber.param(IndexProfile2)],xyz(iNode,2));
                    if sum(camberlocs==IndexProfileFix1)==0
                        morph.camber.dparamddv(iNode,camberlocs==IndexProfileFix2) = 1;
                    elseif sum(camberlocs==IndexProfileFix2)==0
                        morph.camber.dparamddv(iNode,camberlocs==IndexProfileFix1) = 1;
                    else
                        morph.camber.dparamddv(iNode,camberlocs==IndexProfileFix1) = (ProfilesLocation(ProfilesFixInd(IndexProfileFix2))-xyz(iNode,2))/(ProfilesLocation(ProfilesFixInd(IndexProfileFix2))-ProfilesLocation(ProfilesFixInd(IndexProfileFix1)));
                        morph.camber.dparamddv(iNode,camberlocs==IndexProfileFix2) = (xyz(iNode,2)-ProfilesLocation(ProfilesFixInd(IndexProfileFix1)))/(ProfilesLocation(ProfilesFixInd(IndexProfileFix2))-ProfilesLocation(ProfilesFixInd(IndexProfileFix1)));
                    end
                    
                    morph.camber.ini(iNode) = interp1([ProfilesLocation(IndexProfile1),ProfilesLocation(IndexProfile2)],[CamberIni(IndexProfile1),CamberIni(IndexProfile2)],xyz(iNode,2));
%                     morph.camber.low(iNode) = interp1([ProfilesLocation(IndexProfile1),ProfilesLocation(IndexProfile2)],[morph.inp.camber.low(IndexProfile1),morph.inp.camber.low(IndexProfile2)],xyz(iNode,2));
%                     morph.camber.high(iNode) = interp1([ProfilesLocation(IndexProfile1),ProfilesLocation(IndexProfile2)],[morph.inp.camber.high(IndexProfile1),morph.inp.camber.high(IndexProfile2)],xyz(iNode,2));
                    ProfileNode{iNode}(:,1) = interp1(morph.camber.level{iNode},airfdat(:,1:2:end)',morph.camber.ini(iNode))';
                    ProfileNode{iNode}(:,2) = interp1(morph.camber.level{iNode},airfdat(:,2:2:end)',morph.camber.ini(iNode))';
                    airfin = [ProfileNode{iNode}(N:-1:1,1), 0.5*(ProfileNode{iNode}(N:-1:1,2) + ProfileNode{iNode}(N+1:end,2))];
                    shiftin = interp1(airfin(:,1),airfin(:,2),morph.camber.axis(iNode));
                    for i=1:size(airf_coord{IndexProfile1},2)/2
                        airf = (airf_coord{IndexProfile2}(:,2*(i-1)+(1:2)) - airf_coord{IndexProfile1}(:,2*(i-1)+(1:2)))...
                            /(ProfilesLocation(IndexProfile2)-ProfilesLocation(IndexProfile1))  * (xyz(iNode,2)-ProfilesLocation(IndexProfile1)) + airf_coord{IndexProfile1}(:,2*(i-1)+(1:2));
                        camberdat = [airf(N:-1:1,1), 0.5*(airf(N:-1:1,2) + airf(N+1:end,2))];
                        shift = interp1(camberdat(:,1),camberdat(:,2),morph.camber.axis(iNode))-shiftin;
                        morph.camber.NodesUpper{iNode}(:,i) = reshape([airf(N:-1:1,1), airf(N:-1:1,2)-shift]',[],1);
                        morph.camber.NodesLower{iNode}(:,i) = reshape([airf(N:-1:1,1), airf(N+1:end,2)-shift]',[],1);
                        morph.camber.data{iNode}(:,i) = reshape([airf(N:-1:1,1), 0.5*(airf(N:-1:1,2) + airf(N+1:end,2))-shift]',[],1);
                    end
                elseif CamberAirfoil(IndexProfile1) == 1
                    for i=1:size(airf_coord{IndexProfile1},2)/2
                        airfdat(:,2*(i-1)+(1:2)) = (airf_coord{IndexProfile2} - airf_coord{IndexProfile1}(:,2*(i-1)+(1:2)))...
                            /(ProfilesLocation(IndexProfile2)-ProfilesLocation(IndexProfile1))  * (xyz(iNode,2)-ProfilesLocation(IndexProfile1)) + airf_coord{IndexProfile1}(:,2*(i-1)+(1:2));
                    end
                    morph.camber.level{iNode} = camber_level{IndexProfile1};
                    morph.camber.loc(iNode) = 1;
%                     morph.camber.param(iNode) = morph.inp.camber.param(IndexProfile1);
                    morph.camber.dparamddv(iNode,camberlocs==IndexProfileFix1) = 1;
                    
                    morph.camber.ini(iNode) = CamberIni(IndexProfile1);
%                     morph.camber.low(iNode) = morph.inp.camber.low(IndexProfile1);
%                     morph.camber.high(iNode) = morph.inp.camber.high(IndexProfile1);
                    ProfileNode{iNode}(:,1) = interp1(morph.camber.level{iNode},airfdat(:,1:2:end)',morph.camber.ini(iNode))';
                    ProfileNode{iNode}(:,2) = interp1(morph.camber.level{iNode},airfdat(:,2:2:end)',morph.camber.ini(iNode))';
                    airfin = [ProfileNode{iNode}(N:-1:1,1), 0.5*(ProfileNode{iNode}(N:-1:1,2) + ProfileNode{iNode}(N+1:end,2))];
                    shiftin = interp1(airfin(:,1),airfin(:,2),morph.camber.axis(iNode));
                    for i=1:size(airf_coord{IndexProfile1},2)/2
                        airf = (airf_coord{IndexProfile2} - airf_coord{IndexProfile1}(:,2*(i-1)+(1:2)))...
                            /(ProfilesLocation(IndexProfile2)-ProfilesLocation(IndexProfile1))  * (xyz(iNode,2)-ProfilesLocation(IndexProfile1)) + airf_coord{IndexProfile1}(:,2*(i-1)+(1:2));
                        camberdat = [airf(N:-1:1,1), 0.5*(airf(N:-1:1,2) + airf(N+1:end,2))];
                        shift = interp1(camberdat(:,1),camberdat(:,2),morph.camber.axis(iNode))-shiftin;
                        morph.camber.NodesUpper{iNode}(:,i) = reshape([airf(N:-1:1,1), airf(N:-1:1,2)-shift]',[],1);
                        morph.camber.NodesLower{iNode}(:,i) = reshape([airf(N:-1:1,1), airf(N+1:end,2)-shift]',[],1);
                        morph.camber.data{iNode}(:,i) = reshape([airf(N:-1:1,1), 0.5*(airf(N:-1:1,2) + airf(N+1:end,2))-shift]',[],1);
                    end
                elseif CamberAirfoil(IndexProfile2) == 1
                    for i=1:size(airf_coord{IndexProfile2},2)/2
                        airfdat(:,2*(i-1)+(1:2)) = (airf_coord{IndexProfile2}(:,2*(i-1)+(1:2)) - airf_coord{IndexProfile1})...
                            /(ProfilesLocation(IndexProfile2)-ProfilesLocation(IndexProfile1))  * (xyz(iNode,2)-ProfilesLocation(IndexProfile1)) + airf_coord{IndexProfile1};
                    end
                    morph.camber.level{iNode} = camber_level{IndexProfile2};
                    morph.camber.loc(iNode) = 1;
%                     morph.camber.param(iNode) = morph.inp.camber.param(IndexProfile2);
                    morph.camber.dparamddv(iNode,camberlocs==IndexProfileFix2) = 1;
                    
                    morph.camber.ini(iNode) = CamberIni(IndexProfile2);
%                     morph.camber.low(iNode) = morph.inp.camber.low(IndexProfile2);
%                     morph.camber.high(iNode) = morph.inp.camber.high(IndexProfile2);
                    ProfileNode{iNode}(:,1) = interp1(morph.camber.level{iNode},airfdat(:,1:2:end)',morph.camber.ini(iNode))';
                    ProfileNode{iNode}(:,2) = interp1(morph.camber.level{iNode},airfdat(:,2:2:end)',morph.camber.ini(iNode))';
                    airfin = [ProfileNode{iNode}(N:-1:1,1), 0.5*(ProfileNode{iNode}(N:-1:1,2) + ProfileNode{iNode}(N+1:end,2))];
                    shiftin = interp1(airfin(:,1),airfin(:,2),morph.camber.axis(iNode));
                    for i=1:size(airf_coord{IndexProfile2},2)/2
                        airf = (airf_coord{IndexProfile2}(:,2*(i-1)+(1:2)) - airf_coord{IndexProfile1})...
                            /(ProfilesLocation(IndexProfile2)-ProfilesLocation(IndexProfile1))  * (xyz(iNode,2)-ProfilesLocation(IndexProfile1)) + airf_coord{IndexProfile1};
                        camberdat = [airf(N:-1:1,1), 0.5*(airf(N:-1:1,2) + airf(N+1:end,2))];
                        shift = interp1(camberdat(:,1),camberdat(:,2),morph.camber.axis(iNode))-shiftin;
                        morph.camber.NodesUpper{iNode}(:,i) = reshape([airf(N:-1:1,1), airf(N:-1:1,2)-shift]',[],1);
                        morph.camber.NodesLower{iNode}(:,i) = reshape([airf(N:-1:1,1), airf(N+1:end,2)-shift]',[],1);
                        morph.camber.data{iNode}(:,i) = reshape([airf(N:-1:1,1), 0.5*(airf(N:-1:1,2) + airf(N+1:end,2))-shift]',[],1);
                    end
                else
                    ProfileNode{iNode} = (airf_coord{IndexProfile2} - airf_coord{IndexProfile1})...
                                /(ProfilesLocation(IndexProfile2)-ProfilesLocation(IndexProfile1))  * (xyz(iNode,2)-ProfilesLocation(IndexProfile1)) + airf_coord{IndexProfile1};
                end
            else
                ProfileNode{iNode} = (airf_coord{IndexProfile2} - airf_coord{IndexProfile1})...
                                /(ProfilesLocation(IndexProfile2)-ProfilesLocation(IndexProfile1))  * (xyz(iNode,2)-ProfilesLocation(IndexProfile1)) + airf_coord{IndexProfile1};
            end
        end
        ProfileNodeUpper{iNode} = ProfileNode{iNode}(1:N,:);
        ProfileNodeLower{iNode} = ProfileNode{iNode}(N+1:end,:);
    end
    
    if 0 && GRAPHFLAG == 1 % plot Aerofoil profile interpolated (visual check)
        figure (2)
        for iaerofoil = 1: size(xyz,1)
            xlabel('y')
            ylabel('z')
            hold off
            plot3(ProfileNode{iaerofoil}(:,1),xyz(iaerofoil,2)*ones(size(ProfileNode{iaerofoil},1),1),ProfileNode{iaerofoil}(:,2),'--')
            if camberflag == 1
                hold on
                if morph.camber.loc(iaerofoil) == 1
                    for i=1:size(morph.camber.data{iaerofoil},2)
                        plot3(morph.camber.data{iaerofoil}(1:2:end,i),xyz(iaerofoil,2)*ones(size(morph.camber.data{iaerofoil}(1:2:end,i),1),1),morph.camber.data{iaerofoil}(2:2:end,i),'-.')
                    end
                end
            end
            pause
        end
    end
end

% ---
% --- Calculate aerofoils camber lines
camber = nan*ones(N*2,size(xyz,1));
thickness = zeros(size(xyz,1),1);
for ii=1:size(xyz,1)
    camber(:,ii) = reshape([ProfileNode{ii}(N:-1:1,1), 0.5*(ProfileNode{ii}(N:-1:1,2) + ProfileNode{ii}(N+1:end,2))]',[],1); % Camberline
    thickness(ii,1) = max(ProfileNodeUpper{ii}(end:-1:1,2)-ProfileNodeLower{ii}(:,2));
    
    if 0 && GRAPHFLAG == 1    % (visual check)
        figure(3 + ii)
        hold on
        plot(ProfileNode{ii}(:,1),ProfileNode{ii}(:,2),'blue')
        plot(camber(1:2:end,ii),camber(2:2:end,ii),'red')
    end
end
% ---

% --- Calculate scaling factor for chord to keep length of the camber line
% equal
if camberflag == 1 
    morph.camber.scale = cell(1,size(xyz,1));
    for ii=1:size(xyz,1)
        if morph.camber.loc(ii) == 1
            for i = 1:size(morph.camber.data{ii},2)
                options = optimset('Display','off');
                morph.camber.scale{ii}(i,1) = fsolve(@(x) scale_fun(camber(:,ii),morph.camber.data{ii}(:,i),x),0.99,options);
            end
        end
    end
end
% ---

if camberflag == 1
    varargout{1}=morph;
end

constant.inp.Aerofoil.NodeProfilesUpper = ProfileNodeUpper;               % Save aerofoil profiles
constant.inp.Aerofoil.NodeProfilesLower = ProfileNodeLower;               % Save aerofoil profiles
constant.inp.Aerofoil.NodeCamber        = camber;                         % Save aerofoil profiles
constant.inp.Aerofoil.Thickness         = thickness;                      % Save aerofoil profiles