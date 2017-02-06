function [cross,constant,panels] = GenerateWingbox(wing_data,constant,str,cross,xyz,xyz_LE,theta,xref,xsref,c,Nlam,stiff,FLAGSPAR,FLAGSTRING,varargin)

if length(varargin) == 2
    morph = varargin{1};
    fixnodes = varargin{2};
    morphflag = 1;
else
    morphflag = 0;
end

global GRAPHFLAG

nairf = 20; % Number of data points per airfoil surface
nsec = length(str.xyz)/3-1;

wb_le_right = wing_data(:,9);    % in % of chord
wb_le_left  = wb_le_right(2:end);
wb_le       = [flipud(wb_le_left); wb_le_right];

wb_te_right = wing_data(:,10);   % in % of chord
wb_te_left  = wb_te_right(2:end);
wb_te       = [flipud(wb_te_left); wb_te_right];

wb_le_str = interp1(xyz(:,2),wb_le.*c + xyz_LE(:,1),str.xyz(2:3:end),'linear');
wb_te_str = interp1(xyz(:,2),wb_te.*c + xyz_LE(:,1),str.xyz(2:3:end),'linear');

sparloc_right = wing_data(:,12:end);
sparloc_left  = sparloc_right(2:end,:);
sparloc_vec   = [flipud(sparloc_left); sparloc_right];
 
% c_str_0 = (wb_te - wb_le).*c;
% c_str   = interp1(xyz(:,2),c_str_0,str.xyz(2:3:end),'linear');
%  
% % store str chord values
% constant.inp.cbox = c_str;

% In case of span extension
if morphflag == 1 && morph.inp.span.flag == 1
    xyz_fix = xyz(fixnodes==1,:);
    spanflag = 1;
    inbsection = find(str.xyz(2:3:end)>=xyz_fix(morph.inp.span.fixedsec,2) & str.xyz(2:3:end)<xyz_fix(morph.inp.span.fixedsec+1,2));
    doublesection = find(str.xyz(2:3:end)>=xyz_fix(morph.inp.span.doublesec,2) & str.xyz(2:3:end)<xyz_fix(morph.inp.span.doublesec+1,2));
    outbsection = find(str.xyz(2:3:end)>=xyz_fix(morph.inp.span.extsec,2) & str.xyz(2:3:end)<xyz_fix(morph.inp.span.extsec+1,2));
else
    spanflag = 0;
end

countstringer = 0;

% Define the wingbox
for i=1:nsec
    if spanflag == 1 && ~isempty(find(inbsection==i, 1))
        y_sec      = (xyz_fix(morph.inp.span.fixedsec,2)+xyz_fix(morph.inp.span.fixedsec+1,2))/2;
    elseif spanflag == 1 && ~isempty(find(outbsection==i, 1))
        y_sec      = (xyz_fix(morph.inp.span.extsec,2)+xyz_fix(morph.inp.span.extsec+1,2))/2;
    else
        y_sec      = (str.xyz(3*(i-1)+2)+str.xyz(3*i+2))/2;
    end
    if spanflag == 0 || isempty(find(doublesection==i, 1))
        theta_sec  = interp1(xyz(:,2),theta,y_sec);
        xref_sec   = interp1(xyz(:,2),xref,y_sec);
        xsref_sec  = interp1(xyz(:,2),xsref,y_sec);
%         wb_le_sec  = interp1(xyz(:,2),wing_data(:,9),y_sec);
        wb_le_sec  = interp1(xyz(:,2),wb_le,y_sec);
%         wb_te_sec  = interp1(xyz(:,2),wing_data(:,10),y_sec);
        wb_te_sec  = interp1(xyz(:,2),wb_te,y_sec);
        c_sec      = interp1(xyz(:,2),c,y_sec);
        lam_sec    = find(y_sec<=Nlam.xyz_lam(:,2),1)-1;
        if FLAGSTRING
            stiff_pitch_sec = interp1(stiff.yloc,stiff.pitch,y_sec);
            stiff_height_sec = interp1(stiff.yloc,stiff.height,y_sec);
            stiff_EA_sec = interp1(stiff.yloc,stiff.EA,y_sec);
            stiff_mA_sec = interp1(stiff.yloc,stiff.mA,y_sec);
        end
        if FLAGSPAR
%             sparloc = interp1(xyz(:,2),wing_data(:,12:end),y_sec);
            sparloc = interp1(xyz(:,2),sparloc_vec,y_sec);
            nspars = length(sparloc)-sum(isnan(sparloc));
            sparloc(isnan(sparloc))=[];
        end
        
        IndexProfile1 = find(y_sec>=xyz(:,2),1,'last');
        IndexProfile2 = find(y_sec<=xyz(:,2),1);
        
        ProfileWingbox{i} = ([constant.inp.Aerofoil.NodeProfilesUpper{IndexProfile2};constant.inp.Aerofoil.NodeProfilesLower{IndexProfile2}] - [constant.inp.Aerofoil.NodeProfilesUpper{IndexProfile1};constant.inp.Aerofoil.NodeProfilesLower{IndexProfile1}])...
            /(xyz(IndexProfile2,2)-xyz(IndexProfile1,2)+eps)  * (y_sec-xyz(IndexProfile1,2)) + [constant.inp.Aerofoil.NodeProfilesUpper{IndexProfile1};constant.inp.Aerofoil.NodeProfilesLower{IndexProfile1}];
        
        % Rotate the airfoil (apply twist)
        R          = eye(2);%[cos(theta_sec) sin(theta_sec); -sin(theta_sec) cos(theta_sec)];
        x_local    =  ProfileWingbox{i}(:,1);                                      % normalised section profile
        z_local    =  ProfileWingbox{i}(:,2);
        x_local    =  x_local - xref_sec+xsref_sec;                                           % twisted about the quarter chord ( this is an assumption !!! )
        twisted_xz =  R*[x_local' ; z_local'];                                  % Apply rotation
        x_twisted  = (twisted_xz(1,:)' + xref_sec-xsref_sec);                               % remove quarter chord shift
        z_twisted  =  twisted_xz(2,:)';
        
        if FLAGSPAR
            if Nlam.Chord(lam_sec) == 1
                lamloc = [wb_le_sec,wb_te_sec];
                skinloc = unique([wb_le_sec,sparloc,wb_te_sec]);
                if FLAGSTRING
                    nstring = ceil((wb_te_sec-wb_le_sec)/(stiff_pitch_sec/c_sec));
                    stringloc = interpfix(skinloc,nstring);
                    WB_X{i} = interpfix(stringloc,nairf);
                else
                    WB_X{i} = interpfix(skinloc,nairf);
                end
            else
                skinloc = unique([wb_le_sec,sparloc,wb_te_sec]);
                if mod(Nlam.Chord(lam_sec),(length(skinloc)-1))~=0
                    error(['Please specify a number of chordwise laminates equal to 1 or a multiple of the number of laminate bays for section ',num2str(lam_sec)])
                else
                    lamloc = [];
                    for ispar = 1:length(skinloc)-1
                        lamloc = [lamloc(1:end-1),linspace(skinloc(ispar),skinloc(ispar+1),(Nlam.Chord(lam_sec)/(length(skinloc)-1))+1)];
                    end
                    if FLAGSTRING
                        nstring = ceil((wb_te_sec-wb_le_sec)/(stiff_pitch_sec/c_sec));
                        stringloc = interpfix(lamloc,nstring);
                        WB_X{i} = interpfix(stringloc,nairf);
                    else
                        WB_X{i} = interpfix(lamloc,nairf);
                    end
                end
            end
        else
            if Nlam.Chord(lam_sec) == 1
                lamloc = [wb_le_sec,wb_te_sec];
                if FLAGSTRING
                    nstring = ceil((wb_te_sec-wb_le_sec)/(stiff_pitch_sec/c_sec));
                    stringloc = interpfix(lamloc,nstring);
                    WB_X{i} = interpfix(stringloc,nairf);
                else
                    WB_X{i} = interpfix(lamloc,nairf);
                end
            else
                lamloc = linspace(wb_le_sec,wb_te_sec,Nlam.Chord(lam_sec)+1);
                if FLAGSTRING
                    nstring = ceil((wb_te_sec-wb_le_sec)/(stiff_pitch_sec/c_sec));
                    stringloc = interpfix(lamloc,nstring);
                    WB_X{i} = interpfix(stringloc,nairf);
                else
                    WB_X{i} = interpfix(lamloc,nairf);
                end
            end
        end
        
        Nupper = size(constant.inp.Aerofoil.NodeProfilesUpper{IndexProfile1},1);
        WB_Z_Upper{i} = interp1(x_twisted(1:Nupper),z_twisted(1:Nupper),WB_X{i},'PCHIP');             % normalised WB upper Z coord.
        WB_Z_Lower{i} = interp1(x_twisted(Nupper+1:end),z_twisted(Nupper+1:end),flipud(WB_X{i}),'PCHIP');     % normalised WB lower Z coord.
        
        
        cross.yzlocal{i}{1} = unique([WB_X{i},WB_Z_Upper{i};flipud(WB_X{i}),WB_Z_Lower{i}],'rows','stable');
        
        % Define cross-section connectivity
        % Top skin
        countelm = 0;
        cross.elmID{i}{1} = [];
        cross.elmloc{i}{1} = [];
        cross.lam{i}{1} = [];
        for ilam=1:Nlam.Chord(lam_sec)
            node1 = find(cross.yzlocal{i}{1}(:,1)==lamloc(ilam),1);
            node2 = find(cross.yzlocal{i}{1}(:,1)==lamloc(ilam+1),1);
            cross.elmID{i}{1} = [cross.elmID{i}{1};countelm+(1:node2-node1)'];
            countelm = countelm + node2-node1;
            cross.elmloc{i}{1} = [cross.elmloc{i}{1};[(node1:node2-1)',(node1+1:node2)']];
            cross.lam{i}{1} = [cross.lam{i}{1};constant.lam.TopID{lam_sec}(ilam)*ones(node2-node1,1)];
        end
        
        % Define buckling panels for the top skin
        if FLAGSTRING
            [~,panelnodes] = ismember(stringloc,cross.yzlocal{i}{1}(1:end/2,1));
            panels.yzlocal{i}{1} = cross.yzlocal{i}{1}(panelnodes,:);
            panels.elmloc{i}{1} = [(1:length(panelnodes)-1)',(2:length(panelnodes))'];
            % Identify element number and corresponding laminate by looking for elm where first node corresponds to the first panelnode
            [~,elmnum] = ismember(panelnodes(1:end-1),cross.elmloc{i}{1}(:,1));
            panels.lam{i}{1} = cross.lam{i}{1}(elmnum);
            panels.type{i}{1} = ones(size(cross.lam{i}{1}(elmnum))); %1: Top skin
            % Identify all cross-sectional elements that belong to a certain
            % panel
            for j=1:length(panelnodes)-1
                elm1 = find(cross.elmloc{i}{1}(:,1)==panelnodes(j));
                elm2 = find(cross.elmloc{i}{1}(:,2)==panelnodes(j+1));
                panels.elmnum{i}{1}{j} = (elm1:elm2)';
            end
        else
            if FLAGSPAR
                if Nlam.Chord(lam_sec) == 1
                    [~,panelnodes] = ismember(skinloc,cross.yzlocal{i}{1}(1:end/2,1));
                    panels.yzlocal{i}{1} = cross.yzlocal{i}{1}(panelnodes,:);
                    panels.elmloc{i}{1} = [(1:length(panelnodes)-1)',(2:length(panelnodes))'];
                    % Identify element number and corresponding laminate by looking for elm where first node corresponds to the first panelnode
                    [~,elmnum] = ismember(panelnodes(1:end-1),cross.elmloc{i}{1}(:,1));
                    panels.lam{i}{1} = cross.lam{i}{1}(elmnum);
                    panels.type{i}{1} = ones(size(cross.lam{i}{1}(elmnum))); %1: Top skin
                    % Identify all cross-sectional elements that belong to a certain
                    % panel
                    for j=1:length(panelnodes)-1
                        elm1 = find(cross.elmloc{i}{1}(:,1)==panelnodes(j));
                        elm2 = find(cross.elmloc{i}{1}(:,2)==panelnodes(j+1));
                        panels.elmnum{i}{1}{j} = (elm1:elm2)';
                    end
                else
                    [~,panelnodes] = ismember(lamloc,cross.yzlocal{i}{1}(1:end/2,1));
                    panels.yzlocal{i}{1} = cross.yzlocal{i}{1}(panelnodes,:);
                    panels.elmloc{i}{1} = [(1:length(panelnodes)-1)',(2:length(panelnodes))'];
                    % Identify element number and corresponding laminate by looking for elm where first node corresponds to the first panelnode
                    [~,elmnum] = ismember(panelnodes(1:end-1),cross.elmloc{i}{1}(:,1));
                    panels.lam{i}{1} = cross.lam{i}{1}(elmnum);
                    panels.type{i}{1} = ones(size(cross.lam{i}{1}(elmnum))); %1: Top skin
                    % Identify all cross-sectional elements that belong to a certain
                    % panel
                    for j=1:length(panelnodes)-1
                        elm1 = find(cross.elmloc{i}{1}(:,1)==panelnodes(j));
                        elm2 = find(cross.elmloc{i}{1}(:,2)==panelnodes(j+1));
                        panels.elmnum{i}{1}{j} = (elm1:elm2)';
                    end
                end
            end
        end
        
        % Bottom skin
        for ilam=Nlam.Chord(lam_sec):-1:1
            node1 = find(cross.yzlocal{i}{1}(:,1)==lamloc(ilam),1,'last');
            node2 = find(cross.yzlocal{i}{1}(:,1)==lamloc(ilam+1),1,'last');
            if node1 == 1
                nelcross = size(cross.yzlocal{i}{1},1);
                node1 = nelcross;
                cross.elmID{i}{1} = [cross.elmID{i}{1};countelm+(1:node1-node2)'];
                countelm = countelm + node1-node2;
                cross.elmloc{i}{1} = [cross.elmloc{i}{1};[(node2:node1-1)',(node2+1:node1)']];
                cross.lam{i}{1} = [cross.lam{i}{1};constant.lam.BotID{lam_sec}(ilam)*ones(node1-node2,1)];
                
                % Add last element to close the cross-section
                cross.elmID{i}{1} = [cross.elmID{i}{1};countelm+1];
                countelm = countelm + 1;
                cross.elmloc{i}{1} = [cross.elmloc{i}{1};[nelcross,1]];
                cross.lam{i}{1} = [cross.lam{i}{1};constant.lam.BotID{lam_sec}(ilam)];
            else
                cross.elmID{i}{1} = [cross.elmID{i}{1};countelm+(1:node1-node2)'];
                countelm = countelm + node1-node2;
                cross.elmloc{i}{1} = [cross.elmloc{i}{1};[(node2:node1-1)',(node2+1:node1)']];
                cross.lam{i}{1} = [cross.lam{i}{1};constant.lam.BotID{lam_sec}(ilam)*ones(node1-node2,1)];
            end
        end
        
        % Define buckling panels for the bottom skin
        if FLAGSTRING
            [~,panelnodes] = ismember(flipud(stringloc),cross.yzlocal{i}{1}(end/2+1:end,1));
            % Shift panelnodes to bottom skin
            panelnodes = panelnodes + size(cross.yzlocal{i}{1},1)/2;
            startpan = size(panels.yzlocal{i}{1},1);
            panels.yzlocal{i}{1} = [panels.yzlocal{i}{1};cross.yzlocal{i}{1}(panelnodes,:)];
            panels.elmloc{i}{1} = [panels.elmloc{i}{1};startpan+[(1:length(panelnodes)-1)',(2:length(panelnodes))']];
            % Identify element number and corresponding laminate by looking for elm where first node corresponds to the first panelnode
            [~,elmnum] = ismember(panelnodes(1:end-1),cross.elmloc{i}{1}(:,1));
            panels.lam{i}{1} = [panels.lam{i}{1};cross.lam{i}{1}(elmnum)];
            panels.type{i}{1} = [panels.type{i}{1};2*ones(size(cross.lam{i}{1}(elmnum)))]; % 2: Bottom skin
            % Identify all cross-sectional elements that belong to a certain
            % panel
            for j=1:length(panelnodes)-1
                elm1 = find(cross.elmloc{i}{1}(:,1)==panelnodes(j));
                elm2 = find(cross.elmloc{i}{1}(:,2)==panelnodes(j+1));
                panels.elmnum{i}{1}{j+startpan-1} = (elm1:elm2)';
            end
        else
            if FLAGSPAR
                if Nlam.Chord(lam_sec) == 1
                    [~,panelnodes] = ismember(flipud(skinloc'),cross.yzlocal{i}{1}(end/2+1:end,1));
                    % Shift panelnodes to bottom skin
                    panelnodes = panelnodes + size(cross.yzlocal{i}{1},1)/2;
                    startpan = size(panels.yzlocal{i}{1},1);
                    panels.yzlocal{i}{1} = [panels.yzlocal{i}{1};cross.yzlocal{i}{1}(panelnodes,:)];
                    panels.elmloc{i}{1} = [panels.elmloc{i}{1};startpan+[(1:length(panelnodes)-1)',(2:length(panelnodes))']];
                    % Identify element number and corresponding laminate by looking for elm where first node corresponds to the first panelnode
                    [~,elmnum] = ismember(panelnodes(1:end-1),cross.elmloc{i}{1}(:,1));
                    panels.lam{i}{1} = [panels.lam{i}{1};cross.lam{i}{1}(elmnum)];
                    panels.type{i}{1} = [panels.type{i}{1};2*ones(size(cross.lam{i}{1}(elmnum)))]; % 2: Bottom skin
                    % Identify all cross-sectional elements that belong to a certain
                    % panel
                    for j=1:length(panelnodes)-1
                        elm1 = find(cross.elmloc{i}{1}(:,1)==panelnodes(j));
                        elm2 = find(cross.elmloc{i}{1}(:,2)==panelnodes(j+1));
                        panels.elmnum{i}{1}{j+startpan-1} = (elm1:elm2)';
                    end
                else
                    [~,panelnodes] = ismember(flipud(lamloc'),cross.yzlocal{i}{1}(end/2+1:end,1));
                    % Shift panelnodes to bottom skin
                    panelnodes = panelnodes + size(cross.yzlocal{i}{1},1)/2;
                    startpan = size(panels.yzlocal{i}{1},1);
                    panels.yzlocal{i}{1} = [panels.yzlocal{i}{1};cross.yzlocal{i}{1}(panelnodes,:)];
                    panels.elmloc{i}{1} = [panels.elmloc{i}{1};startpan+[(1:length(panelnodes)-1)',(2:length(panelnodes))']];
                    % Identify element number and corresponding laminate by looking for elm where first node corresponds to the first panelnode
                    [~,elmnum] = ismember(panelnodes(1:end-1),cross.elmloc{i}{1}(:,1));
                    panels.lam{i}{1} = [panels.lam{i}{1};cross.lam{i}{1}(elmnum)];
                    panels.type{i}{1} = [panels.type{i}{1};2*ones(size(cross.lam{i}{1}(elmnum)))]; % 2: Bottom skin
                    % Identify all cross-sectional elements that belong to a certain
                    % panel
                    for j=1:length(panelnodes)-1
                        elm1 = find(cross.elmloc{i}{1}(:,1)==panelnodes(j));
                        elm2 = find(cross.elmloc{i}{1}(:,2)==panelnodes(j+1));
                        panels.elmnum{i}{1}{j+startpan-1} = (elm1:elm2)';
                    end
                end
            end
        end
        
        % Spars
        if FLAGSPAR
            for ispar = 1:nspars
                node1 = find(cross.yzlocal{i}{1}(:,1)==sparloc(ispar),1,'first');
                node2 = find(cross.yzlocal{i}{1}(:,1)==sparloc(ispar),1,'last');
                cross.elmID{i}{1} = [cross.elmID{i}{1};countelm+1];
                countelm = countelm + 1;
                cross.elmloc{i}{1} = [cross.elmloc{i}{1};[node1,node2]];
                cross.lam{i}{1} = [cross.lam{i}{1};constant.lam.SparID{lam_sec}(ispar)];
                
                startelm = size(panels.elmloc{i}{1},1);
                node1 = find(panels.yzlocal{i}{1}(:,1)==sparloc(ispar),1,'first');
                node2 = find(panels.yzlocal{i}{1}(:,1)==sparloc(ispar),1,'last');
                panels.elmloc{i}{1} = [panels.elmloc{i}{1};[node1,node2]];
                panels.lam{i}{1} = [panels.lam{i}{1};constant.lam.SparID{lam_sec}(ispar)];
                panels.type{i}{1} = [panels.type{i}{1};3*ones(size(cross.lam{i}{1}(elmnum)))]; % 3: Spars
                panels.elmnum{i}{1}{startelm+1} = countelm;
            end
        end
        
        % Define cross element sections (1: skin+spars, 2: stringer)
        cross.type{i}{1} = ones(size(cross.elmloc{i}{1},1),1);
        
        % Add stringers
        if FLAGSTRING
            if FLAGSPAR
                [~,sparnum] = ismember(skinloc,stringloc);
                stringloc(sparnum) = [];
            else
                stringloc = stringloc(2:end-1); % Remove leading edge and trailing edge stringer
            end
            % Number of stringers per skin
            nstring = length(stringloc);
            % Number of nodes without stringers
            nloc = size(cross.yzlocal{i}{1},1);
            
            % Create stringer nodes
            boxheight = max(WB_Z_Upper{i}-WB_Z_Lower{i});
            stringheight = stiff_height_sec*boxheight;
            [~,stringloc_box_upper] = ismember(stringloc,WB_X{i});
            [~,stringloc_box_lower] = ismember(stringloc,flipud(WB_X{i}));
            string_upper = WB_Z_Upper{i}(stringloc_box_upper)-stringheight;
            string_lower = WB_Z_Lower{i}(stringloc_box_lower)+stringheight;
            
            % Define stringer elements
            clear node1_top node1_bot
            for istring = 1:nstring
                node1_top(istring,1) =  find((cross.yzlocal{i}{1}(:,2)-WB_Z_Upper{i}(stringloc_box_upper(istring))).^2+(cross.yzlocal{i}{1}(:,1)-stringloc(istring)).^2==0);   % Top skin
                node1_bot(istring,1) =  find((cross.yzlocal{i}{1}(:,2)-WB_Z_Lower{i}(stringloc_box_lower(istring))).^2+(cross.yzlocal{i}{1}(:,1)-stringloc(istring)).^2==0);   % Bottom skin
            end
            if nstring == 0
                node1_top = [];
                node1_bot = [];
            end
            node2_top = nloc+(1:nstring)';
            node2_bot = nloc+nstring+(1:nstring)';
            
            % Add stringer nodes to cross-section
            cross.yzlocal{i}{1} = [cross.yzlocal{i}{1};[stringloc,string_upper;stringloc,string_lower]];
            
            % Add stringer elements to cross-section
            cross.elmID{i}{1} = [cross.elmID{i}{1};countelm+(1:2*nstring)'];
            cross.elmloc{i}{1} = [cross.elmloc{i}{1};[node1_top,node2_top;node1_bot,node2_bot]];
            
            % Add stringer laminate
            countstringer = countstringer+1;
            cross.lam{i}{1} = [cross.lam{i}{1};(length(constant.lam.ID)+countstringer)*ones(2*nstring,1)];
            
            constant.stringer.lamID(countstringer,1)=length(constant.lam.ID)+countstringer;
            constant.stringer.matID(countstringer,1) = 2;
            constant.stringer.h(countstringer,1) = stringheight*c_sec;
            constant.stringer.EA(countstringer,1) = stiff_EA_sec;
            constant.stringer.mA(countstringer,1) = stiff_mA_sec;
            
            % Add stringer as type
            cross.type{i}{1}(countelm+(1:2*nstring)) = 2;
        end
        
        % Correct cross.yzlocal for the reference axis location
        cross.yzlocal{i}{1}(:,1) = cross.yzlocal{i}{1}(:,1)-xref_sec+xsref_sec;
        cross.yzlocal{i}{1} = cross.yzlocal{i}{1}*c_sec;
        
        % Correct panels.yzlocal for the reference axis location
        panels.yzlocal{i}{1}(:,1) = panels.yzlocal{i}{1}(:,1)-xref_sec+xsref_sec;
        panels.yzlocal{i}{1} = panels.yzlocal{i}{1}*c_sec;
        
        if 1 && GRAPHFLAG == 1    % (visual check)
            figure(1)
            hold on
            for iplot=1:size(cross.elmloc{i}{1},1)
                Rplot         = [cos(theta_sec) sin(theta_sec); -sin(theta_sec) cos(theta_sec)];
                crossplot = (Rplot*[cross.yzlocal{i}{1}(cross.elmloc{i}{1}(iplot,:),1),cross.yzlocal{i}{1}(cross.elmloc{i}{1}(iplot,:),2)]')';
                plot3(crossplot(:,1)+(str.xyz(3*(i-1)+1)+str.xyz(3*i+1))/2,(str.xyz(3*(i-1)+2)+str.xyz(3*i+2))/2*ones(1,2),crossplot(:,2)+(str.xyz(3*(i-1)+3)+str.xyz(3*i+3))/2,'-r')
            end
        end
        
        % Plot the cross-sections and their corresponding laminates
        if 0
            figure(3+i)
            hold on
            for iplot=1:size(cross.elmloc{i}{1},1)
                plot(cross.yzlocal{i}{1}(cross.elmloc{i}{1}(iplot,:),1),cross.yzlocal{i}{1}(cross.elmloc{i}{1}(iplot,:),2),'-xr')
            end
            xtext = (cross.yzlocal{i}{1}(cross.elmloc{i}{1}(:,1),1)+cross.yzlocal{i}{1}(cross.elmloc{i}{1}(:,2),1))/2;
            ztext = (cross.yzlocal{i}{1}(cross.elmloc{i}{1}(:,1),2)+cross.yzlocal{i}{1}(cross.elmloc{i}{1}(:,2),2))/2;
            for itext = 1:length(xtext)
                text(xtext(itext),ztext(itext),num2str(cross.lam{i}{1}(itext)))
                %             text(xtext(itext),ztext(itext),num2str((itext)))
            end
        end
        
        % Plot the buckling panels and their corresponding laminates
        if 0
            figure(3+i)
            hold on
            for iplot=1:size(panels.elmloc{i}{1},1)
                plot(panels.yzlocal{i}{1}(panels.elmloc{i}{1}(iplot,:),1),panels.yzlocal{i}{1}(panels.elmloc{i}{1}(iplot,:),2),'-xb')
            end
            xtext = (panels.yzlocal{i}{1}(panels.elmloc{i}{1}(:,1),1)+panels.yzlocal{i}{1}(panels.elmloc{i}{1}(:,2),1))/2;
            ztext = (panels.yzlocal{i}{1}(panels.elmloc{i}{1}(:,1),2)+panels.yzlocal{i}{1}(panels.elmloc{i}{1}(:,2),2))/2;
            for itext = 1:length(xtext)
                text(xtext(itext),ztext(itext),num2str(panels.lam{i}{1}(itext)))
            end
        end
        
        % Rotate cross-section back to chord axis
        cross.yzlocal{i}{1} = (R'*cross.yzlocal{i}{1}')';
        
        % Correct for beam axis sweep angle
        % Calculate intersection
        X  = - str.xyz(1:3:end);
        Y  =   str.xyz(2:3:end);
        YX = [Y,X];
        % Relevant points
        B1  = [Y(i)  ,-wb_le_str(i)];
        B2  = [Y(i+1),-wb_le_str(i+1)];
        B3  = [Y(i+1),-wb_te_str(i+1)];
        B4  = [Y(i)  ,-wb_te_str(i)];
        B12 = 0.5*(B1 + B2);
        B43 = 0.5*(B4 + B3);
        B14 = YX(i,:);
        B23 = YX(i+1,:);
        P   = 0.5*(B14 + B23);
        % Calculate line slope
        m1  = (B2(2) - B1(2))/(B2(1) - B1(1));
        m2  = (B4(2) - B3(2))/(B4(1) - B3(1));
        m3  = (B23(2) - B14(2))/(B23(1) - B14(1));
        if m3 == 0
            xI1 = 0.5*(B1(1)+B2(1));
            yI1 = 0.5*(B1(2)+B2(2));
            
            xI2 = 0.5*(B3(1)+B4(1));
            yI2 = 0.5*(B3(2)+B4(2));
        else
            ms  = -1/m3;
            % 1st intersection
            xI1 = (m1*B1(1) - ms*P(1) - B1(2) + P(2))/(m1 - ms);
            yI1 = ms*(xI1 - P(1)) + P(2);
            % 2nd intersection
            xI2 = (m2*B3(1) - ms*P(1) - B3(2) + P(2))/(m2 - ms);
            yI2 = ms*(xI2 - P(1)) + P(2);
            % Length element
            %         l   = sqrt((xI1 - xI2)^2 + (yI1 - yI2)^2);
            %         l0  = norm(B43 - B12);
            %         % Scale cross-sectional length (X coord only)
            %         cross.yzlocal{i}{1}(:,1) = (l/l0)*cross.yzlocal{i}{1}(:,1);
        end
        
        l1   = sqrt((xI1 - P(1))^2 + (yI1 - P(2))^2);
        l01  = norm(P - B12);
        % Scale cross-sectional length (X coord only)
        cross.yzlocal{i}{1}(cross.yzlocal{i}{1}(:,1)<0,1) = (l1/l01)*cross.yzlocal{i}{1}(cross.yzlocal{i}{1}(:,1)<0,1);
        
        l2   = sqrt((P(1) - xI2)^2 + (P(2) - yI2)^2);
        l02  = norm(B43 - P);
        % Scale cross-sectional length (X coord only)
        cross.yzlocal{i}{1}(cross.yzlocal{i}{1}(:,1)>0,1) = (l2/l02)*cross.yzlocal{i}{1}(cross.yzlocal{i}{1}(:,1)>0,1);
        
        
%         sweep = atan(str.e1(3*(i-1)+1)/str.e1(3*(i-1)+2));              % beam axis sweep angle (deg)
%         cross.yzlocal{i}{1}(:,1) = cross.yzlocal{i}{1}(:,1)*cos(sweep);
        cross.yzlocal{i}{1}(:,2) = -cross.yzlocal{i}{1}(:,2); % Flip to convert cross-section to local coordinates (3 is down)
        
        
        %%% Buckling panels
        
        % Rotate cross-section back to chord axis
        panels.yzlocal{i}{1} = (R'*panels.yzlocal{i}{1}')';
        
        % Correct for beam axis sweep angle
%         panels.yzlocal{i}{1}(:,1) = panels.yzlocal{i}{1}(:,1)*cos(sweep);
        panels.yzlocal{i}{1}(panels.yzlocal{i}{1}(:,1)<0,1) = (l1/l01)*panels.yzlocal{i}{1}(panels.yzlocal{i}{1}(:,1)<0,1);
        panels.yzlocal{i}{1}(panels.yzlocal{i}{1}(:,1)>0,1) = (l2/l02)*panels.yzlocal{i}{1}(panels.yzlocal{i}{1}(:,1)>0,1);
        panels.yzlocal{i}{1}(:,2) = -panels.yzlocal{i}{1}(:,2); % Flip to convert cross-section to local coordinates (3 is down)
    end
end

% Include double section
if spanflag == 1
    for i=1:length(doublesection)
        cross.yzlocal{doublesection(i)}{1} = cross.yzlocal{inbsection(1)}{1};
        cross.yzlocal{doublesection(i)}{2} = cross.yzlocal{outbsection(1)}{1};
        cross.elmID{doublesection(i)}{1} = cross.elmID{inbsection(1)}{1};
        cross.elmID{doublesection(i)}{2} = cross.elmID{outbsection(1)}{1};
        cross.elmloc{doublesection(i)}{1} = cross.elmloc{inbsection(1)}{1};
        cross.elmloc{doublesection(i)}{2} = cross.elmloc{outbsection(1)}{1};
        cross.lam{doublesection(i)}{1} = cross.lam{inbsection(1)}{1};
        cross.lam{doublesection(i)}{2} = cross.lam{outbsection(1)}{1};
        cross.type{doublesection(i)}{1} = cross.type{inbsection(1)}{1};
        cross.type{doublesection(i)}{2} = cross.type{outbsection(1)}{1};
        
        panels.yzlocal{doublesection(i)}{1} = panels.yzlocal{inbsection(1)}{1};
        panels.yzlocal{doublesection(i)}{2} = panels.yzlocal{outbsection(1)}{1};
        panels.elmloc{doublesection(i)}{1} = panels.elmloc{inbsection(1)}{1};
        panels.elmloc{doublesection(i)}{2} = panels.elmloc{outbsection(1)}{1};
        panels.elmnum{doublesection(i)}{1} = panels.elmnum{inbsection(1)}{1};
        panels.elmnum{doublesection(i)}{2} = panels.elmnum{outbsection(1)}{1};
        panels.lam{doublesection(i)}{1} = panels.lam{inbsection(1)}{1};
        panels.lam{doublesection(i)}{2} = panels.lam{outbsection(1)}{1};
        panels.type{doublesection(i)}{1} = panels.type{inbsection(1)}{1};
        panels.type{doublesection(i)}{2} = panels.type{outbsection(1)}{1};
    end
end

% Pre-process the cross-section for the cross-sectional modeller
[cross,panels.elmconv] = prepr_cross(cross);