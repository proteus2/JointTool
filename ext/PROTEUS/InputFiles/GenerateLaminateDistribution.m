function [constant,Nlam] = GenerateLaminateDistribution(wing_data,constant,xyz,theta,c,Nlam,FLAGSPAR,varargin)

if length(varargin) == 2
    morph = varargin{1};
    fixnodes = varargin{2};
    morphflag = 1;
else
    morphflag = 0;
end

wb_le_right = wing_data(:,9);    % in % of chord
wb_le_left  = wb_le_right(2:end);
wb_le_vec   = [flipud(wb_le_left); wb_le_right];

wb_te_right = wing_data(:,10);   % in % of chord
wb_te_left  = wb_te_right(2:end);
wb_te_vec   = [flipud(wb_te_left); wb_te_right];

sparloc_right = wing_data(:,12:end);
sparloc_left  = sparloc_right(2:end,:);
sparloc_vec   = [flipud(sparloc_left); sparloc_right];

if morphflag == 1 && morph.inp.span.flag == 1
        
    xyz_fix = xyz(fixnodes==1,:);
        
    counttop = 0;
    countbot = 0;
    nlamtop = sum(Nlam.Chord)-1; % Subtract one to account for double section that will not get a unique laminate
    for i=1:size(Nlam.xyz_lam,1)-1
        if FLAGSPAR
            sparloc(i,:)    = interp1(xyz(:,2),wing_data(:,12:end),(Nlam.xyz_lam(i,2)+Nlam.xyz_lam(i+1,2))/2);
            Nlam.Nspars(i)  = length(sparloc(i,:))-sum(isnan(sparloc(i,:)));
            Nlam.SparNum{i} = find(~isnan(sparloc(i,:)));
        end
        if Nlam.xyz_lam(i,2)<xyz_fix(morph.inp.span.doublesec,2) || Nlam.xyz_lam(i,2)>xyz_fix(morph.inp.span.doublesec,2)
            Nlam.TopID{i} = counttop+(1:Nlam.Chord(i))';
            Nlam.BotID{i} = nlamtop+countbot+(1:Nlam.Chord(i))';
            counttop = counttop+Nlam.Chord(i);
            countbot = countbot+Nlam.Chord(i);
        else
            Nlam.TopID{i} = [];
            Nlam.BotID{i} = [];
        end
    end
else
    counttop = 0;
    countbot = 0;
    nlamtop = sum(Nlam.Chord);
    for i=1:size(Nlam.xyz_lam,1)-1
        if FLAGSPAR
%             sparloc(i,:)    = interp1(xyz(:,2),wing_data(:,12:end),(Nlam.xyz_lam(i,2)+Nlam.xyz_lam(i+1,2))/2);
            sparloc(i,:)    = interp1(xyz(:,2),sparloc_vec,(Nlam.xyz_lam(i,2)+Nlam.xyz_lam(i+1,2))/2);
            Nlam.Nspars(i)  = length(sparloc(i,:))-sum(isnan(sparloc(i,:)));
            Nlam.SparNum{i} = find(~isnan(sparloc(i,:)));
        end
        Nlam.TopID{i} = counttop+(1:Nlam.Chord(i))';
        Nlam.BotID{i} = nlamtop+countbot+(1:Nlam.Chord(i))';
        counttop = counttop+Nlam.Chord(i);
        countbot = countbot+Nlam.Chord(i);
    end
end


if FLAGSPAR
    if morphflag == 1 && morph.inp.span.flag == 1
        countspar = 0;
        for i=1:size(sparloc,2)
            for j=1:size(Nlam.xyz_lam,1)-1
                if ~isnan(sparloc(j,i))
                    if ~isfield(Nlam, 'SparID') || length(Nlam.SparID)<j
                        if Nlam.xyz_lam(j,2)<xyz_fix(morph.inp.span.doublesec,2) || Nlam.xyz_lam(j,2)>xyz_fix(morph.inp.span.doublesec,2)
                            Nlam.SparID{j}(1,1) = 2*nlamtop + countspar+1;
                            countspar = countspar+1;
                        else
                            Nlam.SparID{j} = [];
                        end
                    else
                        if Nlam.xyz_lam(j,2)<xyz_fix(morph.inp.span.doublesec,2) || Nlam.xyz_lam(j,2)>xyz_fix(morph.inp.span.doublesec,2)
                            Nlam.SparID{j}(end+1,1) = 2*nlamtop + countspar+1;
                            countspar = countspar+1;
                        else
                            Nlam.SparID{j} = [];
                        end
                    end
                    
                end
            end
        end
    else
        countspar = 0;
        for i=1:size(sparloc,2)
            for j=1:size(Nlam.xyz_lam,1)-1
                if ~isnan(sparloc(j,i))
                    if ~isfield(Nlam, 'SparID') || length(Nlam.SparID)<j
                        Nlam.SparID{j}(1,1) = 2*nlamtop + countspar+1;
                    else
                        Nlam.SparID{j}(end+1,1) = 2*nlamtop + countspar+1;
                    end
                    countspar = countspar+1;
                end
            end
        end
    end
end

if FLAGSPAR
    Nlam.Total = counttop+countbot+countspar;
else
    Nlam.Total = counttop+countbot;
end


% Define the laminate locations
for i=1:size(Nlam.xyz_lam,1)-1
    if morphflag == 1 && morph.inp.span.flag == 1
        if Nlam.xyz_lam(i,2)<xyz_fix(morph.inp.span.doublesec,2) || Nlam.xyz_lam(i,2)>xyz_fix(morph.inp.span.doublesec,2)
            secflag = 1;
        else
            secflag = 0;
        end
    else
        secflag = 1;
    end
    if FLAGSPAR && secflag == 1
        
        ID1 = find(xyz(:,2)==Nlam.xyz_lam(i,2));
        ID2 = find(xyz(:,2)==Nlam.xyz_lam(i+1,2));
        
        ylam1 = Nlam.xyz_lam(i,2);
        ylam2 = Nlam.xyz_lam(i+1,2);
        
        if isempty(ID1)
%             sparloc1 = interp1(xyz(:,2),wing_data(:,12:end),ylam1);
            sparloc1 = interp1(xyz(:,2),sparloc_vec,ylam1);
%             wb_le1  = interp1(xyz(:,2),wing_data(:,9),ylam1);
            wb_le1  = interp1(xyz(:,2),wb_le_vec,ylam1);
%             wb_te1  = interp1(xyz(:,2),wing_data(:,10),ylam1);
            wb_te1  = interp1(xyz(:,2),wb_te_vec,ylam1);
            c1 = interp1(xyz(:,2),c,ylam1);
            theta1 = interp1(xyz(:,2),theta,ylam1);
            wingLEx1 = interp1(xyz(:,2),constant.Coord3D.Wing_LE_xyz(:,1),ylam1);
            IndexProfile1 = find(ylam1>=xyz(:,2),1,'last');
            IndexProfile2 = find(ylam1<=xyz(:,2),1);
            
            TopSkin1 = (constant.Coord3D.TopSkin_xyz{IndexProfile2} - constant.Coord3D.TopSkin_xyz{IndexProfile1})...
                /(xyz(IndexProfile2,2)-xyz(IndexProfile1,2)+eps)  * (ylam1-xyz(IndexProfile1,2)) + constant.Coord3D.TopSkin_xyz{IndexProfile1};
            BotSkin1 = (constant.Coord3D.BotSkin_xyz{IndexProfile2} - constant.Coord3D.BotSkin_xyz{IndexProfile1})...
                /(xyz(IndexProfile2,2)-xyz(IndexProfile1,2)+eps)  * (ylam1-xyz(IndexProfile1,2)) + constant.Coord3D.BotSkin_xyz{IndexProfile1};
        else
%             sparloc1 = wing_data(ID1,12:end);
            sparloc1 = sparloc_vec(ID1,:);
%             wb_le1  = wing_data(ID1,9);
            wb_le1  = wb_le_vec(ID1);
%             wb_te1  = wing_data(ID1,10);
            wb_te1  = wb_te_vec(ID1);
            c1 = c(ID1);
            theta1 = theta(ID1);
            wingLEx1 = constant.Coord3D.Wing_LE_xyz(ID1,1);
            TopSkin1 = constant.Coord3D.TopSkin_xyz{ID1};
            BotSkin1 = constant.Coord3D.BotSkin_xyz{ID1};
        end
        
        if isempty(ID2)
%             sparloc2 = interp1(xyz(:,2),wing_data(:,12:end),ylam2);
            sparloc2 = interp1(xyz(:,2),sparloc_vec,ylam2);
%             wb_le2   = interp1(xyz(:,2),wing_data(:,9),ylam2);
            wb_le2   = interp1(xyz(:,2),wb_le_vec,ylam2);
%             wb_te2   = interp1(xyz(:,2),wing_data(:,10),ylam2);
            wb_te2   = interp1(xyz(:,2),wb_te_vec,ylam2);
            c2       = interp1(xyz(:,2),c,ylam2);
            theta2   = interp1(xyz(:,2),theta,ylam2);
            wingLEx2 = interp1(xyz(:,2),constant.Coord3D.Wing_LE_xyz(:,1),ylam2);
            IndexProfile1 = find(ylam2>=xyz(:,2),1,'last');
            IndexProfile2 = find(ylam2<=xyz(:,2),1);
            
            TopSkin2 = (constant.Coord3D.TopSkin_xyz{IndexProfile2} - constant.Coord3D.TopSkin_xyz{IndexProfile1})...
                /(xyz(IndexProfile2,2)-xyz(IndexProfile1,2)+eps)  * (ylam2-xyz(IndexProfile1,2)) + constant.Coord3D.TopSkin_xyz{IndexProfile1};
            BotSkin2 = (constant.Coord3D.BotSkin_xyz{IndexProfile2} - constant.Coord3D.BotSkin_xyz{IndexProfile1})...
                /(xyz(IndexProfile2,2)-xyz(IndexProfile1,2)+eps)  * (ylam2-xyz(IndexProfile1,2)) + constant.Coord3D.BotSkin_xyz{IndexProfile1};
            
            
        else
%             sparloc2 = wing_data(ID2,12:end);
            sparloc2 = sparloc_vec(ID2,:);
%             wb_le2  = wing_data(ID2,9);
            wb_le2  = wb_le_vec(ID2);
%             wb_te2  = wing_data(ID2,10);
            wb_te2  = wb_te_vec(ID2);
            c2 = c(ID2);
            theta2 = theta(ID2);
            wingLEx2 = constant.Coord3D.Wing_LE_xyz(ID2,1);
            TopSkin2 = constant.Coord3D.TopSkin_xyz{ID2};
            BotSkin2 = constant.Coord3D.BotSkin_xyz{ID2};
        end
        
        ind = unique([find(isnan(sparloc1)),find(isnan(sparloc2))]);
        sparnum       = find(ones(size(sparloc1)));
        sparnum(ind)  = [];
        sparloc1(ind) = [];
        sparloc2(ind) = [];
        
        if Nlam.Chord(i) == 1
            lamloc1 = [wb_le1,wb_te1];
            lamloc2 = [wb_le2,wb_te2];
        else
            skinloc1 = unique([wb_le1,sparloc1,wb_te1]);
            skinloc2 = unique([wb_le2,sparloc2,wb_te2]);
            if mod(Nlam.Chord(i),(length(skinloc1)-1))~=0
                error(['Please specify a number of chordwise laminates equal to 1 or a multiple of the number of laminate bays for section ',num2str(i)])
            else
                lamloc1 = [];
                lamloc2 = [];
                for ispar = 1:length(skinloc1)-1
                    lamloc1 = [lamloc1(1:end-1),linspace(skinloc1(ispar),skinloc1(ispar+1),(Nlam.Chord(i)/(length(skinloc1)-1))+1)];
                    lamloc2 = [lamloc2(1:end-1),linspace(skinloc2(ispar),skinloc2(ispar+1),(Nlam.Chord(i)/(length(skinloc2)-1))+1)];
                end
            end
        end
        
        lamloc1denorm = lamloc1*c1*cos(theta1)+wingLEx1;
        lamloc2denorm = lamloc2*c2*cos(theta2)+wingLEx2;
        
        
        % Top skin
        for iskin=1:Nlam.Chord(i)
            if isnan(interp1(TopSkin1(:,1),TopSkin1(:,3),lamloc1denorm(iskin)))
                loc1 = [interp1(TopSkin1(:,1),TopSkin1(:,1),lamloc1denorm(iskin),'nearest','extrap'),ylam1,interp1(TopSkin1(:,1),TopSkin1(:,3),lamloc1denorm(iskin),'nearest','extrap')];
            else
                loc1 = [lamloc1denorm(iskin),ylam1,interp1(TopSkin1(:,1),TopSkin1(:,3),lamloc1denorm(iskin))];
            end
            if isnan(interp1(TopSkin1(:,1),TopSkin1(:,3),lamloc1denorm(iskin+1)))
                loc2 = [interp1(TopSkin1(:,1),TopSkin1(:,1),lamloc1denorm(iskin+1),'nearest','extrap'),ylam1,interp1(TopSkin1(:,1),TopSkin1(:,3),lamloc1denorm(iskin+1),'nearest','extrap')];
            else
                loc2 = [lamloc1denorm(iskin+1),ylam1,interp1(TopSkin1(:,1),TopSkin1(:,3),lamloc1denorm(iskin+1))];
            end
            if isnan(interp1(TopSkin2(:,1),TopSkin2(:,3),lamloc2denorm(iskin+1)))
                loc3 = [interp1(TopSkin2(:,1),TopSkin2(:,1),lamloc2denorm(iskin+1),'nearest','extrap'),ylam2,interp1(TopSkin2(:,1),TopSkin2(:,3),lamloc2denorm(iskin+1),'nearest','extrap')];
            else
                loc3 = [lamloc2denorm(iskin+1),ylam2,interp1(TopSkin2(:,1),TopSkin2(:,3),lamloc2denorm(iskin+1))];
            end
            if isnan(interp1(TopSkin2(:,1),TopSkin2(:,3),lamloc2denorm(iskin)))
                loc4 = [interp1(TopSkin2(:,1),TopSkin2(:,1),lamloc2denorm(iskin),'nearest','extrap'),ylam2,interp1(TopSkin2(:,1),TopSkin2(:,3),lamloc2denorm(iskin),'nearest','extrap')];
            else
                loc4 = [lamloc2denorm(iskin),ylam2,interp1(TopSkin2(:,1),TopSkin2(:,3),lamloc2denorm(iskin))];
            end
            constant.lam.coord{Nlam.TopID{i}(iskin)} = [loc1;loc2;loc3;loc4];

            if 0
                figure(10)
                hold on
                fill3(constant.lam.coord{Nlam.TopID{i}(iskin)}(:,1),...
                    constant.lam.coord{Nlam.TopID{i}(iskin)}(:,2),...
                    constant.lam.coord{Nlam.TopID{i}(iskin)}(:,3), [1 0.5 1])
                
                text (mean(constant.lam.coord{Nlam.TopID{i}(iskin)}(:,1)),...
                    mean(constant.lam.coord{Nlam.TopID{i}(iskin)}(:,2)),...
                    mean(constant.lam.coord{Nlam.TopID{i}(iskin)}(:,3)) + 0.1,num2str(Nlam.TopID{i}(iskin)));
            end
        end
        
        % Bottom skin
        for iskin=1:Nlam.Chord(i)
            if isnan(interp1(BotSkin1(:,1),BotSkin1(:,3),lamloc1denorm(iskin)))
                loc1 = [interp1(BotSkin1(:,1),BotSkin1(:,1),lamloc1denorm(iskin),'nearest','extrap'),ylam1,interp1(BotSkin1(:,1),BotSkin1(:,3),lamloc1denorm(iskin),'nearest','extrap')];
            else
                loc1 = [lamloc1denorm(iskin),ylam1,interp1(BotSkin1(:,1),BotSkin1(:,3),lamloc1denorm(iskin))];
            end
            
            if isnan(interp1(BotSkin1(:,1),BotSkin1(:,3),lamloc1denorm(iskin+1)))
                loc2 = [interp1(BotSkin1(:,1),BotSkin1(:,1),lamloc1denorm(iskin+1),'nearest','extrap'),ylam1,interp1(BotSkin1(:,1),BotSkin1(:,3),lamloc1denorm(iskin+1),'nearest','extrap')];
            else
                loc2 = [lamloc1denorm(iskin+1),ylam1,interp1(BotSkin1(:,1),BotSkin1(:,3),lamloc1denorm(iskin+1))];
            end
            
            if isnan(interp1(BotSkin2(:,1),BotSkin2(:,3),lamloc2denorm(iskin+1)))
                loc3 = [interp1(BotSkin2(:,1),BotSkin2(:,1),lamloc2denorm(iskin+1),'nearest','extrap'),ylam2,interp1(BotSkin2(:,1),BotSkin2(:,3),lamloc2denorm(iskin+1),'nearest','extrap')];
            else
                loc3 = [lamloc2denorm(iskin+1),ylam2,interp1(BotSkin2(:,1),BotSkin2(:,3),lamloc2denorm(iskin+1))];
            end
            
            if isnan(interp1(BotSkin2(:,1),BotSkin2(:,3),lamloc2denorm(iskin)))
                loc4 = [interp1(BotSkin2(:,1),BotSkin2(:,1),lamloc2denorm(iskin),'nearest','extrap'),ylam2,interp1(BotSkin2(:,1),BotSkin2(:,3),lamloc2denorm(iskin),'nearest','extrap')];
            else
                loc4 = [lamloc2denorm(iskin),ylam2,interp1(BotSkin2(:,1),BotSkin2(:,3),lamloc2denorm(iskin))];
            end
            constant.lam.coord{Nlam.BotID{i}(iskin)} = [loc1;loc2;loc3;loc4];

            if 0
                figure(10)
                fill3(constant.lam.coord{Nlam.BotID{i}(iskin)}(:,1),...
                    constant.lam.coord{Nlam.BotID{i}(iskin)}(:,2),...
                    constant.lam.coord{Nlam.BotID{i}(iskin)}(:,3), [1 0.5 1])
                
                text (mean(constant.lam.coord{Nlam.BotID{i}(iskin)}(:,1)),...
                    mean(constant.lam.coord{Nlam.BotID{i}(iskin)}(:,2)),...
                    mean(constant.lam.coord{Nlam.BotID{i}(iskin)}(:,3)) - 0.1,num2str(Nlam.BotID{i}(iskin)));
                
            end
        end
     
            
        % Spars
        for ispar=1:length(sparnum)
            ID1 = find(constant.Coord3D.Spars_xyzTop{sparnum(ispar)}(:,2)==Nlam.xyz_lam(i,2));
            ID2 = find(constant.Coord3D.Spars_xyzTop{sparnum(ispar)}(:,2)==Nlam.xyz_lam(i+1,2));
            if isempty(ID1)
                loc1 = [interp1(constant.Coord3D.Spars_xyzTop{sparnum(ispar)}(:,2),constant.Coord3D.Spars_xyzTop{sparnum(ispar)}(:,1),ylam1),...
                        ylam1,...
                        interp1(constant.Coord3D.Spars_xyzTop{sparnum(ispar)}(:,2),constant.Coord3D.Spars_xyzTop{sparnum(ispar)}(:,3),ylam1)];
                loc2 = [interp1(constant.Coord3D.Spars_xyzBot{sparnum(ispar)}(:,2),constant.Coord3D.Spars_xyzBot{sparnum(ispar)}(:,1),ylam1),...
                        ylam1,...
                        interp1(constant.Coord3D.Spars_xyzBot{sparnum(ispar)}(:,2),constant.Coord3D.Spars_xyzBot{sparnum(ispar)}(:,3),ylam1)];
            else
                loc1 = constant.Coord3D.Spars_xyzTop{sparnum(ispar)}(ID1,:);
                loc2 = constant.Coord3D.Spars_xyzBot{sparnum(ispar)}(ID1,:);
            end
            if isempty(ID2)
                loc4 = [interp1(constant.Coord3D.Spars_xyzTop{sparnum(ispar)}(:,2),constant.Coord3D.Spars_xyzTop{sparnum(ispar)}(:,1),ylam2),...
                        ylam2,...
                        interp1(constant.Coord3D.Spars_xyzTop{sparnum(ispar)}(:,2),constant.Coord3D.Spars_xyzTop{sparnum(ispar)}(:,3),ylam2)];
                loc3 = [interp1(constant.Coord3D.Spars_xyzBot{sparnum(ispar)}(:,2),constant.Coord3D.Spars_xyzBot{sparnum(ispar)}(:,1),ylam2),...
                        ylam2,...
                        interp1(constant.Coord3D.Spars_xyzBot{sparnum(ispar)}(:,2),constant.Coord3D.Spars_xyzBot{sparnum(ispar)}(:,3),ylam2)];
            else
                loc4 = constant.Coord3D.Spars_xyzTop{sparnum(ispar)}(ID2,:);
                loc3 = constant.Coord3D.Spars_xyzBot{sparnum(ispar)}(ID2,:);
            end
            
                constant.lam.coord{Nlam.SparID{i}(ispar)} = [loc1;loc2;loc3;loc4];
            
            if 0
                figure(10)
                fill3(constant.lam.coord{Nlam.SparID{i}(ispar)}(:,1),...
                      constant.lam.coord{Nlam.SparID{i}(ispar)}(:,2),...
                      constant.lam.coord{Nlam.SparID{i}(ispar)}(:,3), [1 1 0])
                  
               text (mean(constant.lam.coord{Nlam.SparID{i}(ispar)}(:,1)) - 0.1,...
                    mean(constant.lam.coord{Nlam.SparID{i}(ispar)}(:,2)),...
                    mean(constant.lam.coord{Nlam.SparID{i}(ispar)}(:,3)),num2str(Nlam.SparID{i}(ispar)));
            end
        end
  
    elseif secflag == 1 % No spars
        ID1 = find(xyz(:,2)==Nlam.xyz_lam(i));
        ID2 = find(xyz(:,2)==Nlam.xyz_lam(i+1));
        
        ylam1 = Nlam.xyz_lam(i,2);
        ylam2 = Nlam.xyz_lam(i+1,2);
        
        if isempty(ID1)
            wb_le1  = interp1(xyz(:,2),wing_data(:,9),ylam1);
            wb_te1  = interp1(xyz(:,2),wing_data(:,10),ylam1);
            c1 = interp1(xyz(:,2),c,ylam1);
            theta1 = interp1(xyz(:,2),theta,ylam1);
            wingLEx1 = interp1(xyz(:,2),constant.Coord3D.Wing_LE_xyz(:,1),ylam1);
            IndexProfile1 = find(ylam1>=xyz(:,2),1,'last');
            IndexProfile2 = find(ylam1<=xyz(:,2),1);
            
            TopSkin1 = (constant.Coord3D.TopSkin_xyz{IndexProfile2} - constant.Coord3D.TopSkin_xyz{IndexProfile1})...
                /(xyz(IndexProfile2,2)-xyz(IndexProfile1,2)+eps)  * (ylam1-xyz(IndexProfile1,2)) + constant.Coord3D.TopSkin_xyz{IndexProfile1};
            BotSkin1 = (constant.Coord3D.BotSkin_xyz{IndexProfile2} - constant.Coord3D.BotSkin_xyz{IndexProfile1})...
                /(xyz(IndexProfile2,2)-xyz(IndexProfile1,2)+eps)  * (ylam1-xyz(IndexProfile1,2)) + constant.Coord3D.BotSkin_xyz{IndexProfile1};
        else
            wb_le1  = wing_data(ID1,9);
            wb_te1  = wing_data(ID1,10);
            c1 = c(ID1);
            theta1 = theta(ID1);
            wingLEx1 = constant.Coord3D.Wing_LE_xyz(ID1,1);
            TopSkin1 = constant.Coord3D.TopSkin_xyz{ID1};
            BotSkin1 = constant.Coord3D.BotSkin_xyz{ID1};
        end
        
        if isempty(ID2)
            wb_le2  = interp1(xyz(:,2),wing_data(:,9),ylam2);
            wb_te2  = interp1(xyz(:,2),wing_data(:,10),ylam2);
            c2 = interp1(xyz(:,2),c,ylam2);
            theta2 = interp1(xyz(:,2),theta,ylam2);
            wingLEx2 = interp1(xyz(:,2),constant.Coord3D.Wing_LE_xyz(:,1),ylam2);
            IndexProfile1 = find(ylam2>=xyz(:,2),1,'last');
            IndexProfile2 = find(ylam2<=xyz(:,2),1);
            
            TopSkin2 = (constant.Coord3D.TopSkin_xyz{IndexProfile2} - constant.Coord3D.TopSkin_xyz{IndexProfile1})...
                /(xyz(IndexProfile2,2)-xyz(IndexProfile1,2)+eps)  * (ylam2-xyz(IndexProfile1,2)) + constant.Coord3D.TopSkin_xyz{IndexProfile1};
            BotSkin2 = (constant.Coord3D.BotSkin_xyz{IndexProfile2} - constant.Coord3D.BotSkin_xyz{IndexProfile1})...
                /(xyz(IndexProfile2,2)-xyz(IndexProfile1,2)+eps)  * (ylam2-xyz(IndexProfile1,2)) + constant.Coord3D.BotSkin_xyz{IndexProfile1};
            
            
        else
            wb_le2  = wing_data(ID2,9);
            wb_te2  = wing_data(ID2,10);
            c2 = c(ID2);
            theta2 = theta(ID2);
            wingLEx2 = constant.Coord3D.Wing_LE_xyz(ID2,1);
            TopSkin2 = constant.Coord3D.TopSkin_xyz{ID2};
            BotSkin2 = constant.Coord3D.BotSkin_xyz{ID2};
        end
        
        if Nlam.Chord(i) == 1
            lamloc1 = [wb_le1,wb_te1];
            lamloc2 = [wb_le2,wb_te2];
        else
            lamloc1 = linspace(wb_le1,wb_te1,(Nlam.Chord(i))+1);
            lamloc2 = linspace(wb_le2,wb_te2,(Nlam.Chord(i))+1);
        end
        
        lamloc1denorm = lamloc1*c1*cos(theta1)+wingLEx1;
        lamloc2denorm = lamloc2*c2*cos(theta2)+wingLEx2;
        
        % Top skin
        for iskin=1:Nlam.Chord(i)
            loc1 = [lamloc1denorm(iskin),ylam1,interp1(TopSkin1(:,1),TopSkin1(:,3),lamloc1denorm(iskin))];
            loc2 = [lamloc1denorm(iskin+1),ylam1,interp1(TopSkin1(:,1),TopSkin1(:,3),lamloc1denorm(iskin+1))];
            loc3 = [lamloc2denorm(iskin+1),ylam2,interp1(TopSkin2(:,1),TopSkin2(:,3),lamloc2denorm(iskin+1))];
            loc4 = [lamloc2denorm(iskin),ylam2,interp1(TopSkin2(:,1),TopSkin2(:,3),lamloc2denorm(iskin))];
            constant.lam.coord{Nlam.TopID{i}(iskin)} = [loc1;loc2;loc3;loc4];

%             figure(1)
%             fill3(constant.lam.coord{Nlam.TopID{i}(iskin)}(:,1),...
%                     constant.lam.coord{Nlam.TopID{i}(iskin)}(:,2),...
%                     constant.lam.coord{Nlam.TopID{i}(iskin)}(:,3), [1 0.5 1])
        end
        
        % Bottom skin
        for iskin=1:Nlam.Chord(i)
            loc1 = [lamloc1denorm(iskin),ylam1,interp1(BotSkin1(:,1),BotSkin1(:,3),lamloc1denorm(iskin))];
            loc2 = [lamloc1denorm(iskin+1),ylam1,interp1(BotSkin1(:,1),BotSkin1(:,3),lamloc1denorm(iskin+1))];
            loc3 = [lamloc2denorm(iskin+1),ylam2,interp1(BotSkin2(:,1),BotSkin2(:,3),lamloc2denorm(iskin+1))];
            loc4 = [lamloc2denorm(iskin),ylam2,interp1(BotSkin2(:,1),BotSkin2(:,3),lamloc2denorm(iskin))];
            constant.lam.coord{Nlam.BotID{i}(iskin)} = [loc1;loc2;loc3;loc4];

%             figure(1)
%             fill3(constant.lam.coord{Nlam.BotID{i}(iskin)}(:,1),...
%                     constant.lam.coord{Nlam.BotID{i}(iskin)}(:,2),...
%                     constant.lam.coord{Nlam.BotID{i}(iskin)}(:,3), [1 0.5 1])
        end
    end
end

constant.lam.ID = (1:Nlam.Total)';
% Assumption: all laminates have the same material, can be changed
% according to the user needs
constant.lam.matID = ones(Nlam.Total,1);
constant.lam.xyzb = Nlam.xyz_lam;
constant.lam.TopID = Nlam.TopID;
constant.lam.BotID = Nlam.BotID;
if FLAGSPAR
    constant.lam.SparID = Nlam.SparID;
    constant.lam.SparNum = Nlam.SparNum;
end

