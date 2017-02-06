function [morph] = AddMorphing(xlsFileName,numlc)

[TwistMorphData,TwistMorphDataStr] = xlsread(xlsFileName,'Twist morphing');

% Twist (defined in degrees)
morph.inp.twist.sec = TwistMorphData(1,:); % Define morphable sections as defined by wing geometry; fixed nodes
if sum(morph.inp.twist.sec~=0)
    twistlocs = find(morph.inp.twist.sec==1);
    morph.inp.twist.angleini = zeros(sum(morph.inp.twist.sec),numlc);
    morph.inp.twist.anglefin = zeros(sum(morph.inp.twist.sec),numlc);
    morph.inp.twist.anglelcinp = zeros(sum(morph.inp.twist.sec),numlc);
    for i = 1:numlc
        angleini = TwistMorphData(4+(2*(i-1)+1),1:length(morph.inp.twist.sec));
        morph.inp.twist.angleini(:,i) = angleini(morph.inp.twist.sec==1)';
        for j=1:sum(morph.inp.twist.sec)
            if isnan(morph.inp.twist.angleini(j,i))
                morph.inp.twist.anglelcinp(j,i) = str2double(TwistMorphDataStr{5+(2*(i-1)+1),twistlocs(j)+1}(3));% Define from which loadcases the initial data needs to be read
            end
        end
        anglefin = TwistMorphData(5+(2*(i-1)+1),1:length(morph.inp.twist.sec));
        morph.inp.twist.anglefin(:,i) = anglefin(morph.inp.twist.sec==1)';
    end
    morph.inp.twist.low = TwistMorphData(2,morph.inp.twist.sec==1); % Twist lower limits
    morph.inp.twist.high = TwistMorphData(3,morph.inp.twist.sec==1);   % Twist upper limits
end

[ShearMorphData,ShearMorphDataStr] = xlsread(xlsFileName,'Shear morphing');

% Shear (defined in degrees)
morph.inp.shear.sec = ShearMorphData(1,:); % Define morphable sections as defined by wing geometry; fixed nodes
if sum(morph.inp.shear.sec~=0)
    shearlocs = find(morph.inp.shear.sec==1);
    morph.inp.shear.angleini = zeros(sum(morph.inp.shear.sec),numlc);
    morph.inp.shear.anglefin = zeros(sum(morph.inp.shear.sec),numlc);
    morph.inp.shear.anglelcinp = zeros(sum(morph.inp.shear.sec),numlc);
    for i = 1:numlc
        angleini = ShearMorphData(4+(2*(i-1)+1),1:length(morph.inp.shear.sec));
        morph.inp.shear.angleini(:,i) = angleini(morph.inp.shear.sec==1)';
        for j=1:sum(morph.inp.shear.sec)
            if isnan(morph.inp.shear.angleini(j,i))
                morph.inp.shear.anglelcinp(j,i) = str2double(ShearMorphDataStr{5+(2*(i-1)+1),shearlocs(j)+1}(3));% Define from which loadcases the initial data needs to be read
            end
        end
        anglefin = ShearMorphData(5+(2*(i-1)+1),1:length(morph.inp.shear.sec));
        morph.inp.shear.anglefin(:,i) = anglefin(morph.inp.shear.sec==1)';
    end
    morph.inp.shear.low = ShearMorphData(2,morph.inp.shear.sec==1); % Shear lower limits
    morph.inp.shear.high = ShearMorphData(3,morph.inp.shear.sec==1);   % Shear upper limits
end

[CamberMorphData,CamberMorphDataStr] = xlsread(xlsFileName,'Camber morphing');

% Camber
morph.inp.camber.loc = CamberMorphData(5,:); % Define morphable locations as defined by wing geometry; xyz
if sum(morph.inp.camber.loc~=0)
    morph.camber.nbcf = CamberMorphData(1,1);
    morph.camber.nbcdisp = CamberMorphData(2,1);
    
    camberlocs = find(morph.inp.camber.loc==1);
    morph.inp.camber.ini = zeros(1,length(morph.inp.camber.loc));
    morph.inp.camber.ini(camberlocs) = CamberMorphData(6,camberlocs);
    morph.inp.camber.axis = CamberMorphData(7,1:length(morph.inp.camber.loc)); % Axis about which camber occurs
    morph.inp.camber.paramini = zeros(sum(morph.inp.camber.loc),numlc);
    morph.inp.camber.paramfin = zeros(sum(morph.inp.camber.loc),numlc);
    morph.inp.camber.paramlcinp = zeros(sum(morph.inp.camber.loc),numlc);
    for i = 1:numlc
        paramini = CamberMorphData(10+(2*(i-1)+1),1:length(morph.inp.camber.loc));
        morph.inp.camber.paramini(:,i) = paramini(morph.inp.camber.loc==1)';
        for j=1:sum(morph.inp.camber.loc)
            if isnan(morph.inp.camber.paramini(j,i))
                morph.inp.camber.paramlcinp(j,i) = str2double(CamberMorphDataStr{11+(2*(i-1)+1),camberlocs(j)+1}(3));% Define from which loadcases the initial data needs to be read
            end
        end
        paramfin = CamberMorphData(11+(2*(i-1)+1),1:length(morph.inp.camber.loc));
        morph.inp.camber.paramfin(:,i) = paramfin(morph.inp.camber.loc==1)';
    end
    morph.inp.camber.low = CamberMorphData(8,morph.inp.camber.loc==1); % Camber lower limits
    morph.inp.camber.high = CamberMorphData(9,morph.inp.camber.loc==1);   % Camber upper limits
end

[SpanMorphData,SpanMorphDataStr] = xlsread(xlsFileName,'Span morphing');

% Span extension
morph.inp.span.flag = SpanMorphData(1,1);
if morph.inp.span.flag == 1
    morph.inp.span.fixedsec = SpanMorphData(2,1); % Define the inner fixed section for span extension
    morph.inp.span.lamfixed = SpanMorphData(2,3:end);
    morph.inp.span.lamfixed(isnan(morph.inp.span.lamfixed))=[];
    morph.inp.span.doublesec = SpanMorphData(3,1); % Define the overlap section for span extension
    morph.inp.span.extsec = SpanMorphData(4,1); % Define the overlap section for span extension
    morph.inp.span.lamext = SpanMorphData(4,3:end);
    morph.inp.span.lamext(isnan(morph.inp.span.lamext))=[];
    morph.inp.span.extlcinp = zeros(1,numlc);
    morph.inp.span.extlcinpfin = zeros(1,numlc);
    for i = 1:numlc
        try 
            morph.inp.span.extini(1,i) = SpanMorphData(7+(2*(i-1)+1),1);
        catch err
            morph.inp.span.extini(1,i) = NaN;
        end
        if isnan(morph.inp.span.extini(1,i))
            morph.inp.span.extlcinp(1,i) = str2double(SpanMorphDataStr{8+(2*(i-1)+1),2}(3));% Define from which loadcases the initial data needs to be read
            if strcmp(SpanMorphDataStr{8+(2*(i-1)+1),2},'ini')
                morph.inp.span.extlcinp(1,i) = -1; % -1 indicates that it should be replaced by the undeformed input
            end
        end
        try 
            morph.inp.span.extfin(1,i) = SpanMorphData(7+(2*(i-1)+2),1);
        catch err
            morph.inp.span.extfin(1,i) = NaN;
        end
        if isnan(morph.inp.span.extfin(1,i))
            morph.inp.span.extlcinpfin(1,i) = str2double(SpanMorphDataStr{8+(2*(i-1)+2),2}(3));% Define from which loadcases the initial data needs to be read
            if strcmp(SpanMorphDataStr{8+(2*(i-1)+2),2},'ini')
                morph.inp.span.extlcinpfin(1,i) = -1; % -1 indicates that it should be replaced by the undeformed input
            end
        end
    end
    morph.inp.span.low = SpanMorphData(5,1); % Span extension lower limits
    morph.inp.span.high = SpanMorphData(6,1);   % Span extension upper limits, defined as fraction of overlap section length (<=0.95, for structural and numerical stability)
end

[FoldMorphData,FoldMorphDataStr] = xlsread(xlsFileName,'Fold morphing');

% Fold (defined in degrees)
morph.inp.fold.sec = FoldMorphData(1,:); % Define morphable sections as defined by wing geometry; fixed nodes
if sum(morph.inp.fold.sec~=0)
    foldlocs = find(morph.inp.fold.sec==1);
    morph.inp.fold.angleini = zeros(sum(morph.inp.fold.sec),numlc);
    morph.inp.fold.anglefin = zeros(sum(morph.inp.fold.sec),numlc);
    morph.inp.fold.anglelcinp = zeros(sum(morph.inp.fold.sec),numlc);
    for i = 1:numlc
        angleini = FoldMorphData(4+(2*(i-1)+1),1:length(morph.inp.fold.sec));
        morph.inp.fold.angleini(:,i) = angleini(morph.inp.fold.sec==1)';
        for j=1:sum(morph.inp.fold.sec)
            if isnan(morph.inp.fold.angleini(j,i))
                morph.inp.fold.anglelcinp(j,i) = str2double(FoldMorphDataStr{5+(2*(i-1)+1),foldlocs(j)+1}(3));% Define from which loadcases the initial data needs to be read
            end
        end
        anglefin = FoldMorphData(5+(2*(i-1)+1),1:length(morph.inp.fold.sec));
        morph.inp.fold.anglefin(:,i) = anglefin(morph.inp.fold.sec==1)';
    end
    morph.inp.fold.low = FoldMorphData(2,morph.inp.fold.sec==1); % Fold lower limits
    morph.inp.fold.high = FoldMorphData(3,morph.inp.fold.sec==1);   % Fold upper limits
end
