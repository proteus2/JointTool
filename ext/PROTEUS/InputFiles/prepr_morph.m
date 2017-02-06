function [morph] = prepr_morph(str,morph)

% Twist
if sum(morph.inp.twist.sec)~=0
    morph.twist.sec = zeros(str.Ns,1);
%     morph.twist.angle = zeros(str.Ns,1);
    morph.twist.dangleddv = zeros(str.Ns,sum(morph.inp.twist.sec));
%     morph.twist.low = zeros(str.Ns,1);
%     morph.twist.high = zeros(str.Ns,1);
    twistlocs = find(morph.inp.twist.sec==1);
    
    
    for i=1:length(morph.inp.twist.sec)
        if morph.inp.twist.sec(i) == 1
            elm = find(str.sec==i);
            % Determine the overall beam length of a section
            xyznode1 = [str.xyz(3*(elm-1)+1),str.xyz(3*(elm-1)+2),str.xyz(3*(elm-1)+3)];
            xyznode2 = [str.xyz(3*(elm)+1),str.xyz(3*(elm)+2),str.xyz(3*(elm)+3)];
            elmlength = sqrt(sum((xyznode2-xyznode1).^2,2));
           
            morph.twist.sec(elm) = 1;
%             morph.twist.angle(elm) = deg2rad(morph.inp.twist.angle(i).*(elmlength/sum(elmlength)));
            morph.twist.dangleddv(elm,twistlocs==i) = (elmlength/sum(elmlength));
%             morph.twist.low(elm) = deg2rad(morph.inp.twist.low(i));
%             morph.twist.high(elm) = deg2rad(morph.inp.twist.high(i));
        end
    end
end

% Shear
if sum(morph.inp.shear.sec)~=0
    morph.shear.sec = zeros(str.Ns,1);
%     morph.shear.angle = zeros(str.Ns,1);
    morph.shear.dangleddv = zeros(str.Ns,sum(morph.inp.shear.sec));
%     morph.shear.low = zeros(str.Ns,1);
%     morph.shear.high = zeros(str.Ns,1);
    shearlocs = find(morph.inp.shear.sec==1);
    
    for i=1:length(morph.inp.shear.sec)
        if morph.inp.shear.sec(i) == 1
            elm = find(str.sec==i);
            morph.shear.sec(elm) = 1;
%             morph.shear.angle(elm) = deg2rad(morph.inp.shear.angle(i));
            morph.shear.dangleddv(elm,shearlocs==i) = 1;
%             morph.shear.low(elm) = deg2rad(morph.inp.shear.low(i));
%             morph.shear.high(elm) = deg2rad(morph.inp.shear.high(i));
        end
    end
end

% Span
if morph.inp.span.flag == 1
    morph.span.sec = zeros(str.Ns,1);
%     morph.span.ext = zeros(str.Ns,1);
    morph.span.dextddv = zeros(str.Ns,1);
%     morph.span.low = zeros(str.Ns,1);
%     morph.span.high = zeros(str.Ns,1);
    morph.span.lelmmin = zeros(str.Ns,1);
    morph.span.dl = zeros(str.Ns,1);
    morph.span.ext0 = zeros(str.Ns,1);
    
    elmf = find(str.sec==morph.inp.span.fixedsec);
    elmd = find(str.sec==morph.inp.span.doublesec);
    elme = find(str.sec==morph.inp.span.extsec);
    
    lelm = sqrt((str.xyz(1:3:end-3)-str.xyz(4:3:end)).^2+(str.xyz(2:3:end-3)-str.xyz(5:3:end)).^2+(str.xyz(3:3:end-3)-str.xyz(6:3:end)).^2);
    
    lfixed = sum(lelm([elmf;elmd])); % Total length of the fixed section
    lout = sum(lelm([elmd;elme])); % Total length of the moving section
    ldmax = 0.05/1.05*(lfixed+lout); % At least 5% overlap
    
    ltotmax = lfixed+lout-ldmax;
    ltotmin = max([lfixed+ldmax,lout+ldmax]); % At minimum, the minimum length of the fixed and outboard sections is at least ldmin
    
    Nelf = length(elmf);
    Neld = length(elmd);
    Nele = length(elme);
    
    % Minimum configuration
    lfixedmin = ltotmin - lout;
    loutmin = ltotmin - lfixed;
    ldmin = ltotmin - lfixedmin - loutmin;
    
    morph.span.lelmmin(elmf) = lfixedmin/Nelf;
    morph.span.lelmmin(elmd) = ldmin/Neld;
    morph.span.lelmmin(elme) = loutmin/Nele;
    
    % Maximum configuration
    lfixedmax = ltotmax - lout;
    loutmax = ltotmax - lfixed;
    
    % Extension length per element
    morph.span.dl(elmf) = (lfixedmax-lfixedmin)/Nelf;
    morph.span.dl(elmd) = (ldmax-ldmin)/Neld;
    morph.span.dl(elme) = (loutmax-loutmin)/Nele;
    
    % Initial extension
    morph.span.ext0(elmf) = (lelm(elmf)-morph.span.lelmmin(elmf))./morph.span.dl(elmf);
    morph.span.ext0(elmd) = (lelm(elmd)-morph.span.lelmmin(elmd))./morph.span.dl(elmd);
    morph.span.ext0(elme) = (lelm(elme)-morph.span.lelmmin(elme))./morph.span.dl(elme);
    
    morph.span.sec([elmf;elmd;elme]) = 1;

%     morph.span.ext([elmf;elmd;elme]) = morph.inp.span.ext;
    morph.span.dextddv([elmf;elmd;elme]) = 1;
    
%     morph.span.low([elmf;elmd;elme]) = morph.inp.span.low;
%     morph.span.high([elmf;elmd;elme]) = morph.inp.span.high;
end

% Fold
if sum(morph.inp.fold.sec)~=0
    morph.fold.sec = zeros(str.Ns,1);
%     morph.fold.angle = zeros(str.Ns,1);
    morph.fold.dangleddv = zeros(str.Ns,sum(morph.inp.fold.sec));
%     morph.fold.low = zeros(str.Ns,1);
%     morph.fold.high = zeros(str.Ns,1);
    foldlocs = find(morph.inp.fold.sec==1);
    
    for i=1:length(morph.inp.fold.sec)
        if morph.inp.fold.sec(i) == 1
            elm = find(str.sec==i);
            morph.fold.sec(elm) = 1;
%             morph.fold.angle(elm) = deg2rad(morph.inp.fold.angle(i));
            morph.fold.dangleddv(elm,foldlocs==i) = 1;
%             morph.fold.low(elm) = deg2rad(morph.inp.fold.low(i));
%             morph.fold.high(elm) = deg2rad(morph.inp.fold.high(i));
        end
    end
end
