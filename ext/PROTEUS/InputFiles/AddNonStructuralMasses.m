function [lumped] = AddNonStructuralMasses(xyzshift,constant,lumped,wingID)

global xlsFileName

% =====                                                              ==== %
%                   Non-structural masses [14/04/2015]                    %
%                                                                         %
%  Define the location and mass of non-structural masses (Kg)
% =====                                                              ==== %

Wing_LE_xyz   = constant.Coord3D.Wing_LE_xyz;
Wing_LE_xyz(:,1) = Wing_LE_xyz(:,1)+xyzshift(1);
Wing_LE_xyz(:,2) = Wing_LE_xyz(:,2)+xyzshift(2);
Wing_LE_xyz(:,3) = Wing_LE_xyz(:,3)+xyzshift(3);
Wing_TE_xyz   = constant.Coord3D.Wing_TE_xyz;
Wing_TE_xyz(:,1) = Wing_TE_xyz(:,1)+xyzshift(1);
Wing_TE_xyz(:,2) = Wing_TE_xyz(:,2)+xyzshift(2);
Wing_TE_xyz(:,3) = Wing_TE_xyz(:,3)+xyzshift(3);

[Num,txt] = xlsread(xlsFileName,'NonStructuralMasses');
txt       = txt(2:end,1);

for i = 1:length(txt)
    
    if ~isfield(lumped,'type') || isempty(find(strcmp(txt{i},lumped.type)))
        if isfield(lumped,'type')
            EntryNum                  = length(lumped.type)+1;
        else
            EntryNum                  = 1;
        end
        lumped.type{EntryNum}     = txt{i};
        massflag                  = 1;
    else
        EntryNum = find(strcmp(txt{i},lumped.type));
        massflag = 0;
    end
    
    if massflag == 1
        lumped.mass{EntryNum}     = Num(i,1);
    else
        try
            lumped.mass{EntryNum}     = [lumped.mass{EntryNum};Num(i,1)];
        catch
            keyboard
        end
    end
    
    Xloc = Num(i,2);
    if strcmp(wingID,'right')
        Yloc = Num(i,3);
    elseif strcmp(wingID,'left')
        Yloc = -Num(i,3); 
    end
    Zloc = Num(i,4);
    
    if isnan(Zloc)
        XZ_LElocal = interp1(Wing_LE_xyz(:,2),Wing_LE_xyz(:,[1 3]),Yloc);
        XZ_TElocal = interp1(Wing_TE_xyz(:,2),Wing_TE_xyz(:,[1 3]),Yloc);
        Zloc = interp1([XZ_LElocal(1) XZ_TElocal(1)] , [XZ_LElocal(2) XZ_TElocal(2)], Xloc);
    end
    if massflag == 1
        lumped.location{EntryNum} = [Xloc Yloc Zloc];
    else
        lumped.location{EntryNum} = [lumped.location{EntryNum};Xloc Yloc Zloc];
    end
    
    IG                        = [Num(i,5),Num(i,8),Num(i,9)
                                 Num(i,8),Num(i,6),Num(i,10)
                                 Num(i,9),Num(i,10),Num(i,7)];  
    if massflag == 1
        lumped.IG{EntryNum}       = IG;                     % specify the mass inertia matrix % A 3x3 matrix per mass in input
    else
        lumped.IG{EntryNum}       = [lumped.IG{EntryNum},IG];                     % specify the mass inertia matrix % A 3x3 matrix per mass in input
    end
end
end
%---
