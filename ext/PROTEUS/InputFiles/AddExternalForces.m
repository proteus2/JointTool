function [fext] = AddExternalForces(fext,xyzshift,constant,wingID)

global xlsFileName

% =====                                                              ==== %
%                        External Forces [14/04/2015]                     %
%                                                                         %
%    here are defined the external forces acting on the wing
%    (e.g. weight of the engine, pylon, etc ... ) 
% =====                                                              ==== %

Wing_LE_xyz      = constant.Coord3D.Wing_LE_xyz;
Wing_LE_xyz(:,1) = Wing_LE_xyz(:,1)+xyzshift(1);
Wing_LE_xyz(:,2) = Wing_LE_xyz(:,2)+xyzshift(2);
Wing_LE_xyz(:,3) = Wing_LE_xyz(:,3)+xyzshift(3);
Wing_TE_xyz      = constant.Coord3D.Wing_TE_xyz;
Wing_TE_xyz(:,1) = Wing_TE_xyz(:,1)+xyzshift(1);
Wing_TE_xyz(:,2) = Wing_TE_xyz(:,2)+xyzshift(2);
Wing_TE_xyz(:,3) = Wing_TE_xyz(:,3)+xyzshift(3);

[Num,txt] = xlsread(xlsFileName,'ExternalForces');
txt       = txt(2:end,1);

for i = 1:length(txt)
    EntryNum                  = i;
    fext.type{EntryNum}       = txt{i};
    fext.magnitude{EntryNum}  = Num(i,1:6); 
    Xloc = Num(i,7);
    if strcmp(wingID,'right')
        Yloc = Num(i,8);
    elseif strcmp(wingID,'left')
        Yloc = -Num(i,8);
    end
    Zloc = Num(i,9);
    if ~isnan(Zloc)
        fext.location{EntryNum} = [Xloc Yloc Zloc];
    else
        XZ_LElocal = interp1(Wing_LE_xyz(:,2),Wing_LE_xyz(:,[1 3]),Yloc);
        XZ_TElocal = interp1(Wing_TE_xyz(:,2),Wing_TE_xyz(:,[1 3]),Yloc);
        Zloc = interp1([XZ_LElocal(1) XZ_TElocal(1)] , [XZ_LElocal(2) XZ_TElocal(2)], Xloc);
        fext.location{EntryNum} = [Xloc Yloc Zloc];
    end
    fext.follower{EntryNum}  = Num(i,10);
    fext.alphaflag{EntryNum} = Num(i,11);
end
end
% ---
