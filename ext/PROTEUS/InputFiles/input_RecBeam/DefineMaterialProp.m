function [constant] = DefineMaterialProp(constant)

% =====                                                              ==== %
%                        MaterialProp [17/07/2015]                        %
%                                                                         %
%  Material properties are defined and assigned in this .m file
%
%  The material properties are stored in
%       constant.mat
% =====                                                              ==== %

% ---
if 1    % Aluminium
    cprintf('green', 'Aluminium \n');
    % 7075T73
    constant.mat.ID(1) = 1;
    constant.mat.rho(1) = 2810;
    constant.mat.E11(1) = 70e9;
    constant.mat.E22(1) = 70e9;
    constant.mat.G12(1) = 26.6e9;
    constant.mat.nu12(1) = 0.3;
%     constant.mat.exmax(1) = 3146e-6;
%     constant.mat.gmax(1) = 6738e-6;
    % Convervative allowables
    constant.mat.Xt(1) = 2280e6*0.8*0.65*0.8;
    constant.mat.Xc(1) = 1725e6*0.8*0.65*0.8;
    constant.mat.Yt(1)  = 57e6*0.8*0.65*0.8;
    constant.mat.Yc(1)  = 228e6*0.8*0.65*0.8;
    constant.mat.S(1)  = 76e6*0.8*0.65*0.8;
    
    constant.mat.exmax(1) = 2280e6/constant.mat.E11(1)*0.8*0.65*0.8;
    constant.mat.exmin(1) = -1725e6/constant.mat.E11(1)*0.8*0.65*0.8;
    constant.mat.gmax(1)  = 76e6/constant.mat.G12(1)*0.8*0.65*0.8;
    
    constant.mat.ID(2)    = 2;
    constant.mat.E11(2)   = 68.95e9;
    constant.mat.E22(2)   = 68.95e9;
    constant.mat.G12(2)   = 68.95e9/(2*(1+0.3));
    constant.mat.nu12(2)  = 0.3;
end
% ---

% ---
if 0    % Composite T300/5208
    cprintf('green', 'Composite \n');
    constant.mat.ID(1)    = 1;
    constant.mat.rho(1)   = 1452;
    constant.mat.E11(1)   = 8.30e10;
    constant.mat.E22(1)   = 8.50e9;
    constant.mat.G12(1)   = 4.20e9;
    constant.mat.nu12(1)  = 0.35;
    % Convervative allowables
    constant.mat.exmax(1) = 4500e-6;
    constant.mat.exmin(1) = -4500e-6;
    constant.mat.gmax(1)  = 7000e-6;
    
    % Stringer material properties
    constant.mat.ID(2)    = 1;
    constant.mat.E11(2)   = 68.95e9;
    constant.mat.E22(2)   = 68.95e9;
    constant.mat.G12(2)   = 68.95e9/(2*(1+0.3));
    constant.mat.nu12(2)  = 0.3;
end
% ---

if 0    % Composite AS4/3501-6
    cprintf('green', 'Composite \n');
    constant.mat.ID(1)    = 1;
    constant.mat.rho(1)   = 1600;
    constant.mat.E11(1)   = 147e9;
    constant.mat.E22(1)   = 10.3e9;
    constant.mat.G12(1)   = 7.0e9;
    constant.mat.nu12(1)  = 0.27;
    % Convervative allowables
    constant.mat.Xt(1) = 2280e6*0.8*0.65*0.8;
    constant.mat.Xc(1) = 1725e6*0.8*0.65*0.8;
    constant.mat.Yt(1)  = 57e6*0.8*0.65*0.8;
    constant.mat.Yc(1)  = 228e6*0.8*0.65*0.8;
    constant.mat.S(1)  = 76e6*0.8*0.65*0.8;
    
    constant.mat.exmax(1) = 2280e6/constant.mat.E11(1)*0.8*0.65*0.8;
    constant.mat.exmin(1) = -1725e6/constant.mat.E11(1)*0.8*0.65*0.8;
    constant.mat.gmax(1)  = 76e6/constant.mat.G12(1)*0.8*0.65*0.8;
    
    
    % Stringer material properties
    constant.mat.ID(2)    = 1;
    constant.mat.E11(2)   = 68.95e9;
    constant.mat.E22(2)   = 68.95e9;
    constant.mat.G12(2)   = 68.95e9/(2*(1+0.3));
    constant.mat.nu12(2)  = 0.3;
end