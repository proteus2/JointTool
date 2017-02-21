function [constant] = GenerateLaminateIniGuess(constant,Nlam,FLAGSPAR)

% =====                                                              ==== %
%                       LaminatesIniGuess [14/04/2015]                    %
%                                                                         %
%   generates an initial stacking sequence guess for all laminates
%   the initial guess is quasi-isotropic
%
%   All initial design data is stored in
%       constant.ini
% =====                                                              ==== %

% For sake of clarity, laminates fiber angles are defined w.r.t. the global y-axis, positive aft and down. 
% Thickness Percentages 
perc_pm45 = 50;
perc_90   = 25;
perc_0    = 25;

% Initial spanwise thickness distribution based on number of spanwise laminates
% thicknessDistGuess = [linspace(0.01,0.02,sum(Nlam.Span(1:2))),linspace(0.02,0.005,Nlam.Span(3))];
thicknessDistGuess = 0.005*ones(sum(Nlam.Span),1);


%% ----- Automated processing to distribute the laminates
% Initial Guess Design
constant.lam.ini.ID    = constant.lam.ID;
constant.lam.ini.matID = constant.lam.matID;

sweep = atand(diff(Nlam.xyz_lam(:,1))./diff(Nlam.xyz_lam(:,2)));
dihedral = atand(diff(Nlam.xyz_lam(:,3))./diff(Nlam.xyz_lam(:,2)));

for i=1:length(constant.lam.TopID)
    constant.lam.ini.t(constant.lam.TopID{i},1)   = thicknessDistGuess(i);
    constant.lam.ini.t(constant.lam.BotID{i},1)   = thicknessDistGuess(i);
    if FLAGSPAR
        constant.lam.ini.t(constant.lam.SparID{i},1)   = thicknessDistGuess(i);
    end
    constant.lam.OrientationShift(constant.lam.TopID{i},1) = -sweep(i);
    constant.lam.OrientationShift(constant.lam.BotID{i},1) = sweep(i);
    if FLAGSPAR
        constant.lam.OrientationShift(constant.lam.SparID{i},1) = dihedral(i);
    end
    
    % Top skin layup
    for ilam=1:length(constant.lam.TopID{i})
        Layup   = [-45*ones(ceil(perc_pm45/2),1); 45*ones(ceil(perc_pm45/2),1); 90*ones(perc_90,1); 0*ones(perc_0,1)];
        constant.lam.ini.layup{constant.lam.TopID{i}(ilam)} = Layup+ constant.lam.OrientationShift(constant.lam.TopID{i}(ilam));
    end
    
    % Bottom skin layup
    for ilam=1:length(constant.lam.BotID{i})
        % Note minus sign because bottom skin is trailing edge to leading
        % edge positive
        Layup   = -[-45*ones(ceil(perc_pm45/2),1); 45*ones(ceil(perc_pm45/2),1); 90*ones(perc_90,1); 0*ones(perc_0,1)];
        constant.lam.ini.layup{constant.lam.BotID{i}(ilam)} = Layup+ constant.lam.OrientationShift(constant.lam.BotID{i}(ilam));
    end
    
    % Spar skin layup
    if FLAGSPAR
        for ilam=1:length(constant.lam.SparID{i})
            Layup   = [45*ones(ceil(perc_pm45/2),1); -45*ones(ceil(perc_pm45/2),1); 90*ones(perc_90,1); 0*ones(perc_0,1)];
            constant.lam.ini.layup{constant.lam.SparID{i}(ilam)} = Layup+ constant.lam.OrientationShift(constant.lam.SparID{i}(ilam));
        end
    end
end
