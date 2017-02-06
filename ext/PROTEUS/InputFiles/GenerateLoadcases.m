function loadcase=GenerateLoadcases(ExcelFile)

% global xlsFileName
%% Define load cases
% Fuel Cases
% Tank ID | Pos. ID     | Location [m]   | Mass [%]
% --------------------------------------------
%    1    |   root      |    0 - 5       |  33%
%    2    |  main wing  |    5 - 21      |  62%
%    3    |    tip      |   21 - 29.3    |   5%


LoadcaseExcel = xlsread(ExcelFile,'Loadcases');
loadflag      = 1;%size(LoadcaseExcel,1);
for iLC = 1:loadflag
    M  (iLC,1)  = LoadcaseExcel(iLC,2);
    EAS(iLC,1)  = LoadcaseExcel(iLC,3);
    H  (iLC,1)  = LoadcaseExcel(iLC,4);
    LCload(iLC,1) = LoadcaseExcel(iLC,5);
    nz (iLC,1)  = LoadcaseExcel(iLC,6);
    gustflag(iLC,1) = LoadcaseExcel(iLC,7);
    trim(iLC,1) = LoadcaseExcel(iLC,8);
    alpha0(iLC,1) = LoadcaseExcel(iLC,9);
    fuel_level{iLC,1} = LoadcaseExcel(iLC,10:end);
end



% Load Case Definition
% ---
if 0    % 1 load case
    LoadTable = [ 
%       'Mach(at true airspeed)'   'EAS (m/s)'  'Altitude (m)'  'nz'
%           0.85	83.389  0	     1;     % cruise (RE = 40 million = 1.225*83.389*7.00532/1.789e-5)
%           0.85	140.31  10668	 1;     % cruise (RE = 67 million = 1.225*140.31*7.00532/1.789e-5) at sea level or 46.418 million 0.3796*252.06*7.00532/1.444e-5 at 35000ft
%           0.85	133.73  11277.6	 1;     % cruise (RE = 64 million = 1.225*133.73*7.00532/1.789e-5) at sea level or 42.556 million 0.3483*249.76*7.00532/1.432e-5 at 37000ft
%           0.85  130.51  11582.4  1;     % cruise (RE = 62 million = 1.225*130.51*7.00532/1.789e-5) at sea level or 38.556 million 0.3164*249.76*7.00532/1.432e-5 at 39000ft
          0.85  127.34  11887.2  1;     % cruise (RE = 61 million = 1.225*127.34*7.00532/1.789e-5) at sea level or 40.564 million 0.3320*249.76*7.00532/1.432e-5 at 39000ft
%         0.75	213.5	2913	2.5;   % LC 1    
%         0.75	185.4	5078	1;   
%         0.75	148.5	8252	-1;        
%         0.5   213.5	-4000	2.5;
        ];
  
  LoadTable  = [LoadTable];
  
  for loadcaseID = 1 : size(LoadTable,1)
      M  (loadcaseID,1)  = LoadTable(loadcaseID,1);
      EAS(loadcaseID,1)  = LoadTable(loadcaseID,2);
      H  (loadcaseID,1)  = LoadTable(loadcaseID,3);
      nz (loadcaseID,1)  = LoadTable(loadcaseID,4);
      
      fuel_level{loadcaseID} = 0.1*[1 1 1]
  end
   loadflag   = size(M,1);
end
% ---



% ---
if 0 % static tip load deflection
    loadcaseID      = 1;
    nz(loadcaseID)  = 1;
 
    M(loadcaseID)   = 0;                               % mach 
    EAS(loadcaseID) = 1;                               % EAS
    H(loadcaseID)   = 0;                                    % altitude (m)

    % Fuel settings
    fuel_level{loadcaseID} = 0*[1  1  1];
    loadflag = 1;
end
% ---



% ---
if 0 % Find Critial Load Cases
    
    loadcase.findcritical = true;
    NumberOfLoadcase      = 25;              % input

    % ---
    if 1 % LHS space filling design for loadcases + Extreme ones
        QpsiLimits   = [0.1 7.5];
        AltLimits    = [-4000 12000];
        NZLimits     = [-1 2.5];
        FuelLimits   = [0 1];
        
        SpaceFillingArray = lhsdesign(NumberOfLoadcase*2,4);
        
        H        = SpaceFillingArray(:,1)* sum(abs(AltLimits))   + AltLimits(1);
        nz       = SpaceFillingArray(:,2)* sum(abs(NZLimits))    + NZLimits(1);
        Qpsi     = SpaceFillingArray(:,3)* sum(abs(QpsiLimits))  + QpsiLimits(1);
        Q        = Qpsi * 6894;
        EAS      = sqrt(2*Q/1.225);
        Fuelcoef = SpaceFillingArray(:,4)* sum(abs(FuelLimits)) + FuelLimits(1);
        
        for loadcaseID = 1 : NumberOfLoadcase*2
            [~, M(loadcaseID,1)]   = eas2tas(EAS(loadcaseID) ,H(loadcaseID));
            fuel_level{loadcaseID} = Fuelcoef(loadcaseID)*rand(1,3);
%             fuel_level{loadcaseID} = Fuelcoef(loadcaseID)*[1 1 1];
        end
        
        loadcaseID = NumberOfLoadcase;
        % remove loadcases where Mach # > 0.75
        try
            H          = H(M<0.75);           H          = H(1:NumberOfLoadcase,1);
            nz         = nz(M<0.75);          nz         = nz(1:NumberOfLoadcase,1);
            EAS        = EAS(M<0.75);         EAS        = EAS(1:NumberOfLoadcase,1);
            fuel_level = fuel_level(M<0.75);  fuel_level = fuel_level(1:NumberOfLoadcase);
            M          = M(M<0.75);           M          = M(1:NumberOfLoadcase,1);
        catch
            error('try restarting')
        end
            
        % --- Add extreme cases
        tempo = combvec(QpsiLimits,AltLimits,NZLimits,FuelLimits);
        for i = 1 : size(tempo,2)
            loadcaseID       = loadcaseID + 1;
            H(loadcaseID,1)  = tempo(2,i);
            nz(loadcaseID,1) = tempo(3,i);
            fuel_level{loadcaseID} = tempo(4,i)*[1 1 1];
            
            Qpsi                 = tempo(1,i);
            EAS(loadcaseID,1)    = sqrt(2*Qpsi*6894/1.225);
            [~, M(loadcaseID,1)] = eas2tas(EAS(loadcaseID) ,H(loadcaseID));
        end
        
        H          = H(M<0.75);          
        nz         = nz(M<0.75);         
        EAS        = EAS(M<0.75);         
        fuel_level = fuel_level(M<0.75);    
        M          = M(M<0.75);           
    end
    % ---
    
    % ---
    if 0 % Enumaration
        for Qpsi = 0.1 : 1: 8.1
            Q = Qpsi * 6894; % dynamic pressure in Pascal
            
            %         for iNz = -1 : 1 : 2
            
            for iAlt = -4000 : 16000 : 12000          % Altitude loop (in meter)
                
                Veas        = sqrt(2*Q/1.225) + 0.0001;
                [~, Mach]   = eas2tas(Veas ,iAlt);
                
                if Mach < 0.8
                    
                    iNz = NZArray(randi([1 length(NZArray)],1));
                    
                    loadcaseID = loadcaseID +1;
                    fuel_level{loadcaseID}{1} = FuelArray(randi([1,length(FuelArray)],1))*[1 1 1];
                    numfuelcases(loadcaseID)  = 1;
                    nz(loadcaseID)            = iNz;
                    EAS(loadcaseID)           = Veas;                   % EAS
                    H(loadcaseID)             = iAlt;
                    M(loadcaseID)             = Mach;
                    
                    fuel_config               = 0 : 0.1 : 1;
                    numfuelcases(loadcaseID)  = length(fuel_config);
                    for iFuel = 1: length(fuel_config)
                        fuel_level{loadcaseID}{iFuel} = fuel_config(iFuel)*[1   1    1];                                     % fuel tanks level (in %)
                    end
                end
            end
            %         end
        end
    end
    % ---
    loadcaseID = size(M,1);
    loadflag   = loadcaseID;
else
    loadcase.findcritical = false;
end
% ---


%% Select Loadcase
if loadflag<=length(nz)
    loadcase.EAS        = EAS(1:loadflag);
    loadcase.H          = H(1:loadflag);
    loadcase.M          = M(1:loadflag);
    loadcase.nz         = nz(1:loadflag);
    loadcase.LCload     = LCload(1:loadflag);
    loadcase.gustflag   = gustflag(1:loadflag);
    loadcase.trim       = trim(1:loadflag);
    loadcase.alpha0     = alpha0(1:loadflag);
    loadcase.fuel_level = fuel_level(1:loadflag);
else
   err('Reference to non-existent loadcase. Change loadflag in main.m (1 <= loadflag <= 6)') 
end
