classdef strModelSettings
    % STRUCTURALMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lin  % Linear
        grav % Gravity
        NLCC % Non-Linear Convergence Criterion (Not used a.t.m --> Implement in seperate property class - protected)
        Fext % External forces
    end
    
    methods
       function obj=strModelSettings() % constructor
           obj.lin  = 0;
           obj.grav = 0;
           obj.NLCC = 0.0001;
           obj.Fext = 0;
       end
    end 
end

