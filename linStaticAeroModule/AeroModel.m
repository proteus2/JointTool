classdef AeroModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=Hidden)
        geometry
    end
    
    methods
        function obj=AeroModel(xmlFile)
            addpath('dAEDalusNXT/aerodynamics')
            addpath('dAEDalusNXT/geometry')
            obj.geometry=class_aircraft(xmlFile,1);
        end
        function obj=calculateMatrices()
            
        end
        function obj=calculateForces()
            
        end
        function obj=updateGrid(grid)
            
        end
    end
    
end

