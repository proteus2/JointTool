classdef aeroStrCplSettings
    %AEROSTUCTURALMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        linear
        static
        matrix
        
    end
    
    methods
       function obj=aeroStrCplSettings() % constructor
           obj.linear=1;
           obj.static=1;
           obj.matrix=1;
       end
    end
    
end

