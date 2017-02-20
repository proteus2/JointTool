classdef aeroElasticModelSettings
    %AEROSTUCTURALMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        relConvTol 
    end
    
    methods
       function obj=aeroElasticModelSettings() % constructor 
           %set default values
           obj.relConvTol=10e-3;  
       end
    end
    
end

