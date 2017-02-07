classdef aeroStrModel
    %AEROSTUCTURALMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        AeroModel
        StrModel
        Tas
        Tsa
        cplSettings = aeroStrCplSettings();
    end
    
    methods
       function obj=aeroStrModel(inAeroModel, inStrModel) %constructor
           obj.AeroModel=inAeroModel;
           obj.StrModel=inStrModel;
       end
       function obj=computeCouplingMatrices(obj)
           %depending on settings calculate the matrices
           %
           % TODO
           %
           obj.Tas='';
           obj.Tsa='';
       end
    end
    
end

