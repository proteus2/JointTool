classdef aeroElasticModel
    %AEROSTUCTURALMODEL Summary of this class goes here
    %   Detailed explanation goes here
    properties (Abstract)
    end
    properties 
        aeroModel
        strModel
        settings
    end
    
    methods (Abstract)
        obj=transformDisplacements(obj);
        obj=transformForces(obj);
    end
      
    methods
        function obj=aeroElasticModel(inAeroModel,inStrModel)
            obj.aeroModel=inAeroModel;
            obj.strModel=inStrModel;
            obj.settings=aeroElasticModelSettings();
        end
        function obj=solve(obj)
            err=1;
            while err>obj.settings.relConvTol
                %store current displacements;
                prvDisp=obj.strModel.disp;
                %calculate forces
                obj.aeroModel=obj.aeroModel.solve();
                %transform forces and set them to struct model
                obj=obj.transformForces();
                %calculate displacements
                obj.strModel=obj.strModel.solve();
                %transform grid and set updated grid
                obj=obj.transformDisplacements();
                %evaluate convergence (tbd)
                Disp=obj.strModel.disp;
                err=max((Disp./prvDisp)-1);
             end
        end
        function obj=plotResults(obj)
            figure; 
            obj.aeroModel.plotGrid();
%             obj.strModel.plotGrid();
            obj.aeroModel.plotGridDeflected();     
            obj.strModel.plotGridDeflected();       
        end
    end
    
end

