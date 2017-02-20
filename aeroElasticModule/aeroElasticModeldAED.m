classdef aeroElasticModeldAED < aeroElasticModel
    %LINSTATICAEROSTRUCTURALMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(Hidden)
        dummyDaedStr
    end
    
    methods
        function obj=aeroElasticModeldAED(inAeroModel, inStrModel) %constructor
            obj=obj@aeroElasticModel(inAeroModel,inStrModel); %call constructor of superclass
            %dAEDalus coupling initialization
            obj.aeroModel.geometry.wings(1)=obj.aeroModel.geometry.wings(1).compute_wingbox_coords(2);
            obj.dummyDaedStr=class_beam_collection(class_wing(size(obj.strModel.grid,1)-1,'none'));
            %fill geometry with required information
            [obj.dummyDaedStr.beam(1).beamelement(:).T]=deal(eye(12,12));
            obj.dummyDaedStr.beam(1).isExternalFEM=1;
            obj.dummyDaedStr.beam(1).Kff=zeros(size(obj.strModel.grid,1)*6);
            obj.dummyDaedStr.beam(1).Mff=zeros(size(obj.strModel.grid,1)*6);
            obj.dummyDaedStr.beam(1).update_K=0;
            obj.dummyDaedStr.beam(1).update_M=0;
            obj.dummyDaedStr.beam(1).node_coords=obj.strModel.grid;
            obj.dummyDaedStr.settings.gravity=0;
            obj.dummyDaedStr.settings.fuel_mass=0;
            obj.dummyDaedStr.settings.engine=0;
            obj.dummyDaedStr.settings.nonlinear=0;
            obj.dummyDaedStr.settings.landing_gear=0;
            
            obj.dummyDaedStr.beam(1).epsilon=zeros(size(obj.strModel.grid,1),1);
            %Rectangular
            obj.dummyDaedStr.beam(1).dist_c4_sc=2.5*ones(size(obj.strModel.grid,1),1);
            obj.aeroModel.geometry.wings(1).wingbox_coords=repmat(obj.strModel.grid(31:end,:)',1,1,2);
            obj.aeroModel.geometry.wings(1).wingbox_coords=repmat(obj.strModel.grid(31:end,:)',1,1,2);
            obj.aeroModel.geometry.wings(1).wing_segments(1).wingbox_coords=repmat(obj.strModel.grid(31:end,:)',1,1,2);
             obj.aeroModel.geometry.wings(1).wingbox_rl_coords=[0 cumsum(sqrt(diff(obj.strModel.grid(31:end,1)).^2+diff(obj.strModel.grid(31:end,2)).^2+diff(obj.strModel.grid(31:end,3)).^2))'];
            %CRM 
%             obj.aeroModel.geometry.wings(1).wingbox_coords=repmat(obj.strModel.grid(21:end,:)',1,1,2);
%             obj.aeroModel.geometry.wings(1).wingbox_coords=repmat(obj.strModel.grid(21:end,:)',1,1,2);
%             obj.aeroModel.geometry.wings(1).wing_segments(1).wingbox_coords=repmat(obj.strModel.grid(21:24,:)',1,1,2);
%             obj.aeroModel.geometry.wings(1).wing_segments(2).wingbox_coords=repmat(obj.strModel.grid(24:29,:)',1,1,2);
%             obj.aeroModel.geometry.wings(1).wing_segments(3).wingbox_coords=repmat(obj.strModel.grid(29:31,:)',1,1,2);
%             obj.aeroModel.geometry.wings(1).wing_segments(4).wingbox_coords=repmat(obj.strModel.grid(31:38,:)',1,1,2);
%             obj.aeroModel.geometry.wings(1).wing_segments(5).wingbox_coords=repmat(obj.strModel.grid(38:41,:)',1,1,2);
%             obj.aeroModel.geometry.wings(1).wingbox_rl_coords=[0 cumsum(sqrt(diff(obj.strModel.grid(21:end,1)).^2+diff(obj.strModel.grid(21:end,2)).^2+diff(obj.strModel.grid(21:end,3)).^2))'];
            for i=1:length(obj.dummyDaedStr.beam(1).beamelement(:))
                obj.dummyDaedStr.beam(1).beamelement(i).le=norm(obj.strModel.grid(i,:)'-obj.strModel.grid(1+i,:)');
            end
            %compute force interpolation information
            obj.aeroModel.geometry=obj.aeroModel.geometry.compute_force_interpolation_matrix(obj.dummyDaedStr);
            close
        end
        

        function obj=transformForces(obj)
            obj.aeroModel.geometry=obj.aeroModel.geometry.compute_beam_forces(obj.aeroModel.F, obj.aeroModel.F*0,obj.dummyDaedStr);
            obj.dummyDaedStr.beam(1)=obj.dummyDaedStr.beam(1).f_set_aeroloads(obj.aeroModel.geometry.wings(1));
            obj.dummyDaedStr=obj.dummyDaedStr.f_assemble(1,0);
            obj.strModel.Fs=obj.dummyDaedStr.Ftest;
        end
        function obj=transformDisplacements(obj)
            obj.dummyDaedStr.nodal_deflections=obj.strModel.disp;
            obj.dummyDaedStr=obj.dummyDaedStr.f_postprocess();
            obj.aeroModel.geometry=obj.aeroModel.geometry.compute_deflected_grid(obj.dummyDaedStr.f_get_deflections);
            obj.aeroModel.gridDeflected=obj.aeroModel.geometry.grid_deflected;
        end
    end
    
end

