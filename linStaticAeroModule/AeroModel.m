classdef AeroModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess=private)
        %Geometry class
        geometry
        %Vortex Lattice Kernel
        VLM
    end
    properties 
        %Body grid for coupling
        grid
        %Forces in Body Fixed Frame
        F 
    end
    
    methods
        function obj=AeroModel(xmlFile, aeroState) %constructor
            %read geofile
            obj.geometry=class_aircraft(xmlFile,1);
            %aero grid settings
            obj.geometry.grid_settings.x_max_grid_size= obj.geometry.reference.c_ref/16;
            obj.geometry.grid_settings.y_max_grid_size= obj.geometry.reference.b_ref/80;
            %compute grid
            obj.geometry=obj.geometry.compute_grid();
            %initialize vlm kernel
            obj.VLM=class_VLM_solver(obj.geometry.grid,obj.geometry.te_idx, obj.geometry.panels,aeroState, obj.geometry.reference);
            % use nvec2 for normal vector correction
            obj.VLM=obj.VLM.update_nvec2(obj.geometry);
            %store grid for coupling algorithms
            obj.grid=obj.VLM.grid;
        end
        function obj=calculateForces(obj)
            %this function calculates the forces 
            obj.VLM=obj.VLM.f_solve_std;
            obj.F=obj.VLM.F_body;
        end
        function plotGrid(obj)
            %plots the computed grid
            obj.geometry.plot_grid();
        end
        function plotCp(obj)
            %plots the computed grid
            obj.VLM.plot_cp();
        end
        function obj=updateGrid(obj)
            
        end
    end
    
end

