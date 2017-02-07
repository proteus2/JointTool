classdef aeroModel < handle
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
        % colloc points
        colloc
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
            %store collocation point coordinates
            obj.colloc=obj.VLM.colloc;
        end
        function obj=calculateForces(obj)
            %this function calculates the forces 
            obj.VLM=obj.VLM.f_solve_std;
            obj.F=obj.VLM.F_body;
        end
        function obj=updateGrid(obj,gridIn)
            obj.geometry.grid_deflected=gridIn;
            obj.VLM=class_VLM_solver(obj.geometry.grid_deflected,obj.geometry.te_idx, obj.geometry.panels,obj.VLM.state, obj.geometry.reference);
            % use nvec2 for normal vector correction
            obj.VLM=obj.VLM.update_nvec2(obj.geometry);
            %store grid for coupling algorithms
            obj.grid=obj.VLM.grid;
            %store collocation point coordinates
            obj.colloc=obj.VLM.colloc;
        end
    end
    
    
    
    
    methods %static
        function plotGrid(obj)
            %plots the computed grid
            obj.geometry.plot_grid();
        end
        function plotGeo(obj)
            %plots the geometry
            obj.geometry.plot();
        end
        function plotGridVol(obj)
            %plots the computed volume Grid
            obj.geometry.plot_grid_vol();
        end
        function plotCp(obj)
            %plots the computed grid
            obj.VLM.plot_cp();
        end
    end
    
end

