classdef aeroModel 
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(Hidden) %hidden
        %Geometry class
        geometry
    end
    properties %(Access = private) % private
        %Vortex Lattice Kernel
        VLM
        %Unsteady Vortex Lattice Kernel
        UVLM
    end
    properties %public
        %Body grid for coupling
        gridDeflected
        % aerodynamic state
        state
    end
    properties %(SetAccess = protected)%public protected
        %Undeformed body grid 
        grid
        %Forces in Body Fixed Frame
        F 
        % force vector application points on undeformed body grid
        fvap
    end
    
    methods % constructor
        function obj=aeroModel(xmlFile, aeroState) %constructor
            %read geofile
            obj.geometry=class_aircraft(xmlFile,1);
            %aero grid settings
            obj.geometry.grid_settings.x_max_grid_size= obj.geometry.reference.c_ref/16;
            obj.geometry.grid_settings.y_max_grid_size= obj.geometry.reference.b_ref/40;
            %compute grid
            obj.geometry=obj.geometry.compute_grid();
            %initialize vlm kernel
            obj.VLM=class_VLM_solver(obj.geometry.grid,obj.geometry.te_idx, obj.geometry.panels,aeroState, obj.geometry.reference);
            %store force vector application points
            obj.fvap=obj.VLM.r;
            %store state
            obj.state=aeroState;
        end
    end
    methods % set
        function obj=set.gridDeflected(obj,inGrid)
            disp('grid changed')
            obj.gridDeflected=inGrid;
            obj=obj.updateGrid(); %update kernel on set
        end
        function obj=set.state(obj,inState)
            obj.state=inState;
             obj=obj.updateState(); 
        end
    end
    methods %get
        function grid=get.grid(obj)
            disp('grid output:')
            grid=obj.geometry.grid;
        end
    end
    methods
        function obj=solve(obj)
            %initialize VLM kernel with gridDeflected
            obj.VLM=class_VLM_solver(obj.geometry.grid_deflected,obj.geometry.te_idx, obj.geometry.panels,obj.state, obj.geometry.reference);
            %this function calculates the forces 
            obj.VLM=obj.VLM.f_solve_std;
            obj.F=obj.VLM.F_body;
        end
        function obj=updateGrid(obj)
            obj.geometry.grid_deflected=obj.gridDeflected;
            
        end
        function obj=updateState(obj)
            %update VLM kernel with obj.state
            obj.VLM.state=obj.state;
        end
        function obj=initializeDynSimulation(obj, state,dt)
            obj.geometry.grid_settings.wake=1;
            obj.geometry=obj.geometry.compute_grid();            
            UVLM_settings=class_UVLM_computation_settings();
            UVLM_settings.debug=0;
            UVLM_settings.movie=1;
            UVLM_settings.wakelength_factor=1;
            
            obj.UVLM=class_UVLM_solver(obj.geometry.name,obj.geometry.grid,obj.geometry.is_te,obj.geometry.panels,state,obj.geometry.grid_wake,obj.geometry.panels_wake,obj.geometry.reference,UVLM_settings);
            obj.UVLM=obj.UVLM.initialize_time_domain_solution(dt);
            obj.F=obj.UVLM.F_body;
        end
        function obj=solveTimeStep(obj,state,i)
             obj.UVLM=obj.UVLM.solve_time_domain_aerodynamics(obj.geometry,[0*[0 0 0]' ([-1 0 0;0 1 0; 0 0 -1]*[0 0 0]')],state.V_inf,state.alpha,0,[0 0 0],state.rho_air,i,['results/test']);
             obj.F=obj.UVLM.F_body;
        end
    end
    methods %static
        function plotGrid(obj)
            %plots the computed grid
            obj.geometry.plot_grid();
        end
        function plotGridDeflected(obj)
            %plots the computed grid
            obj.geometry.plot_grid_deflected();
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

