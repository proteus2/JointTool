classdef structuralModel
    properties
        grid
        disp
        Fs
        Fext
        settings = strModelSettings();
    end
    properties (Access = private)
        Ks
        Kfext
        frdof
        fxdof
        dFextdsc
    end
    methods % Constructor
        function obj = structuralModel(model)
            [constant,~,~,~,statics] = genStrModel(model);
            % Extract properties
            obj.grid  = reshape(constant.str.xyz,3,constant.str.Ns+1)';
            obj.Ks    = statics.str.Ks;
            obj.Kfext = statics.str.Kfext;
            obj.Fs    = statics.str.Fs;
            obj.Fext  = statics.str.Fext;
            obj.disp  = statics.str.p;
            obj.frdof = constant.str.frdof;
            obj.fxdof = constant.str.fxdof;
            obj.fxdof = constant.str.fxdof;
            obj.dFextdsc = statics.sens.dFextdsc;
            % Save array
            curdir = cd;
            cd('ext/PROTEUS/results')
            save('constant.mat','constant')
            save('statics.mat','statics')
            cd(curdir)
        end
    end
    methods % Solution
        function obj = solve(obj)
            % Format (Nx, Ny, Nz, Mx, My, Mz)_node
            % Assemble stiffness and forces
            if ~obj.settings.Fext
                J = obj.Ks(obj.frdof,obj.frdof);
                R = obj.Fs(obj.frdof,1);
            else
                J = obj.Ks(obj.frdof,obj.frdof) - obj.Kfext;
                R = obj.Fs(obj.frdof,1) + obj.Fext(obj.frdof,1);
            end
            % Solve structural equation
            if obj.settings.lin
                %==========================================================%
                % LINEAR STRUCTURAL SOLUTION
                %==========================================================%
                % Calculate displacements
                obj.disp = [];
                obj.disp(obj.frdof,1) = J\R;
                % Update forces based on the structural deformations
                obj.Fs   = obj.Ks*obj.disp;
                obj.Fext = obj.Fext + obj.Kfext*obj.disp;
            else
                %==========================================================%
                % NON-LINEAR STRUCTURAL SOLUTION
                %==========================================================%
                % Init. displacements
                obj.disp = zeros(size(obj.disp));
                % Settings
                sc   = 0;
                Find = 0;
                Fi   = 0;
                Ff   = max(abs(R(3:6:end)));
                if Fi == Ff && Ff == 0
                    dsc = 0.25;
                else
                    if Fi == Ff
                        deltaF = min(1,0.1*Ff);
                        Fi = Ff-deltaF;
                    end
                    if Ff>Fi
                        step = ceil(4*(Ff-Fi)/Ff);
                    else
                        step = ceil(4*(Fi-Ff)/Fi);
                    end
                    dsc = 1/step;
                end
                count    = 0;
                exitflag = 0;
                tol1     = 1e-7; % Inner loop tolerance
                tol2     = 1e-4; % Outer loop tolerance
                fprintf('iter    Speed [m/s]    Tip T3 displacement [m]     check\n') 
                % Solve structural equation
                while (sc<1)
                    count = count + 1;
                    scold = sc;
                    Find  = Find + 1;
                    sc    = sc + dsc;
                    Fsq   = Fi^2 + sc*(Ff^2 - Fi^2);
                    F     = sqrt(Fsq);
                    if sc>1
                        F = Ff;
                        sc = 1;
                    end
                    check = 1;
                    iter = 0;
                    % Load structure
                    load('constant.mat')
                    load('statics.mat')
                    dold = statics.str.p;
                    % Predictor step
                    statics.str.p(obj.frdof,1) = statics.str.p(obj.frdof,1) +...
                        J\(statics.str.Fs(obj.frdof,1) + (F/Ff)*R + statics.sens.dFextdsc(obj.frdof,1)*(sc - scold));
                    % Start inner loop
                    while (check>tol1)
                        iter = iter + 1;
                        curdir = cd;
                        cd('ext/PROTEUS/statics/Kernel')
                        statics = structure(constant,statics,0,1,0,0);
                        statics = pgen_str(constant,statics,0,0);
                        statics = fext(constant,statics,sc,0,0,1,statics.str.alpha,0);
                        cd(curdir)
                        if ~obj.settings.Fext
                            R  = (F/Ff)*obj.Fs(obj.frdof,1) - statics.str.Fs(obj.frdof,1);
                            J  = statics.str.Ks(obj.frdof,obj.frdof);
                        else
                            R  = (F/Ff)*obj.Fs + statics.str.Fext(obj.frdof,1) - statics.str.Fs(obj.frdof,1);
                            J  = statics.str.Ks(obj.frdof,obj.frdof) - statics.str.Kfext(obj.frdof,obj.frdof);
                        end
                        dd = J\R;
                        statics.str.p(obj.frdof,1) = statics.str.p(obj.frdof,1) + dd;
                        check = norm(dd,inf)/norm(statics.str.p,inf);
                        if iter>20
                            obj.disp = dold;
                            sc = sc - dsc;
                            dsc = dsc/2;
                            fprintf('Reducing step size to %5.3e \n',dsc)
                            sc = sc + dsc;
                            Fsq = Fi^2 + sc*(Ff^2 - Fi^2);
                            F = sqrt(Fsq);
                            iter = 0;
                        end
                        if abs(dsc)<tol2
                            fprintf('Step size too small, no solution found \n')
                            check = 1E-9;
                            exitflag = 1;
                        end
                    end
                    % End inner loop
                    fprintf(' %2i       %3.1f              %5.3f              %4.2e',iter,F,statics.str.p(end-3),check)
                    fprintf('\n')
                    % Adapt step size beased on convergence
                    if iter<3 && exitflag~=1 && sc~=1
                        dsc = 2*dsc;
                        fprintf('Great convergence, increasing step size to %5.3e \n',dsc)
                    elseif iter>8
                        dsc = 0.5*dsc;
                        fprintf('Descreasing step size to %5.3e \n',dsc)
                    end
                    % Exit the loop if no solution if found
                    if exitflag == 1
                        break
                    end
                end
                obj.disp = statics.str.p;
            end
        end
    end
    methods % Graphics
        function obj = plotGridVol(obj,varargin)
            load('constant');
            if isempty(varargin)
                plotBeamElements(constant{1})
            else
                delta = varargin{1};
                figure()
                ax1 = axes;
                plotBeamElements(constant{1},0.3,0.03)
                colormap(ax1,'gray');
                freezeColors(ax1);
                plotBeamElements(constant{1},delta,1,0.1)
                colormap(ax1,'hot');
            end
        end
        function obj = plotGrid(obj,varargin)
            % Determine figure ID
            load('constant');
            if ~isempty(varargin)
                fID = varargin{1};
                figure(fID)
            end
            hold on
            % Extract properties
            P0 = obj.grid;
            C  = constant.inp.cbox;
            % Beam axis
            plot3(P0(:,1),P0(:,2),P0(:,3),...
                'ro-','LineWidth',2,'MarkerFaceColor',[1 0 0]);
            % Fish bones (visual aid for rotations)
            Null = zeros(constant.str.Ns+1,1);
            xref = constant.inp.xref;
            LE = P0 + [xref   Null Null].*[C Null Null];
            TE = P0 + [xref-1 Null Null].*[C Null Null];
            for i=1:constant.str.Ns+1
                plot3([LE(i,1), P0(i,1)],...
                    [LE(i,2), P0(i,2)],...
                    [LE(i,3), P0(i,3)],'g.-');
                plot3([TE(i,1), P0(i,1)],...
                    [TE(i,2), P0(i,2)],...
                    [TE(i,3), P0(i,3)],'g.-');
            end
            % Axes properties
            xlabel('Chord [m]')
            ylabel('Span [m]')
            zlabel('Height [m]')
            view(-50,25)
            axis equal
        end
        function obj = plotGridDef(obj,varargin)
            % Determine figure ID
            load('constant');
            if ~isempty(varargin)
                fID = varargin{1};
                figure(fID)
            end
            hold on
            % Extract properties
            P0 = obj.grid;
            P  = obj.disp;
            if size(P,2)==1
                P = reshape(P,6,constant.str.Ns+1)';
            end
            C  = constant.inp.cbox;
            % Deformed beam axis
            plot3(P0(:,1)+P(:,1),...
                P0(:,2)+P(:,2),...
                P0(:,3)+P(:,3),'ro-','LineWidth',2,'MarkerFaceColor',[1 0 0]);
            % Fish bones (visual aid for rotations)
            Null = zeros(constant.str.Ns+1,1);
            LE =   0.5*[C Null Null];
            TE = - 0.5*[C Null Null];
            for i=1:constant.str.Ns+1
                R = expon(P(i,4:6));
                LErot = (R*LE(i,:)')' + P0(i,:) + P(i,1:3);
                TErot = (R*TE(i,:)')' + P0(i,:) + P(i,1:3);
                plot3([LErot(:,1), P0(i,1)+P(i,1)],...
                    [LErot(:,2), P0(i,2)+P(i,2)],...
                    [LErot(:,3), P0(i,3)+P(i,3)],'g.-');
                plot3([TErot(:,1), P0(i,1)+P(i,1)],...
                    [TErot(:,2), P0(i,2)+P(i,2)],...
                    [TErot(:,3), P0(i,3)+P(i,3)],'g.-');
            end
            % Axes properties
            xlabel('Chord [m]')
            ylabel('Span [m]')
            zlabel('Height [m]')
            view(-50,25)
            axis equal
        end
    end
end
