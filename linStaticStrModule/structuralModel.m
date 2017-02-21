classdef structuralModel
   properties 
      grid  
      disp
      Fs
      settings = strModelSettings();
   end
   properties (SetAccess = protected)
      Fext
   end
   properties (Access = private)
      Ks
      Kfext
      frdof
      fxdof 
   end
   methods % Constructor
       function obj = structuralModel(model)
           [constant,~,~,~,statics] = genStrModel(model);
           % Extract properties
           obj.grid  = reshape(constant{1}.str.xyz,3,constant{1}.str.Ns+1)';
           obj.Ks    = statics{1}.str.Ks;
           obj.Kfext = statics{1}.str.Kfext;
           obj.Fs    = statics{1}.str.Fs;
           obj.Fext  = statics{1}.str.Fext;
           obj.disp  = statics{1}.str.p;
           obj.frdof = constant{1}.str.frdof;
           obj.fxdof = constant{1}.str.fxdof;
           obj.fxdof = constant{1}.str.fxdof;
           % Save array
           curdir = cd;
           cd('ext/PROTEUS/results')
           save('constant.mat','constant')
           cd(curdir)
       end
   end
   methods % Set
       function obj = set.Fext(obj,Fext)
           obj.Fext = Fext;
       end
   end
   methods % Solution
       function obj = solve(obj)
           % Format 
           % (Nx, Ny, Nz, Mx, My, Mz)_per node
           obj.disp = [];
           obj.disp(obj.frdof,1) = obj.Ks(obj.frdof,obj.frdof)\obj.Fs(obj.frdof,1);
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
           C  = constant{1}.inp.cbox;
           % Beam axis
           plot3(P0(:,1),P0(:,2),P0(:,3),...
                 'ro-','LineWidth',2,'MarkerFaceColor',[1 0 0]);
           % Fish bones (visual aid for rotations)
           Null = zeros(constant{1}.str.Ns+1,1);
           xref = constant{1}.inp.xref;
           LE = P0 + [xref   Null Null].*[C Null Null];
           TE = P0 + [xref-1 Null Null].*[C Null Null];
           for i=1:constant{1}.str.Ns+1
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
              P = reshape(P,6,constant{1}.str.Ns+1)'; 
           end
           C  = constant{1}.inp.cbox;
           % Deformed beam axis
           plot3(P0(:,1)+P(:,1),...
                 P0(:,2)+P(:,2),...
                 P0(:,3)+P(:,3),'ro-','LineWidth',2,'MarkerFaceColor',[1 0 0]);
           % Fish bones (visual aid for rotations)
           Null = zeros(constant{1}.str.Ns+1,1);
           LE =   0.5*[C Null Null];
           TE = - 0.5*[C Null Null];
           for i=1:constant{1}.str.Ns+1
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