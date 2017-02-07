classdef structuralModel < handle
   properties 
      grid  
      Ks
      Kfext
      Fs
      Fext
      p
      frdof
      fxdof
   end
   methods
       function obj = structuralModel(model)
           [constant,~,~,~,statics] = strModule(model);
           % Extract properties
           obj.grid  = reshape(constant{1}.str.xyz,3,constant{1}.str.Ns+1)';
           obj.Ks    = statics{1}.str.Ks;
           obj.Kfext = statics{1}.str.Kfext;
           obj.Fs    = reshape(statics{1}.str.Fs,6,constant{1}.str.Ns+1)';
           obj.Fext  = reshape(statics{1}.str.Fext,6,constant{1}.str.Ns+1)';
           obj.p     = reshape(statics{1}.str.p,6,constant{1}.str.Ns+1)';
           obj.frdof = constant{1}.str.frdof;
           obj.fxdof = constant{1}.str.fxdof;
           obj.fxdof = constant{1}.str.fxdof;
           % Save array
           curdir = cd;
           cd('ext/PROTEUS/results')
           save('constant.mat','constant')
           cd(curdir)
       end
       function obj = solve(obj,F)
           % Format 
           % (Nx, Ny, Nz, Mx, My, Mz)_per node
           obj.p = obj.Ks(obj.frdof,obj.frdof)\F(obj.frdof,1);
       end
       function obj = plot(obj,varargin)
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
       function obj = plotT3(obj)
           figure()
           hold on
           grid on
           plot(obj.grid(:,2),obj.p(:,3),'k');
           xlabel('Span [m]')
           ylabel('Vertical Displacement [m]')
       end
       function obj = plotR2(obj)
           figure()
           hold on
           grid on
           plot(obj.grid(:,2),rad2deg(obj.p(:,5)),'k')
           xlabel('Span [m]')
           ylabel('Torsional Deflection [deg]')
       end
   end
end