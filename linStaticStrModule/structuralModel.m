classdef structuralModel
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
       end
       function obj = solve(obj,F)
           % Format 
           % (Nx, Ny, Nz, Mx, My, Mz)_per node
           obj.p = obj.Ks(obj.frdof,obj.frdof)\F(obj.frdof,1);
       end
   end
end