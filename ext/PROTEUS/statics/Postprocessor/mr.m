function statics = mr(constant,statics,ders,tailflag,morphflag)

R     = statics.str.Fext-statics.str.Fs;

statics.str.Mr = -R(4:6);
% statics.str.L = statics.aero.Lift;
if ders==1 
    if tailflag == 1
        dC_R = statics.sens.dC_Fext-statics.sens.dC_Fs;
        dRdt = statics.sens.dFextdt-statics.sens.dFsdt;
        
        statics.sens.dC_Mr = -dC_R(4:6,:);
        statics.sens.dMrdt = -dRdt(4:6,:);
    end
end