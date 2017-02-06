function [statics,exitflag] = strSolve(constant,statics,Vi,Vf,alphar,ders,lin,tailflag,varargin)

morphflag = 0;

%=========================================================================%
% Initialization
%=========================================================================%
% Initialisation of the NR parameters
Vind     = 0;

if Vi == Vf && Vf == 0
    dsc = 0.25;
else
    if Vi == Vf
        deltaV = min(1,0.1*Vf);
        Vi = Vf-deltaV;
    end
    if Vf>Vi
        step = ceil(4*(Vf-Vi)/Vf);
    else
        step = ceil(4*(Vi-Vf)/Vi);
    end
    dsc = 1/step;
end
sc       = 0;

exitflag = 0;
tol1     = 1e-7; % Inner loop tollerance
tol2     = 1e-4; % Tollerance on outer loop step size
ppecho.print = 0;
ppecho.plot  = 0;

% Extract parameters from structured arrays
frdof = constant.str.frdof;

%=========================================================================%
% AEROELASTIC SOLUTIONS
%=========================================================================%
% NOT TRIMMED
%=========================================================================%
% if ppecho.print==1
    fprintf('Analysis run.\n')
% end
% Start outer loop
if lin == 1
    %=================================================================%
    % LINEAR STRUCTURAL SOLUTION
    %=================================================================%
    statics.str.p = zeros(size(statics.str.p));
    statics = structure(constant,statics,ders,tailflag,morphflag,0);
    statics = pgen_str(constant,statics,ders,morphflag);
    statics = fext (constant,statics,1,ders,0,tailflag,0); % Routine to create stiffness matrix from external forces
    R       = statics.str.Fext; % Add aerodynamic forces
    Rr      = R(frdof);
    J       = statics.str.Ks - statics.str.Kfext;
    Jr      = J(frdof,frdof);
    plin    = Jr\Rr;

    statics.str.p(frdof,1) = plin;

    %%% Update of the structural, aerodynamic and external forces based on the
    %%% structural deformations
    statics.str.Fs   = statics.str.Ks*statics.str.p;
    statics.str.Fext = statics.str.Fext + statics.str.Kfext*statics.str.p;

    % Compute the corresponding internal forces
    statics = structure_lin(constant,statics,ders);

    if ppecho.print==1
        fprintf('Speed [m/s]  Tip displacement [m]\n')
        fprintf('%3.1f              %5.3f',Vf,statics.str.p(constant.str.dof.wing.dz(end))) 
        fprintf('\n')
    end
    

    if ders==1
        for i = 1:numel(statics.str.C)
            dKsdC = reshape(statics.sens.pdC_Ks(:,i),size(statics.str.Ks,2),size(statics.str.Ks,1))';
            dpdC(:,i) = -(statics.str.Ks(frdof,frdof)-statics.str.Kfext(frdof,frdof))\dKsdC(frdof,frdof)*plin;
            dpfulldC = zeros(length(statics.str.p),1);
            dpfulldC(frdof) = dpdC(:,i);
            dFsdC(:,i) = statics.str.Ks*dpfulldC+dKsdC*statics.str.p;
        end

        statics.sens.dC_p  = zeros(length(statics.str.p),numel(statics.str.C));
        statics.sens.dC_p(frdof,:) = dpdC;
        statics.sens.dC_Ks = statics.sens.pdC_Ks;
        statics.sens.dC_Fs = dFsdC;
        statics.sens.dC_Fext = statics.str.Kfext*statics.sens.dC_p;
        statics.sens.dC_re = statics.sens.pdC_re+statics.sens.dp_re*statics.sens.dC_p;

        for i = 1:size(statics.sens.dKfextdt,2)
            dKextdt = reshape(statics.sens.dKfextdt(:,i),size(statics.str.Kfext,2),size(statics.str.Kfext,1))';
            dpdt(:,i) = (statics.str.Ks(frdof,frdof)-statics.str.Kfext(frdof,frdof))\(statics.sens.dFextdt(frdof,i)+dKextdt(frdof,frdof)*plin);

            dpfulldt = zeros(length(statics.str.p),1);
            dpfulldt(frdof) = dpdt(:,i);
            dFextdt(:,i) = statics.sens.dFextdt(:,i)+statics.str.Kfext*dpfulldt+dKextdt*statics.str.p;
        end
        statics.sens.dpdt  = zeros(length(statics.str.p),size(statics.sens.dKfextdt,2));
        statics.sens.dpdt(frdof,:) = dpdt;
        statics.sens.dFsdt = statics.str.Ks*statics.sens.dpdt;
        statics.sens.dredt = statics.sens.dp_re*statics.sens.dpdt;

        statics.sens.dKfextdt = statics.sens.dKfextdt;
        statics.sens.dFextdt = dFextdt;
    end
else
    %=================================================================%
    % NON-LINEAR SOLUTION
    %=================================================================%
    if ppecho.print==1
        fprintf('iter    Speed [m/s]  Tip displacement [m]     check\n')
    end

    %%% Initialize for predictor step
    statics = structure(constant,statics,0,tailflag,morphflag,0);
    statics = pgen_str(constant,statics,0,morphflag);
    statics = fext(constant,statics,sc,0,0,tailflag,alphar,0);
    J  = statics.str.Ks - statics.str.Kfext;
    Jr = J(constant.str.frdof,constant.str.frdof); 
    R  = statics.str.Fext - statics.str.Fs;
    Rr = R(constant.str.frdof);
    count = 0;

    while(sc<1)
        count = count + 1;
        Vind = Vind+1;
        scold = sc;
        sc = sc +dsc;
        Vsq = Vi^2 + sc*(Vf^2-Vi^2);
        V = sqrt(Vsq);
        if sc>1
            V  = Vf;
            sc = 1;
        end

        check = 1;
        iter  = 0;

        % Predictor step
        statics.str.p(constant.str.frdof,1) = statics.str.p(constant.str.frdof,1) +...
                Jr\(Rr + statics.sens.dFextdsc(constant.str.frdof,1)*(sc-scold)); 

        % Start inner loop
        while(check>tol1)
            iter    = iter + 1;
            % visual check
            statics = structure(constant,statics,0,tailflag,morphflag,0);
            statics = pgen_str(constant,statics,0,morphflag);     
            statics = fext(constant,statics,sc,0,0,tailflag,alphar,0);
            R = statics.str.Fext - statics.str.Fs;
            Rr = R(constant.str.frdof);
            J  = statics.str.Ks - statics.str.Kfext; 
            Jr = J(constant.str.frdof,constant.str.frdof);
            dd = Jr\Rr;
            statics.str.p(constant.str.frdof,1) = statics.str.p(constant.str.frdof,1)+dd; 
            check = norm(dd,inf)/norm(statics.str.p,inf);
            % If inner loop takes too long
            if iter>20
                statics.str.p = dold;
                sc = sc - dsc;
                dsc = dsc/2;
                fprintf('Reducing step size to %5.3e \n',dsc)
                sc = sc + dsc;
                Vsq = Vi^2 + sc*(Vf^2-Vi^2);
                V = sqrt(Vsq);
                iter = 0;
            end
            if abs(dsc)<tol2
                fprintf('Step size too small, no solution found \n')
                check    = 1e-9;
                exitflag = 1;
            end
        end

        % End inner loop
        % Plot options
        if ppecho.plot==1
            hold off,pause(.1),plotcfg(elm,elma),pause(.1)
            sb = statusbar('Speedup without trim, iter is %i',iter);
            set(sb(1).ProgressBar,'Visible','on','Value',V/Vf*100);
        elseif ppecho.print==1
            fprintf(' %2i       %3.1f              %5.3f              %4.2e',iter,V,statics.str.p(end-3),check)
            fprintf('\n')
        end

        % Variable step size
        if iter<3 && exitflag~=1 && sc~=1
            dsc = dsc*2;
            fprintf('Great convergence, increasing step size to %5.3e \n',dsc)
        elseif iter>8
            dsc = dsc/2;
            fprintf('Decreasing step size to %5.3e \n',dsc)
        end

        % If no solution is found, exit the loop
        if exitflag==1
            break
        end
        
        if ders==1
            % Symmetric aeroelasticity
            fprintf('Calculating sensitivities wrt C\n')
            statics = structure(constant,statics,1,tailflag,morphflag,0);
            statics = fext(constant,statics,1,ders,1,tailflag,alphar);
            
            dC_p = -(statics.str.Ks(frdof,frdof)-statics.str.Kfext(frdof,frdof))\statics.sens.pdC_Fs(frdof,:);
            dpdt = (statics.str.Ks(frdof,frdof)-statics.str.Kfext(frdof,frdof))\statics.sens.dFextdt(frdof,:);
            
            statics.sens.dC_p  = zeros(length(statics.str.p),numel(statics.str.C));
            statics.sens.dC_p(frdof,:) = dC_p;
            statics.sens.dpdt  = zeros(length(statics.str.p),size(dpdt,2));
            statics.sens.dpdt(frdof,:) = dpdt;
            
            statics.sens.dC_Fs = statics.sens.pdC_Fs+statics.sens.dp_Fs*statics.sens.dC_p;
            statics.sens.dC_Ks = statics.sens.pdC_Ks+statics.sens.dp_Ks*statics.sens.dC_p;
            statics.sens.dC_re = statics.sens.pdC_re+statics.sens.dp_re*statics.sens.dC_p;
            
            statics.sens.dFsdt = statics.sens.dp_Fs*statics.sens.dpdt;
            statics.sens.dKsdt = statics.sens.dp_Ks*statics.sens.dpdt;
            statics.sens.dredt = statics.sens.dp_re*statics.sens.dpdt;
            
            statics.sens.dC_Fext = statics.str.Kfext*statics.sens.dC_p;
            statics.sens.dC_Kfext = statics.sens.dKfextdp*statics.sens.dC_p;
            
            statics.sens.dKfextdt = statics.sens.dKfextdt+statics.sens.dKfextdp*statics.sens.dpdt;
            statics.sens.dFextdt = statics.sens.dFextdt+statics.str.Kfext*statics.sens.dpdt;
            
            statics = pgen_str(constant,statics,1,morphflag);
        else
            statics = structure(constant,statics,0,tailflag,morphflag,0);
            statics = pgen_str(constant,statics,0,morphflag);
            statics = fext(constant,statics,1,0,0,tailflag,alphar,0);
        end
        
    end
end
%=========================================================================%
% END
%=========================================================================%
