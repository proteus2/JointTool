function [constant,lampar,stringer,BlendingCst,crossmod,feas] = analysis_crossmod_preproc(constant,dv,ANALYSISTYPE,ders)

curdir = cd;

if ANALYSISTYPE == 1 || ANALYSISTYPE == 3
    tailflag = 1;
else
    tailflag = 0;
end

if ANALYSISTYPE == 2 || ANALYSISTYPE == 3
    morphflag = 1;
else
    morphflag = 0;
end

%% Lamination Parameters

cd('ext/PROTEUS/lam_param')
feas = 1;
if ders == 1 && tailflag == 1
    lampar.dAddv = sparse(6*length(constant.lam.ID),9*length(constant.lam.ID));
    lampar.dDddv = sparse(6*length(constant.lam.ID),9*length(constant.lam.ID));
    lampar.dtddv = sparse(length(constant.lam.ID),9*length(constant.lam.ID));
    lampar.dg1ddv = sparse(length(constant.lam.ID),9*length(constant.lam.ID));
    lampar.dg2ddv = sparse(length(constant.lam.ID),9*length(constant.lam.ID));
    lampar.dg3ddv = sparse(length(constant.lam.ID),9*length(constant.lam.ID));
    lampar.dg4ddv = sparse(length(constant.lam.ID),9*length(constant.lam.ID));
    lampar.dg5ddv = sparse(length(constant.lam.ID),9*length(constant.lam.ID));
    lampar.dg6ddv = sparse(length(constant.lam.ID),9*length(constant.lam.ID));
end
for i= 1:length(constant.lam.ID)
    lamparcell = lam_param(constant.mat,constant.lam.matID(i),dv(9*i),dv(9*(i-1)+(1:8)),ders);
    lampar.A((i-1)*3+(1:3),(1:3)) = lamparcell.A;
    lampar.D((i-1)*3+(1:3),(1:3)) = lamparcell.D;
    
    if tailflag == 1 
        lampar.g1(i,1) = lamparcell.g1;
        lampar.g2(i,1) = lamparcell.g2;
        lampar.g3(i,1) = lamparcell.g3;
        lampar.g4(i,1) = lamparcell.g4;
        lampar.g5(i,1) = lamparcell.g5;
        lampar.g6(i,1) = lamparcell.g6;
        
        if lampar.g1(i,1)<-1e-12 || lampar.g2(i,1)<-1e-12 || lampar.g3(i,1)<-1e-12 || lampar.g4(i,1)<-1e-12
            display('Infeasible design, switch to linear solver')
            feas = 0;
            lin  = 1;
            if strcmp(constant.opt.Optimiser,'None')
                error(['Infeasible lamination param. detected during analysis, Check Lam. #',num2str(i)])
            end
        end
    end
    lampar.t(i) = dv(9*i);
    if ders == 1 && tailflag == 1
        lampar.dAddv((i-1)*6+(1:6),(i-1)*9+(1:9)) = lamparcell.dAddv;
        lampar.dDddv((i-1)*6+(1:6),(i-1)*9+(1:9)) = lamparcell.dDddv;
        lampar.dtddv(i,(i-1)*9+9) = 1;
        lampar.dg1ddv(i,(i-1)*9+(1:9)) = lamparcell.dg1(1:9);
        lampar.dg2ddv(i,(i-1)*9+(1:9)) = lamparcell.dg2(1:9);
        lampar.dg3ddv(i,(i-1)*9+(1:9)) = lamparcell.dg3(1:9);
        lampar.dg4ddv(i,(i-1)*9+(1:9)) = lamparcell.dg4(1:9);
        lampar.dg5ddv(i,(i-1)*9+(1:9)) = lamparcell.dg5(1:9);
        lampar.dg6ddv(i,(i-1)*9+(1:9)) = lamparcell.dg6(1:9);
    end
end

cd(curdir)
fprintf(' Lamination parameters ;')

%% Stringer laminates

if isfield(constant,'stringer')

    cd('ext/PROTEUS/lam_param')
    for i= 1:length(constant.stringer.lamID)
        stringercell = lam_param_stringer(constant.mat,constant.stringer.matID(i),constant.stringer.EA(i),constant.stringer.h(i),zeros(8,1));
        stringer.A((i-1)*3+(1:3),(1:3)) = stringercell.A;
        stringer.D((i-1)*3+(1:3),(1:3)) = stringercell.D;
        
        stringer.t(i) = stringercell.t;
    end
    
    cd(curdir)
    fprintf(' Stringer laminates ;')
else
    stringer = [];
end

%% Cross-sectional Modeller

cd('ext/PROTEUS/cross_mod')
[crossmod,stringer] = cross_mod(constant,lampar,stringer,1,ders,tailflag,morphflag);
cd(curdir)
fprintf(' Cross-sectional modeller ;')

%% Blending Constraints
if constant.opt.BlendConst && tailflag == 1
    cd Blending
    gblend  = [];
    dgblend = [];

    BlendOptions = constant.opt.BlendOptions;
    
    [BlendConstTopSkin ,BlendConstBotSkin ,BlendConstSpar,...
        DBlendConstTopSkin,DBlendConstBotSkin,DBlendConstSpar] = ComputeBlendingConstraints (constant,dv,BlendOptions);
        
    if BlendOptions.TopSkin
        gblend = [gblend; BlendConstTopSkin.IP(find(BlendConstTopSkin.IP));
            BlendConstTopSkin.OOP(find(BlendConstTopSkin.OOP ))];           %#ok<FNDSB>
    end
    
    if  BlendOptions.BotSkin
        gblend = [gblend; BlendConstBotSkin.IP(find(BlendConstBotSkin.IP));
            BlendConstBotSkin.OOP(find(BlendConstBotSkin.OOP))];            %#ok<FNDSB>
    end
    
    if BlendOptions.Spar
        for iSpar = 1 : length(BlendConstSpar)
            Const = BlendConstSpar{iSpar};
            gblend = [gblend; Const.IP(find(Const.IP)); Const.OOP(find(Const.OOP))];  %#ok<FNDSB>
        end
    end
    
    if ders == 1
        LamTopSkin = cell2mat(constant.lam.TopID(:));          % ID of Unique laminates TopSkin
        LamBotSkin = cell2mat(constant.lam.BotID(:));
        LamSpar    = cell2mat(constant.lam.SparID(:));
        
        if BlendOptions.TopSkin
            dgblend = AddBlending_localFct(LamTopSkin,dgblend, BlendConstTopSkin.IP',  DBlendConstTopSkin.IP,'IP');
            dgblend = AddBlending_localFct(LamTopSkin,dgblend, BlendConstTopSkin.OOP', DBlendConstTopSkin.OOP,'OOP');
        end
        
        if BlendOptions.BotSkin
            dgblend = AddBlending_localFct(LamBotSkin,dgblend, BlendConstBotSkin.IP',  DBlendConstBotSkin.IP,'IP');
            dgblend = AddBlending_localFct(LamBotSkin,dgblend, BlendConstBotSkin.OOP', DBlendConstBotSkin.OOP,'OOP');
        end
        
        if BlendOptions.Spar
            SparNum     = cell2mat(constant.lam.SparNum(:)')';
            SparNumbers = unique(SparNum);
            for iSpar = 1 : length(BlendConstSpar)
                indexes = find(SparNumbers(iSpar)==SparNum);
                dgblend = AddBlending_localFct(LamSpar(indexes),dgblend, BlendConstSpar{iSpar}.IP',  DBlendConstSpar{iSpar}.IP,'IP');
                dgblend = AddBlending_localFct(LamSpar(indexes),dgblend, BlendConstSpar{iSpar}.OOP', DBlendConstSpar{iSpar}.OOP,'OOP');
            end
        end
        BlendingCst.dgblend = dgblend;
    end
    BlendingCst.gblend = gblend;
    
    cd (constant.curdir)
else
    BlendingCst.gblend  = [];
    BlendingCst.dgblend = [];
end


%% Local Blending function
function [dg] = AddBlending_localFct(LamID,dg,Const,DConst,String)
Sizedg = size(dg);
index  = 1;
for j = 1:length(LamID)-1
    index1 = LamID(j);
    T_sec1  = 9*(index1)         ;
    LP_sec1 = 9*(index1-1)+(1:8) ;
    
    for i = j+1 : length(LamID)
        index2  = LamID(i);
        T_sec2  = 9*(index2)         ;
        LP_sec2 = 9*(index2-1)+(1:8) ;
        
        if Const(j,i)~=0
            if strcmp(String,'IP')==1
                dg(Sizedg(1) + index,LP_sec1(1:4))    = DConst{j,i}(1:4);
                dg(Sizedg(1) + index,LP_sec2(1:4))    = DConst{j,i}(5:8);
            end
            if strcmp(String,'OOP')==1
                dg(Sizedg(1) + index,LP_sec1(5:8))    = DConst{j,i}(1:4);
                dg(Sizedg(1) + index,LP_sec2(5:8))    = DConst{j,i}(5:8);
            end
            dg(Sizedg(1) + index,[T_sec1 T_sec2]) = DConst{j,i}(9:10);
            index = index +1;
        end
        
    end
end

