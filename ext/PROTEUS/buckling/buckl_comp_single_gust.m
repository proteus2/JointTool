function [buckl] = buckl_comp_single_gust(constant,lampar,Numcross,strainraw,ders,gustflag,tailflag,varargin)

if length(varargin) == 1
    morph = varargin{1};
    morphflag = 1;
else
    morphflag = 0;
end

buckl_loc = constant.buckl.yloc;

if ders == 1
    if tailflag == 1
        ndA = size(strainraw{1}.dexoutdA{1},2);
        ndD = size(strainraw{1}.dexoutdD{1},2);
        ndt = size(strainraw{1}.dexoutdt{1},2);
        if gustflag == 1
            nddv = size(strainraw{1}.dexoutddv{1},2);
        end
    end
    if morphflag == 1
        if morph.camber
            ndparam = size(strainraw{1}.dexoutdparam{1},2);
        end
        if morph.twist
            ndphi = size(strainraw{1}.dexoutdphi{1},2);
        end
        if morph.fold
            ndtheta = size(strainraw{1}.dexoutdtheta{1},2);
        end
        if morph.shear
            ndpsi = size(strainraw{1}.dexoutdpsi{1},2);
        end
        if morph.span
            ndext = size(strainraw{1}.dexoutdext{1},2);
        end
    end
end

% Note if no gust, Numgust is 1 and will only take the static loads into
% account

lammat = [];
bucklsec = [];
bucklcross = [];
for i=1:length(constant.buckl.lam)
    for j=1:length(constant.buckl.lam{i});
        lammat = [lammat;constant.buckl.lam{i}{j}];
        bucklsec = [bucklsec;i*ones(length(constant.buckl.lam{i}{j}),1)];
        bucklcross = [bucklcross;j*ones(length(constant.buckl.lam{i}{j}),1)];
    end
end

ntot = 0;
for i=1:length(buckl_loc)-1
    
    nsec1 = constant.buckl.nsec1{i};
    nsec2 = constant.buckl.nsec2{i};
    
    if nsec1 == nsec2% Buckling panel covers a single element
        for ncross=1:size(constant.buckl.yzpanel{i},2)
            numcross = Numcross{nsec1}(ncross);
            ntot = ntot + 1;
            for k=1:constant.buckl.Npan{i}(1,ncross)
                
                ABD = [lampar.A(3*(constant.buckl.lam{i}{ncross}(k,1)-1)+(1:3),:) zeros(3); zeros(3) lampar.D(3*(constant.buckl.lam{i}{ncross}(k,1)-1)+(1:3),:)];
                
                if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                    evecelm1 = [0*strainraw{numcross}.exout(constant.buckl.crosselm{i}{ncross}{k},1)';0*strainraw{numcross}.eyout(constant.buckl.crosselm{i}{ncross}{k},1)';strainraw{numcross}.gammaout(constant.buckl.crosselm{i}{ncross}{k},1)'];
                    evecelm2 = [0*strainraw{numcross}.exout(constant.buckl.crosselm{i}{ncross}{k},2)';0*strainraw{numcross}.eyout(constant.buckl.crosselm{i}{ncross}{k},2)';strainraw{numcross}.gammaout(constant.buckl.crosselm{i}{ncross}{k},2)'];
                else
                    evecelm1 = [strainraw{numcross}.exout(constant.buckl.crosselm{i}{ncross}{k},1)';strainraw{numcross}.eyout(constant.buckl.crosselm{i}{ncross}{k},1)';strainraw{numcross}.gammaout(constant.buckl.crosselm{i}{ncross}{k},1)'];
                    evecelm2 = [strainraw{numcross}.exout(constant.buckl.crosselm{i}{ncross}{k},2)';strainraw{numcross}.eyout(constant.buckl.crosselm{i}{ncross}{k},2)';strainraw{numcross}.gammaout(constant.buckl.crosselm{i}{ncross}{k},2)'];
                end
                clear evec1 evec2
                evec1(1,:) = interp1(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1,:);evecelm2(1,:)],buckl_loc(i));
                evec1(2,:) = interp1(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2,:);evecelm2(2,:)],buckl_loc(i));
                evec1(3,:) = interp1(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3,:);evecelm2(3,:)],buckl_loc(i));
                
                evec2(1,:) = interp1(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1,:);evecelm2(1,:)],buckl_loc(i+1));
                evec2(2,:) = interp1(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2,:);evecelm2(2,:)],buckl_loc(i+1));
                evec2(3,:) = interp1(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3,:);evecelm2(3,:)],buckl_loc(i+1));
                
                N1 = ABD(1:3,1:3)*evec1;
                N2 = ABD(1:3,1:3)*evec2;
                
                clear r1 r2 dr1dDpart dr2dDpart dr1dN dr2dN
                for l=1:size(N1,2)
                    [r1(:,l),dr1dDpart{l},dr1dN{l}] = quadbuckl(N1(:,l),ABD(4:6,4:6),2,constant.buckl.stiff{i}{ncross}{k});
                    [r2(:,l),dr2dDpart{l},dr2dN{l}] = quadbuckl(N2(:,l),ABD(4:6,4:6),2,constant.buckl.stiff{i}{ncross}{k});
                end
                
                % Identify critical modes
                if size(r1,2)>1
                    [~,ind1] = sort(r1(1,:),'descend');
                    [~,ind2] = sort(r2(1,:),'descend');
                    
                    r = [r1(:,ind1(1));r1(:,ind1(2));r2(:,ind2(1));r2(:,ind2(2))];
                else
                    r = [r1;r1;r2;r2];
                end
                
                rcell{i}{ncross}(k,:) = r';
                rvec{ntot}(k,:) = r';
                if ders == 1
                    if size(r1,2)>1
                        if tailflag == 1
                            dN1dA = zeros(6,ndA);
                            dN1dD = zeros(6,ndD);
                            dN1dt = zeros(6,ndt);
                            
                            dN2dA = zeros(6,ndA);
                            dN2dD = zeros(6,ndD);
                            dN2dt = zeros(6,ndt);
                            
                            dr1dA = zeros(4,ndA);
                            dr1dD = zeros(4,ndD);
                            dr1dt = zeros(4,ndt);
                            
                            dr2dA = zeros(4,ndA);
                            dr2dD = zeros(4,ndD);
                            dr2dt = zeros(4,ndt);
                            
                            if gustflag == 1
                                dN1ddv = zeros(6,nddv);
                                dN2ddv = zeros(6,nddv);
                                dr1ddv = zeros(4,nddv);
                                dr2ddv = zeros(4,nddv);
                            end
                            
                            dABDdA = zeros(9,ndA);
                            for m=1:2
                                % Loads on first end of the element
                                if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                    devecelm1dA = [0*strainraw{numcross}.dexoutdA{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdA{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdA{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    devecelm2dA = [0*strainraw{numcross}.dexoutdA{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdA{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdA{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    
                                    devecelm1dD = [0*strainraw{numcross}.dexoutdD{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdD{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdD{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    devecelm2dD = [0*strainraw{numcross}.dexoutdD{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdD{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdD{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    
                                    devecelm1dt = [0*strainraw{numcross}.dexoutdt{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdt{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdt{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    devecelm2dt = [0*strainraw{numcross}.dexoutdt{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdt{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdt{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    if gustflag == 1
                                        devecelm1ddv = [0*strainraw{numcross}.dexoutddv{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutddv{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutddv{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                        devecelm2ddv = [0*strainraw{numcross}.dexoutddv{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutddv{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutddv{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    end
                                else
                                    devecelm1dA = [strainraw{numcross}.dexoutdA{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdA{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdA{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    devecelm2dA = [strainraw{numcross}.dexoutdA{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdA{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdA{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    
                                    devecelm1dD = [strainraw{numcross}.dexoutdD{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdD{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdD{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    devecelm2dD = [strainraw{numcross}.dexoutdD{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdD{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdD{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    
                                    devecelm1dt = [strainraw{numcross}.dexoutdt{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdt{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdt{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    devecelm2dt = [strainraw{numcross}.dexoutdt{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdt{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdt{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    if gustflag == 1
                                        devecelm1ddv = [strainraw{numcross}.dexoutddv{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutddv{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutddv{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                        devecelm2ddv = [strainraw{numcross}.dexoutddv{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutddv{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutddv{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    end
                                end
                                
                                
                                [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1,ind1(m));evecelm2(1,ind1(m))],buckl_loc(i),1,2);
                                [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2,ind1(m));evecelm2(2,ind1(m))],buckl_loc(i),1,2);
                                [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3,ind1(m));evecelm2(3,ind1(m))],buckl_loc(i),1,2);
                                
                                devec1dA(1,:) =  devec1devecelm1*[devecelm1dA(1,:);devecelm2dA(1,:)];
                                devec1dA(2,:) =  devec1devecelm2*[devecelm1dA(2,:);devecelm2dA(2,:)];
                                devec1dA(3,:) =  devec1devecelm3*[devecelm1dA(3,:);devecelm2dA(3,:)];
                                
                                devec1dD(1,:) =  devec1devecelm1*[devecelm1dD(1,:);devecelm2dD(1,:)];
                                devec1dD(2,:) =  devec1devecelm2*[devecelm1dD(2,:);devecelm2dD(2,:)];
                                devec1dD(3,:) =  devec1devecelm3*[devecelm1dD(3,:);devecelm2dD(3,:)];
                                
                                devec1dt(1,:) =  devec1devecelm1*[devecelm1dt(1,:);devecelm2dt(1,:)];
                                devec1dt(2,:) =  devec1devecelm2*[devecelm1dt(2,:);devecelm2dt(2,:)];
                                devec1dt(3,:) =  devec1devecelm3*[devecelm1dt(3,:);devecelm2dt(3,:)];
                                
                                if gustflag == 1
                                    devec1ddv(1,:) =  devec1devecelm1*[devecelm1ddv(1,:);devecelm2ddv(1,:)];
                                    devec1ddv(2,:) =  devec1devecelm2*[devecelm1ddv(2,:);devecelm2ddv(2,:)];
                                    devec1ddv(3,:) =  devec1devecelm3*[devecelm1ddv(3,:);devecelm2ddv(3,:)];
                                end
                                
                                row    = [1,1,1,2,2,3];     % Row number of independent A and D matrices elements
                                column = [1,2,3,2,3,3];     % Column number of independent A and D matrices elements
                                
                                for l = 1:6
                                    
                                    dABDdAloc = zeros(3,3);
                                    dABDdAloc(row(l),column(l)) = 1;
                                    dABDdAloc(column(l),row(l)) = 1;
                                    
                                    dABDdA(:,(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = reshape(dABDdAloc',[],1);
                                    
                                    dN1dA(3*(m-1)+(1:3),(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = dABDdAloc*evec1(:,ind1(m));
                                    
                                    if row(l)==column(l)
                                        dr1dD(2*(m-1)+(1:2),(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = [dr1dDpart{ind1(m)}(row(l),column(l),1);dr1dDpart{ind1(m)}(row(l),column(l),2)];
                                    else
                                        dr1dD(2*(m-1)+(1:2),(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = [dr1dDpart{ind1(m)}(row(l),column(l),1);dr1dDpart{ind1(m)}(row(l),column(l),2)]+[dr1dDpart{ind1(m)}(column(l),row(l),1);dr1dDpart{ind1(m)}(column(l),row(l),2)];
                                    end
                                end
                                
                                dN1dA(3*(m-1)+(1:3),:)=dN1dA(3*(m-1)+(1:3),:)+ABD(1:3,1:3)*devec1dA;
                                dN1dD(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1dD;
                                dN1dt(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1dt;
                                
                                if gustflag == 1
                                    dN1ddv(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1ddv;
                                end
                                
                                dr1dA(2*(m-1)+(1:2),:)=dr1dA(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dA(3*(m-1)+(1:3),:);
                                dr1dD(2*(m-1)+(1:2),:)=dr1dD(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dD(3*(m-1)+(1:3),:);
                                dr1dt(2*(m-1)+(1:2),:)=dr1dt(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dt(3*(m-1)+(1:3),:);
                                
                                if gustflag == 1
                                    dr1ddv(2*(m-1)+(1:2),:)=dr1ddv(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1ddv(3*(m-1)+(1:3),:);
                                end
                                
                                % Loads on second end of the element
                                if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                    devecelm1dA = [0*strainraw{numcross}.dexoutdA{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdA{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdA{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    devecelm2dA = [0*strainraw{numcross}.dexoutdA{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdA{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdA{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    
                                    devecelm1dD = [0*strainraw{numcross}.dexoutdD{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdD{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdD{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    devecelm2dD = [0*strainraw{numcross}.dexoutdD{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdD{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdD{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    
                                    devecelm1dt = [0*strainraw{numcross}.dexoutdt{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdt{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdt{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    devecelm2dt = [0*strainraw{numcross}.dexoutdt{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdt{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdt{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    
                                    if gustflag == 1
                                        devecelm1ddv = [0*strainraw{numcross}.dexoutddv{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutddv{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutddv{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                        devecelm2ddv = [0*strainraw{numcross}.dexoutddv{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutddv{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutddv{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    end
                                else
                                    devecelm1dA = [strainraw{numcross}.dexoutdA{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdA{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdA{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    devecelm2dA = [strainraw{numcross}.dexoutdA{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdA{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdA{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    
                                    devecelm1dD = [strainraw{numcross}.dexoutdD{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdD{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdD{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    devecelm2dD = [strainraw{numcross}.dexoutdD{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdD{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdD{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    
                                    devecelm1dt = [strainraw{numcross}.dexoutdt{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdt{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdt{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    devecelm2dt = [strainraw{numcross}.dexoutdt{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdt{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdt{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    
                                    if gustflag == 1
                                        devecelm1ddv = [strainraw{numcross}.dexoutddv{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutddv{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutddv{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                        devecelm2ddv = [strainraw{numcross}.dexoutddv{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutddv{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutddv{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    end
                                end
                                [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1,ind2(m));evecelm2(1,ind2(m))],buckl_loc(i+1),1,2);
                                [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2,ind2(m));evecelm2(2,ind2(m))],buckl_loc(i+1),1,2);
                                [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3,ind2(m));evecelm2(3,ind2(m))],buckl_loc(i+1),1,2);
                                
                                devec2dA(1,:) =  devec2devecelm1*[devecelm1dA(1,:);devecelm2dA(1,:)];
                                devec2dA(2,:) =  devec2devecelm2*[devecelm1dA(2,:);devecelm2dA(2,:)];
                                devec2dA(3,:) =  devec2devecelm3*[devecelm1dA(3,:);devecelm2dA(3,:)];
                                
                                devec2dD(1,:) =  devec2devecelm1*[devecelm1dD(1,:);devecelm2dD(1,:)];
                                devec2dD(2,:) =  devec2devecelm2*[devecelm1dD(2,:);devecelm2dD(2,:)];
                                devec2dD(3,:) =  devec2devecelm3*[devecelm1dD(3,:);devecelm2dD(3,:)];
                                
                                devec2dt(1,:) =  devec2devecelm1*[devecelm1dt(1,:);devecelm2dt(1,:)];
                                devec2dt(2,:) =  devec2devecelm2*[devecelm1dt(2,:);devecelm2dt(2,:)];
                                devec2dt(3,:) =  devec2devecelm3*[devecelm1dt(3,:);devecelm2dt(3,:)];
                                
                                if gustflag == 1
                                    devec2ddv(1,:) =  devec2devecelm1*[devecelm1ddv(1,:);devecelm2ddv(1,:)];
                                    devec2ddv(2,:) =  devec2devecelm2*[devecelm1ddv(2,:);devecelm2ddv(2,:)];
                                    devec2ddv(3,:) =  devec2devecelm3*[devecelm1ddv(3,:);devecelm2ddv(3,:)];
                                end
                                
                                row    = [1,1,1,2,2,3];     % Row number of independent A and D matrices elements
                                column = [1,2,3,2,3,3];     % Column number of independent A and D matrices elements
                                
                                for l = 1:6
                                    
                                    dABDdAloc = zeros(3,3);
                                    dABDdAloc(row(l),column(l)) = 1;
                                    dABDdAloc(column(l),row(l)) = 1;
                                    
                                    dABDdA(:,(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = reshape(dABDdAloc',[],1);
                                    
                                    dN2dA(3*(m-1)+(1:3),(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = dABDdAloc*evec2(:,ind2(m));
                                    
                                    if row(l)==column(l)
                                        dr2dD(2*(m-1)+(1:2),(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = [dr2dDpart{ind2(m)}(row(l),column(l),1);dr2dDpart{ind2(m)}(row(l),column(l),2)];
                                    else
                                        dr2dD(2*(m-1)+(1:2),(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = [dr2dDpart{ind2(m)}(row(l),column(l),1);dr2dDpart{ind2(m)}(row(l),column(l),2)]+[dr2dDpart{ind2(m)}(column(l),row(l),1);dr2dDpart{ind2(m)}(column(l),row(l),2)];
                                    end
                                end
                                
                                dN2dA(3*(m-1)+(1:3),:)=dN2dA(3*(m-1)+(1:3),:)+ABD(1:3,1:3)*devec2dA;
                                dN2dD(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2dD;
                                dN2dt(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2dt;
                                
                                if gustflag == 1
                                    dN2ddv(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2ddv;
                                end
                                
                                dr2dA(2*(m-1)+(1:2),:)=dr2dA(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dA(3*(m-1)+(1:3),:);
                                dr2dD(2*(m-1)+(1:2),:)=dr2dD(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dD(3*(m-1)+(1:3),:);
                                dr2dt(2*(m-1)+(1:2),:)=dr2dt(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dt(3*(m-1)+(1:3),:);
                                
                                if gustflag == 1
                                    dr2ddv(2*(m-1)+(1:2),:)=dr2ddv(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2ddv(3*(m-1)+(1:3),:);
                                end
                            end
                            
                            drdA = [dr1dA;dr2dA];
                            drdD = [dr1dD;dr2dD];
                            drdt = [dr1dt;dr2dt];
                            
                            if gustflag == 1
                                drddv = [dr1ddv;dr2ddv];
                            end
                            
                            %                                 drcelldA{i}{ncross}(8*(k-1)+(1:8),:) = drdA;
                            %                                 drcelldD{i}{ncross}(8*(k-1)+(1:8),:) = drdD;
                            %                                 drcelldt{i}{ncross}(8*(k-1)+(1:8),:) = drdt;
                            %
                            %                                 if gustflag == 1
                            %                                     drcellddv{i}{ncross}(8*(k-1)+(1:8),:) = drddv;
                            %                                 end
                            drvecdA{ntot}(8*(k-1)+(1:8),:) = drdA;
                            drvecdD{ntot}(8*(k-1)+(1:8),:) = drdD;
                            drvecdt{ntot}(8*(k-1)+(1:8),:) = drdt;
                            
                            if gustflag == 1
                                drvecddv{ntot}(8*(k-1)+(1:8),:) = drddv;
                            end
                        end
                        if morphflag == 1
                            if morph.camber
                                dN1dparam = zeros(6,ndparam);
                                dN2dparam = zeros(6,ndparam);
                                dr1dparam = zeros(4,ndparam);
                                dr2dparam = zeros(4,ndparam);
                                
                                for m=1:2
                                    % Loads on first end of the element
                                    if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                        devecelm1dparam = [0*strainraw{numcross}.dexoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                        devecelm2dparam = [0*strainraw{numcross}.dexoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    else
                                        devecelm1dparam = [strainraw{numcross}.dexoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                        devecelm2dparam = [strainraw{numcross}.dexoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    end
                                    
                                    
                                    [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1,ind1(m));evecelm2(1,ind1(m))],buckl_loc(i),1,2);
                                    [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2,ind1(m));evecelm2(2,ind1(m))],buckl_loc(i),1,2);
                                    [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3,ind1(m));evecelm2(3,ind1(m))],buckl_loc(i),1,2);
                                    
                                    
                                    devec1dparam(1,:) =  devec1devecelm1*[devecelm1dparam(1,:);devecelm2dparam(1,:)];
                                    devec1dparam(2,:) =  devec1devecelm2*[devecelm1dparam(2,:);devecelm2dparam(2,:)];
                                    devec1dparam(3,:) =  devec1devecelm3*[devecelm1dparam(3,:);devecelm2dparam(3,:)];
                                    
                                    
                                    dN1dparam(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1dparam;
                                    
                                    dr1dparam(2*(m-1)+(1:2),:)=dr1dparam(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dparam(3*(m-1)+(1:3),:);
                                    
                                    % Loads on second end of the element
                                    if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                        devecelm1dparam = [0*strainraw{numcross}.dexoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                        devecelm2dparam = [0*strainraw{numcross}.dexoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    else
                                        devecelm1dparam = [strainraw{numcross}.dexoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                        devecelm2dparam = [strainraw{numcross}.dexoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    end
                                    [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1,ind2(m));evecelm2(1,ind2(m))],buckl_loc(i+1),1,2);
                                    [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2,ind2(m));evecelm2(2,ind2(m))],buckl_loc(i+1),1,2);
                                    [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3,ind2(m));evecelm2(3,ind2(m))],buckl_loc(i+1),1,2);
                                    
                                    devec2dparam(1,:) =  devec2devecelm1*[devecelm1dparam(1,:);devecelm2dparam(1,:)];
                                    devec2dparam(2,:) =  devec2devecelm2*[devecelm1dparam(2,:);devecelm2dparam(2,:)];
                                    devec2dparam(3,:) =  devec2devecelm3*[devecelm1dparam(3,:);devecelm2dparam(3,:)];
                                    
                                    dN2dparam(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2dparam;
                                    
                                    dr2dparam(2*(m-1)+(1:2),:)=dr2dparam(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dparam(3*(m-1)+(1:3),:);
                                end
                                
                                drdparam = [dr1dparam;dr2dparam];
                                
                                
                                drvecdparam{ntot}(8*(k-1)+(1:8),:) = drdparam;
                            end
                            if morph.twist
                                dN1dphi = zeros(6,ndphi);
                                dN2dphi = zeros(6,ndphi);
                                dr1dphi = zeros(4,ndphi);
                                dr2dphi = zeros(4,ndphi);
                                
                                for m=1:2
                                    % Loads on first end of the element
                                    if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                        devecelm1dphi = [0*strainraw{numcross}.dexoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                        devecelm2dphi = [0*strainraw{numcross}.dexoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    else
                                        devecelm1dphi = [strainraw{numcross}.dexoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                        devecelm2dphi = [strainraw{numcross}.dexoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    end
                                    
                                    
                                    [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1,ind1(m));evecelm2(1,ind1(m))],buckl_loc(i),1,2);
                                    [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2,ind1(m));evecelm2(2,ind1(m))],buckl_loc(i),1,2);
                                    [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3,ind1(m));evecelm2(3,ind1(m))],buckl_loc(i),1,2);
                                    
                                    
                                    devec1dphi(1,:) =  devec1devecelm1*[devecelm1dphi(1,:);devecelm2dphi(1,:)];
                                    devec1dphi(2,:) =  devec1devecelm2*[devecelm1dphi(2,:);devecelm2dphi(2,:)];
                                    devec1dphi(3,:) =  devec1devecelm3*[devecelm1dphi(3,:);devecelm2dphi(3,:)];
                                    
                                    
                                    dN1dphi(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1dphi;
                                    
                                    dr1dphi(2*(m-1)+(1:2),:)=dr1dphi(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dphi(3*(m-1)+(1:3),:);
                                    
                                    % Loads on second end of the element
                                    if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                        devecelm1dphi = [0*strainraw{numcross}.dexoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                        devecelm2dphi = [0*strainraw{numcross}.dexoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    else
                                        devecelm1dphi = [strainraw{numcross}.dexoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                        devecelm2dphi = [strainraw{numcross}.dexoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    end
                                    [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1,ind2(m));evecelm2(1,ind2(m))],buckl_loc(i+1),1,2);
                                    [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2,ind2(m));evecelm2(2,ind2(m))],buckl_loc(i+1),1,2);
                                    [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3,ind2(m));evecelm2(3,ind2(m))],buckl_loc(i+1),1,2);
                                    
                                    devec2dphi(1,:) =  devec2devecelm1*[devecelm1dphi(1,:);devecelm2dphi(1,:)];
                                    devec2dphi(2,:) =  devec2devecelm2*[devecelm1dphi(2,:);devecelm2dphi(2,:)];
                                    devec2dphi(3,:) =  devec2devecelm3*[devecelm1dphi(3,:);devecelm2dphi(3,:)];
                                    
                                    dN2dphi(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2dphi;
                                    
                                    dr2dphi(2*(m-1)+(1:2),:)=dr2dphi(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dphi(3*(m-1)+(1:3),:);
                                end
                                
                                drdphi = [dr1dphi;dr2dphi];
                                
                                
                                drvecdphi{ntot}(8*(k-1)+(1:8),:) = drdphi;
                            end
                            if morph.fold
                                dN1dtheta = zeros(6,ndtheta);
                                dN2dtheta = zeros(6,ndtheta);
                                dr1dtheta = zeros(4,ndtheta);
                                dr2dtheta = zeros(4,ndtheta);
                                
                                for m=1:2
                                    % Loads on first end of the element
                                    if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                        devecelm1dtheta = [0*strainraw{numcross}.dexoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                        devecelm2dtheta = [0*strainraw{numcross}.dexoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    else
                                        devecelm1dtheta = [strainraw{numcross}.dexoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                        devecelm2dtheta = [strainraw{numcross}.dexoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    end
                                    
                                    
                                    [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1,ind1(m));evecelm2(1,ind1(m))],buckl_loc(i),1,2);
                                    [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2,ind1(m));evecelm2(2,ind1(m))],buckl_loc(i),1,2);
                                    [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3,ind1(m));evecelm2(3,ind1(m))],buckl_loc(i),1,2);
                                    
                                    
                                    devec1dtheta(1,:) =  devec1devecelm1*[devecelm1dtheta(1,:);devecelm2dtheta(1,:)];
                                    devec1dtheta(2,:) =  devec1devecelm2*[devecelm1dtheta(2,:);devecelm2dtheta(2,:)];
                                    devec1dtheta(3,:) =  devec1devecelm3*[devecelm1dtheta(3,:);devecelm2dtheta(3,:)];
                                    
                                    
                                    dN1dtheta(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1dtheta;
                                    
                                    dr1dtheta(2*(m-1)+(1:2),:)=dr1dtheta(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dtheta(3*(m-1)+(1:3),:);
                                    
                                    % Loads on second end of the element
                                    if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                        devecelm1dtheta = [0*strainraw{numcross}.dexoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                        devecelm2dtheta = [0*strainraw{numcross}.dexoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    else
                                        devecelm1dtheta = [strainraw{numcross}.dexoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                        devecelm2dtheta = [strainraw{numcross}.dexoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    end
                                    [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1,ind2(m));evecelm2(1,ind2(m))],buckl_loc(i+1),1,2);
                                    [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2,ind2(m));evecelm2(2,ind2(m))],buckl_loc(i+1),1,2);
                                    [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3,ind2(m));evecelm2(3,ind2(m))],buckl_loc(i+1),1,2);
                                    
                                    devec2dtheta(1,:) =  devec2devecelm1*[devecelm1dtheta(1,:);devecelm2dtheta(1,:)];
                                    devec2dtheta(2,:) =  devec2devecelm2*[devecelm1dtheta(2,:);devecelm2dtheta(2,:)];
                                    devec2dtheta(3,:) =  devec2devecelm3*[devecelm1dtheta(3,:);devecelm2dtheta(3,:)];
                                    
                                    dN2dtheta(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2dtheta;
                                    
                                    dr2dtheta(2*(m-1)+(1:2),:)=dr2dtheta(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dtheta(3*(m-1)+(1:3),:);
                                end
                                
                                drdtheta = [dr1dtheta;dr2dtheta];
                                
                                
                                drvecdtheta{ntot}(8*(k-1)+(1:8),:) = drdtheta;
                            end
                            if morph.shear
                                dN1dpsi = zeros(6,ndpsi);
                                dN2dpsi = zeros(6,ndpsi);
                                dr1dpsi = zeros(4,ndpsi);
                                dr2dpsi = zeros(4,ndpsi);
                                
                                for m=1:2
                                    % Loads on first end of the element
                                    if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                        devecelm1dpsi = [0*strainraw{numcross}.dexoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                        devecelm2dpsi = [0*strainraw{numcross}.dexoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    else
                                        devecelm1dpsi = [strainraw{numcross}.dexoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                        devecelm2dpsi = [strainraw{numcross}.dexoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    end
                                    
                                    
                                    [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1,ind1(m));evecelm2(1,ind1(m))],buckl_loc(i),1,2);
                                    [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2,ind1(m));evecelm2(2,ind1(m))],buckl_loc(i),1,2);
                                    [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3,ind1(m));evecelm2(3,ind1(m))],buckl_loc(i),1,2);
                                    
                                    
                                    devec1dpsi(1,:) =  devec1devecelm1*[devecelm1dpsi(1,:);devecelm2dpsi(1,:)];
                                    devec1dpsi(2,:) =  devec1devecelm2*[devecelm1dpsi(2,:);devecelm2dpsi(2,:)];
                                    devec1dpsi(3,:) =  devec1devecelm3*[devecelm1dpsi(3,:);devecelm2dpsi(3,:)];
                                    
                                    
                                    dN1dpsi(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1dpsi;
                                    
                                    dr1dpsi(2*(m-1)+(1:2),:)=dr1dpsi(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dpsi(3*(m-1)+(1:3),:);
                                    
                                    % Loads on second end of the element
                                    if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                        devecelm1dpsi = [0*strainraw{numcross}.dexoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                        devecelm2dpsi = [0*strainraw{numcross}.dexoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    else
                                        devecelm1dpsi = [strainraw{numcross}.dexoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                        devecelm2dpsi = [strainraw{numcross}.dexoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    end
                                    [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1,ind2(m));evecelm2(1,ind2(m))],buckl_loc(i+1),1,2);
                                    [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2,ind2(m));evecelm2(2,ind2(m))],buckl_loc(i+1),1,2);
                                    [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3,ind2(m));evecelm2(3,ind2(m))],buckl_loc(i+1),1,2);
                                    
                                    devec2dpsi(1,:) =  devec2devecelm1*[devecelm1dpsi(1,:);devecelm2dpsi(1,:)];
                                    devec2dpsi(2,:) =  devec2devecelm2*[devecelm1dpsi(2,:);devecelm2dpsi(2,:)];
                                    devec2dpsi(3,:) =  devec2devecelm3*[devecelm1dpsi(3,:);devecelm2dpsi(3,:)];
                                    
                                    dN2dpsi(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2dpsi;
                                    
                                    dr2dpsi(2*(m-1)+(1:2),:)=dr2dpsi(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dpsi(3*(m-1)+(1:3),:);
                                end
                                
                                drdpsi = [dr1dpsi;dr2dpsi];
                                
                                
                                drvecdpsi{ntot}(8*(k-1)+(1:8),:) = drdpsi;
                            end
                            if morph.span
                                dN1dext = zeros(6,ndext);
                                dN2dext = zeros(6,ndext);
                                dr1dext = zeros(4,ndext);
                                dr2dext = zeros(4,ndext);
                                
                                for m=1:2
                                    % Loads on first end of the element
                                    if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                        devecelm1dext = [0*strainraw{numcross}.dexoutdext{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdext{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdext{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                        devecelm2dext = [0*strainraw{numcross}.dexoutdext{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);0*strainraw{numcross}.deyoutdext{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdext{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    else
                                        devecelm1dext = [strainraw{numcross}.dexoutdext{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdext{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdext{1}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                        devecelm2dext = [strainraw{numcross}.dexoutdext{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.deyoutdext{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:);strainraw{numcross}.dgammaoutdext{2}(constant.buckl.crosselm{i}{ncross}{k}(ind1(m)),:)];
                                    end
                                    
                                    
                                    [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1,ind1(m));evecelm2(1,ind1(m))],buckl_loc(i),1,2);
                                    [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2,ind1(m));evecelm2(2,ind1(m))],buckl_loc(i),1,2);
                                    [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3,ind1(m));evecelm2(3,ind1(m))],buckl_loc(i),1,2);
                                    
                                    
                                    devec1dext(1,:) =  devec1devecelm1*[devecelm1dext(1,:);devecelm2dext(1,:)];
                                    devec1dext(2,:) =  devec1devecelm2*[devecelm1dext(2,:);devecelm2dext(2,:)];
                                    devec1dext(3,:) =  devec1devecelm3*[devecelm1dext(3,:);devecelm2dext(3,:)];
                                    
                                    
                                    dN1dext(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1dext;
                                    
                                    dr1dext(2*(m-1)+(1:2),:)=dr1dext(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dext(3*(m-1)+(1:3),:);
                                    
                                    % Loads on second end of the element
                                    if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                        devecelm1dext = [0*strainraw{numcross}.dexoutdext{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdext{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdext{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                        devecelm2dext = [0*strainraw{numcross}.dexoutdext{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);0*strainraw{numcross}.deyoutdext{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdext{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    else
                                        devecelm1dext = [strainraw{numcross}.dexoutdext{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdext{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdext{1}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                        devecelm2dext = [strainraw{numcross}.dexoutdext{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.deyoutdext{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:);strainraw{numcross}.dgammaoutdext{2}(constant.buckl.crosselm{i}{ncross}{k}(ind2(m)),:)];
                                    end
                                    [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1,ind2(m));evecelm2(1,ind2(m))],buckl_loc(i+1),1,2);
                                    [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2,ind2(m));evecelm2(2,ind2(m))],buckl_loc(i+1),1,2);
                                    [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3,ind2(m));evecelm2(3,ind2(m))],buckl_loc(i+1),1,2);
                                    
                                    devec2dext(1,:) =  devec2devecelm1*[devecelm1dext(1,:);devecelm2dext(1,:)];
                                    devec2dext(2,:) =  devec2devecelm2*[devecelm1dext(2,:);devecelm2dext(2,:)];
                                    devec2dext(3,:) =  devec2devecelm3*[devecelm1dext(3,:);devecelm2dext(3,:)];
                                    
                                    dN2dext(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2dext;
                                    
                                    dr2dext(2*(m-1)+(1:2),:)=dr2dext(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dext(3*(m-1)+(1:3),:);
                                end
                                
                                drdext = [dr1dext;dr2dext];
                                
                                
                                drvecdext{ntot}(8*(k-1)+(1:8),:) = drdext;
                            end
                        end
                    else
                        if tailflag == 1
                            dN1dA = zeros(3,ndA);
                            dN1dD = zeros(3,ndD);
                            dN1dt = zeros(3,ndt);
                            
                            dN2dA = zeros(3,ndA);
                            dN2dD = zeros(3,ndD);
                            dN2dt = zeros(3,ndt);
                            
                            dr1dA = zeros(2,ndA);
                            dr1dD = zeros(2,ndD);
                            dr1dt = zeros(2,ndt);
                            
                            dr2dA = zeros(2,ndA);
                            dr2dD = zeros(2,ndD);
                            dr2dt = zeros(2,ndt);
                            
                            if gustflag == 1
                                dN1ddv = zeros(3,nddv);
                                dN2ddv = zeros(3,nddv);
                                dr1ddv = zeros(2,nddv);
                                dr2ddv = zeros(2,nddv);
                            end
                            
                            dABDdA = zeros(9,ndA);
                            
                            % Loads on first end of the element
                            if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                devecelm1dA = [0*strainraw{numcross}.dexoutdA{1}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdA{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdA{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                devecelm2dA = [0*strainraw{numcross}.dexoutdA{2}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdA{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdA{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                
                                devecelm1dD = [0*strainraw{numcross}.dexoutdD{1}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdD{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdD{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                devecelm2dD = [0*strainraw{numcross}.dexoutdD{2}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdD{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdD{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                
                                devecelm1dt = [0*strainraw{numcross}.dexoutdt{1}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdt{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdt{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                devecelm2dt = [0*strainraw{numcross}.dexoutdt{2}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdt{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdt{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                
                                if gustflag == 1
                                    devecelm1ddv = [0*strainraw{numcross}.dexoutddv{1}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutddv{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutddv{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                    devecelm2ddv = [0*strainraw{numcross}.dexoutddv{2}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutddv{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutddv{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                end
                            else
                                devecelm1dA = [strainraw{numcross}.dexoutdA{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdA{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdA{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                devecelm2dA = [strainraw{numcross}.dexoutdA{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdA{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdA{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                
                                devecelm1dD = [strainraw{numcross}.dexoutdD{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdD{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdD{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                devecelm2dD = [strainraw{numcross}.dexoutdD{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdD{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdD{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                
                                devecelm1dt = [strainraw{numcross}.dexoutdt{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdt{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdt{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                devecelm2dt = [strainraw{numcross}.dexoutdt{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdt{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdt{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                
                                if gustflag == 1
                                    devecelm1ddv = [strainraw{numcross}.dexoutddv{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutddv{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutddv{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                    devecelm2ddv = [strainraw{numcross}.dexoutddv{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutddv{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutddv{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                end
                            end
                            [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i),1,2);
                            [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i),1,2);
                            [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i),1,2);
                            
                            devec1dA(1,:) =  devec1devecelm1*[devecelm1dA(1,:);devecelm2dA(1,:)];
                            devec1dA(2,:) =  devec1devecelm2*[devecelm1dA(2,:);devecelm2dA(2,:)];
                            devec1dA(3,:) =  devec1devecelm3*[devecelm1dA(3,:);devecelm2dA(3,:)];
                            
                            devec1dD(1,:) =  devec1devecelm1*[devecelm1dD(1,:);devecelm2dD(1,:)];
                            devec1dD(2,:) =  devec1devecelm2*[devecelm1dD(2,:);devecelm2dD(2,:)];
                            devec1dD(3,:) =  devec1devecelm3*[devecelm1dD(3,:);devecelm2dD(3,:)];
                            
                            devec1dt(1,:) =  devec1devecelm1*[devecelm1dt(1,:);devecelm2dt(1,:)];
                            devec1dt(2,:) =  devec1devecelm2*[devecelm1dt(2,:);devecelm2dt(2,:)];
                            devec1dt(3,:) =  devec1devecelm3*[devecelm1dt(3,:);devecelm2dt(3,:)];
                            
                            if gustflag == 1
                                devec1ddv(1,:) =  devec1devecelm1*[devecelm1ddv(1,:);devecelm2ddv(1,:)];
                                devec1ddv(2,:) =  devec1devecelm2*[devecelm1ddv(2,:);devecelm2ddv(2,:)];
                                devec1ddv(3,:) =  devec1devecelm3*[devecelm1ddv(3,:);devecelm2ddv(3,:)];
                            end
                            
                            row    = [1,1,1,2,2,3];     % Row number of independent A and D matrices elements
                            column = [1,2,3,2,3,3];     % Column number of independent A and D matrices elements
                            
                            for l = 1:6
                                
                                dABDdAloc = zeros(3,3);
                                dABDdAloc(row(l),column(l)) = 1;
                                dABDdAloc(column(l),row(l)) = 1;
                                
                                dABDdA(:,(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = reshape(dABDdAloc',[],1);
                                
                                dN1dA(:,(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = dABDdAloc*evec1;
                                
                                if row(l)==column(l)
                                    dr1dD(:,(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = [dr1dDpart{1}(row(l),column(l),1);dr1dDpart{1}(row(l),column(l),2)];
                                else
                                    dr1dD(:,(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = [dr1dDpart{1}(row(l),column(l),1);dr1dDpart{1}(row(l),column(l),2)]+[dr1dDpart{1}(column(l),row(l),1);dr1dDpart{1}(column(l),row(l),2)];
                                end
                            end
                            
                            dN1dA=dN1dA+ABD(1:3,1:3)*devec1dA;
                            dN1dD=ABD(1:3,1:3)*devec1dD;
                            dN1dt=ABD(1:3,1:3)*devec1dt;
                            
                            if gustflag == 1
                                dN1ddv=ABD(1:3,1:3)*devec1ddv;
                            end
                            
                            dr1dA=dr1dA+dr1dN{1}'*dN1dA;
                            dr1dD=dr1dD+dr1dN{1}'*dN1dD;
                            dr1dt=dr1dt+dr1dN{1}'*dN1dt;
                            
                            if gustflag == 1
                                dr1ddv=dr1ddv+dr1dN{1}'*dN1ddv;
                            end
                            
                            % Loads on second end of the element
                            [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i+1),1,2);
                            [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i+1),1,2);
                            [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i+1),1,2);
                            
                            devec2dA(1,:) =  devec2devecelm1*[devecelm1dA(1,:);devecelm2dA(1,:)];
                            devec2dA(2,:) =  devec2devecelm2*[devecelm1dA(2,:);devecelm2dA(2,:)];
                            devec2dA(3,:) =  devec2devecelm3*[devecelm1dA(3,:);devecelm2dA(3,:)];
                            
                            devec2dD(1,:) =  devec2devecelm1*[devecelm1dD(1,:);devecelm2dD(1,:)];
                            devec2dD(2,:) =  devec2devecelm2*[devecelm1dD(2,:);devecelm2dD(2,:)];
                            devec2dD(3,:) =  devec2devecelm3*[devecelm1dD(3,:);devecelm2dD(3,:)];
                            
                            devec2dt(1,:) =  devec2devecelm1*[devecelm1dt(1,:);devecelm2dt(1,:)];
                            devec2dt(2,:) =  devec2devecelm2*[devecelm1dt(2,:);devecelm2dt(2,:)];
                            devec2dt(3,:) =  devec2devecelm3*[devecelm1dt(3,:);devecelm2dt(3,:)];
                            
                            if gustflag == 1
                                devec2ddv(1,:) =  devec2devecelm1*[devecelm1ddv(1,:);devecelm2ddv(1,:)];
                                devec2ddv(2,:) =  devec2devecelm2*[devecelm1ddv(2,:);devecelm2ddv(2,:)];
                                devec2ddv(3,:) =  devec2devecelm3*[devecelm1ddv(3,:);devecelm2ddv(3,:)];
                            end
                            
                            row    = [1,1,1,2,2,3];     % Row number of independent A and D matrices elements
                            column = [1,2,3,2,3,3];     % Column number of independent A and D matrices elements
                            
                            for l = 1:6
                                
                                dABDdAloc = zeros(3,3);
                                dABDdAloc(row(l),column(l)) = 1;
                                dABDdAloc(column(l),row(l)) = 1;
                                
                                dABDdA(:,(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = reshape(dABDdAloc',[],1);
                                
                                dN2dA(:,(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = dABDdAloc*evec2;
                                
                                if row(l)==column(l)
                                    dr2dD(:,(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = [dr2dDpart{1}(row(l),column(l),1);dr2dDpart{1}(row(l),column(l),2)];
                                else
                                    dr2dD(:,(6*(constant.buckl.lam{i}{ncross}(k,1)-1)+l)) = [dr2dDpart{1}(row(l),column(l),1);dr2dDpart{1}(row(l),column(l),2)]+[dr2dDpart{1}(column(l),row(l),1);dr2dDpart{1}(column(l),row(l),2)];
                                end
                            end
                            
                            dN2dA=dN2dA+ABD(1:3,1:3)*devec2dA;
                            dN2dD=ABD(1:3,1:3)*devec2dD;
                            dN2dt=ABD(1:3,1:3)*devec2dt;
                            
                            if gustflag == 1
                                dN2ddv=ABD(1:3,1:3)*devec2ddv;
                            end
                            
                            dr2dA=dr2dN{1}'*dN2dA;
                            dr2dD=dr2dD+dr2dN{1}'*dN2dD;
                            dr2dt=dr2dN{1}'*dN2dt;
                            
                            if gustflag == 1
                                dr2ddv=dr2dN{1}'*dN2ddv;
                            end
                            
                            drdA = [dr1dA;dr1dA;dr2dA;dr2dA];
                            drdD = [dr1dD;dr1dD;dr2dD;dr2dD];
                            drdt = [dr1dt;dr1dt;dr2dt;dr2dt];
                            
                            if gustflag == 1
                                drddv = [dr1ddv;dr1ddv;dr2ddv;dr2ddv];
                            end
                            
                            drcelldA{i}{ncross}(8*(k-1)+(1:8),:) = drdA;
                            drcelldD{i}{ncross}(8*(k-1)+(1:8),:) = drdD;
                            drcelldt{i}{ncross}(8*(k-1)+(1:8),:) = drdt;
                            
                            if gustflag == 1
                                drcellddv{i}{ncross}(8*(k-1)+(1:8),:) = drddv;
                            end
                            
                            drvecdA{ntot}(8*(k-1)+(1:8),:) = drdA;
                            drvecdD{ntot}(8*(k-1)+(1:8),:) = drdD;
                            drvecdt{ntot}(8*(k-1)+(1:8),:) = drdt;
                            
                            if gustflag == 1
                                drvecddv{ntot}(8*(k-1)+(1:8),:) = drddv;
                            end
                        end
                        if morphflag == 1
                            if morph.camber
                                dN1dparam = zeros(3,ndparam);
                                dN2dparam = zeros(3,ndparam);
                                dr1dparam = zeros(2,ndparam);
                                dr2dparam = zeros(2,ndparam);
                                
                                
                                % Loads on first end of the element
                                if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                    devecelm1dparam = [0*strainraw{numcross}.dexoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                    devecelm2dparam = [0*strainraw{numcross}.dexoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                else
                                    devecelm1dparam = [strainraw{numcross}.dexoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                    devecelm2dparam = [strainraw{numcross}.dexoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                end
                                [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i),1,2);
                                [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i),1,2);
                                [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i),1,2);
                                
                                
                                devec1dparam(1,:) =  devec1devecelm1*[devecelm1dparam(1,:);devecelm2dparam(1,:)];
                                devec1dparam(2,:) =  devec1devecelm2*[devecelm1dparam(2,:);devecelm2dparam(2,:)];
                                devec1dparam(3,:) =  devec1devecelm3*[devecelm1dparam(3,:);devecelm2dparam(3,:)];
                                
                                dN1dparam=ABD(1:3,1:3)*devec1dparam;
                                
                                dr1dparam=dr1dparam+dr1dN{1}'*dN1dparam;
                                
                                % Loads on second end of the element
                                [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i+1),1,2);
                                [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i+1),1,2);
                                [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i+1),1,2);
                                
                                
                                devec2dparam(1,:) =  devec2devecelm1*[devecelm1dparam(1,:);devecelm2dparam(1,:)];
                                devec2dparam(2,:) =  devec2devecelm2*[devecelm1dparam(2,:);devecelm2dparam(2,:)];
                                devec2dparam(3,:) =  devec2devecelm3*[devecelm1dparam(3,:);devecelm2dparam(3,:)];
                                
                                dN2dparam=ABD(1:3,1:3)*devec2dparam;
                                
                                dr2dparam=dr2dN{1}'*dN2dparam;
                                
                                drdparam = [dr1dparam;dr1dparam;dr2dparam;dr2dparam];
                                
                                drcelldparam{i}{ncross}(8*(k-1)+(1:8),:) = drdparam;
                                
                                drvecdparam{ntot}(8*(k-1)+(1:8),:) = drdparam;
                            end
                            if morph.twist
                                dN1dphi = zeros(3,ndphi);
                                dN2dphi = zeros(3,ndphi);
                                dr1dphi = zeros(2,ndphi);
                                dr2dphi = zeros(2,ndphi);
                                
                                
                                % Loads on first end of the element
                                if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                    devecelm1dphi = [0*strainraw{numcross}.dexoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                    devecelm2dphi = [0*strainraw{numcross}.dexoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                else
                                    devecelm1dphi = [strainraw{numcross}.dexoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                    devecelm2dphi = [strainraw{numcross}.dexoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                end
                                [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i),1,2);
                                [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i),1,2);
                                [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i),1,2);
                                
                                
                                devec1dphi(1,:) =  devec1devecelm1*[devecelm1dphi(1,:);devecelm2dphi(1,:)];
                                devec1dphi(2,:) =  devec1devecelm2*[devecelm1dphi(2,:);devecelm2dphi(2,:)];
                                devec1dphi(3,:) =  devec1devecelm3*[devecelm1dphi(3,:);devecelm2dphi(3,:)];
                                
                                dN1dphi=ABD(1:3,1:3)*devec1dphi;
                                
                                dr1dphi=dr1dphi+dr1dN{1}'*dN1dphi;
                                
                                % Loads on second end of the element
                                [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i+1),1,2);
                                [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i+1),1,2);
                                [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i+1),1,2);
                                
                                
                                devec2dphi(1,:) =  devec2devecelm1*[devecelm1dphi(1,:);devecelm2dphi(1,:)];
                                devec2dphi(2,:) =  devec2devecelm2*[devecelm1dphi(2,:);devecelm2dphi(2,:)];
                                devec2dphi(3,:) =  devec2devecelm3*[devecelm1dphi(3,:);devecelm2dphi(3,:)];
                                
                                dN2dphi=ABD(1:3,1:3)*devec2dphi;
                                
                                dr2dphi=dr2dN{1}'*dN2dphi;
                                
                                drdphi = [dr1dphi;dr1dphi;dr2dphi;dr2dphi];
                                
                                drcelldphi{i}{ncross}(8*(k-1)+(1:8),:) = drdphi;
                                
                                drvecdphi{ntot}(8*(k-1)+(1:8),:) = drdphi;
                            end
                            if morph.fold
                                dN1dtheta = zeros(3,ndtheta);
                                dN2dtheta = zeros(3,ndtheta);
                                dr1dtheta = zeros(2,ndtheta);
                                dr2dtheta = zeros(2,ndtheta);
                                
                                
                                % Loads on first end of the element
                                if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                    devecelm1dtheta = [0*strainraw{numcross}.dexoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                    devecelm2dtheta = [0*strainraw{numcross}.dexoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                else
                                    devecelm1dtheta = [strainraw{numcross}.dexoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                    devecelm2dtheta = [strainraw{numcross}.dexoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                end
                                [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i),1,2);
                                [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i),1,2);
                                [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i),1,2);
                                
                                
                                devec1dtheta(1,:) =  devec1devecelm1*[devecelm1dtheta(1,:);devecelm2dtheta(1,:)];
                                devec1dtheta(2,:) =  devec1devecelm2*[devecelm1dtheta(2,:);devecelm2dtheta(2,:)];
                                devec1dtheta(3,:) =  devec1devecelm3*[devecelm1dtheta(3,:);devecelm2dtheta(3,:)];
                                
                                dN1dtheta=ABD(1:3,1:3)*devec1dtheta;
                                
                                dr1dtheta=dr1dtheta+dr1dN{1}'*dN1dtheta;
                                
                                % Loads on second end of the element
                                [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i+1),1,2);
                                [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i+1),1,2);
                                [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i+1),1,2);
                                
                                
                                devec2dtheta(1,:) =  devec2devecelm1*[devecelm1dtheta(1,:);devecelm2dtheta(1,:)];
                                devec2dtheta(2,:) =  devec2devecelm2*[devecelm1dtheta(2,:);devecelm2dtheta(2,:)];
                                devec2dtheta(3,:) =  devec2devecelm3*[devecelm1dtheta(3,:);devecelm2dtheta(3,:)];
                                
                                dN2dtheta=ABD(1:3,1:3)*devec2dtheta;
                                
                                dr2dtheta=dr2dN{1}'*dN2dtheta;
                                
                                drdtheta = [dr1dtheta;dr1dtheta;dr2dtheta;dr2dtheta];
                                
                                drcelldtheta{i}{ncross}(8*(k-1)+(1:8),:) = drdtheta;
                                
                                drvecdtheta{ntot}(8*(k-1)+(1:8),:) = drdtheta;
                            end
                            if morph.shear
                                dN1dpsi = zeros(3,ndpsi);
                                dN2dpsi = zeros(3,ndpsi);
                                dr1dpsi = zeros(2,ndpsi);
                                dr2dpsi = zeros(2,ndpsi);
                                
                                
                                % Loads on first end of the element
                                if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                    devecelm1dpsi = [0*strainraw{numcross}.dexoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                    devecelm2dpsi = [0*strainraw{numcross}.dexoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                else
                                    devecelm1dpsi = [strainraw{numcross}.dexoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                    devecelm2dpsi = [strainraw{numcross}.dexoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                end
                                [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i),1,2);
                                [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i),1,2);
                                [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i),1,2);
                                
                                
                                devec1dpsi(1,:) =  devec1devecelm1*[devecelm1dpsi(1,:);devecelm2dpsi(1,:)];
                                devec1dpsi(2,:) =  devec1devecelm2*[devecelm1dpsi(2,:);devecelm2dpsi(2,:)];
                                devec1dpsi(3,:) =  devec1devecelm3*[devecelm1dpsi(3,:);devecelm2dpsi(3,:)];
                                
                                dN1dpsi=ABD(1:3,1:3)*devec1dpsi;
                                
                                dr1dpsi=dr1dpsi+dr1dN{1}'*dN1dpsi;
                                
                                % Loads on second end of the element
                                [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i+1),1,2);
                                [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i+1),1,2);
                                [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i+1),1,2);
                                
                                
                                devec2dpsi(1,:) =  devec2devecelm1*[devecelm1dpsi(1,:);devecelm2dpsi(1,:)];
                                devec2dpsi(2,:) =  devec2devecelm2*[devecelm1dpsi(2,:);devecelm2dpsi(2,:)];
                                devec2dpsi(3,:) =  devec2devecelm3*[devecelm1dpsi(3,:);devecelm2dpsi(3,:)];
                                
                                dN2dpsi=ABD(1:3,1:3)*devec2dpsi;
                                
                                dr2dpsi=dr2dN{1}'*dN2dpsi;
                                
                                drdpsi = [dr1dpsi;dr1dpsi;dr2dpsi;dr2dpsi];
                                
                                drcelldpsi{i}{ncross}(8*(k-1)+(1:8),:) = drdpsi;
                                
                                drvecdpsi{ntot}(8*(k-1)+(1:8),:) = drdpsi;
                            end
                            if morph.span
                                dN1dext = zeros(3,ndext);
                                dN2dext = zeros(3,ndext);
                                dr1dext = zeros(2,ndext);
                                dr2dext = zeros(2,ndext);
                                
                                
                                % Loads on first end of the element
                                if constant.buckl.type{i}{ncross}(k) == 3 % Spars, only account for shear load
                                    devecelm1dext = [0*strainraw{numcross}.dexoutdext{1}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdext{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdext{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                    devecelm2dext = [0*strainraw{numcross}.dexoutdext{2}(constant.buckl.crosselm{i}{ncross}{k},:);0*strainraw{numcross}.deyoutdext{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdext{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                else
                                    devecelm1dext = [strainraw{numcross}.dexoutdext{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdext{1}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdext{1}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                    devecelm2dext = [strainraw{numcross}.dexoutdext{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.deyoutdext{2}(constant.buckl.crosselm{i}{ncross}{k},:);strainraw{numcross}.dgammaoutdext{2}(constant.buckl.crosselm{i}{ncross}{k},:)];
                                end
                                [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i),1,2);
                                [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i),1,2);
                                [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i),1,2);
                                
                                
                                devec1dext(1,:) =  devec1devecelm1*[devecelm1dext(1,:);devecelm2dext(1,:)];
                                devec1dext(2,:) =  devec1devecelm2*[devecelm1dext(2,:);devecelm2dext(2,:)];
                                devec1dext(3,:) =  devec1devecelm3*[devecelm1dext(3,:);devecelm2dext(3,:)];
                                
                                dN1dext=ABD(1:3,1:3)*devec1dext;
                                
                                dr1dext=dr1dext+dr1dN{1}'*dN1dext;
                                
                                % Loads on second end of the element
                                [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i+1),1,2);
                                [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i+1),1,2);
                                [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec1-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i+1),1,2);
                                
                                
                                devec2dext(1,:) =  devec2devecelm1*[devecelm1dext(1,:);devecelm2dext(1,:)];
                                devec2dext(2,:) =  devec2devecelm2*[devecelm1dext(2,:);devecelm2dext(2,:)];
                                devec2dext(3,:) =  devec2devecelm3*[devecelm1dext(3,:);devecelm2dext(3,:)];
                                
                                dN2dext=ABD(1:3,1:3)*devec2dext;
                                
                                dr2dext=dr2dN{1}'*dN2dext;
                                
                                drdext = [dr1dext;dr1dext;dr2dext;dr2dext];
                                
                                drcelldext{i}{ncross}(8*(k-1)+(1:8),:) = drdext;
                                
                                drvecdext{ntot}(8*(k-1)+(1:8),:) = drdext;
                            end
                        end
                    end
                end
            end
            panelm(i,ncross)=constant.buckl.Npan{i}(1,ncross);
        end
    else % Buckling panel covers multiple elements
        Npan = 0;
        ntotstart = ntot;
        for j=nsec1:nsec2
            nsec = j;
            for ncross=1:size(constant.buckl.yzpanel{i},2)
                if constant.buckl.Npan{i}(j-nsec1+1,ncross) ~= 0
                    numcross = Numcross{nsec}(ncross);
                end
                if j == nsec1 || length(Npan)<ncross
                    Npan(ncross) = 0;
                end
                ntot = ntotstart+ncross;
                for k=1:constant.buckl.Npan{i}(j-nsec1+1,ncross)
                    try
                        ABD = [lampar.A(3*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+(1:3),:) zeros(3); zeros(3) lampar.D(3*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+(1:3),:)];
                    catch err
                        keyboard
                    end
                    
                    if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                        evecelm1 = [0*strainraw{numcross}.exout(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},1)';0*strainraw{numcross}.eyout(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},1)';strainraw{numcross}.gammaout(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},1)'];
                        evecelm2 = [0*strainraw{numcross}.exout(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},2)';0*strainraw{numcross}.eyout(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},2)';strainraw{numcross}.gammaout(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},2)'];
                    else
                        evecelm1 = [strainraw{numcross}.exout(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},1)';strainraw{numcross}.eyout(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},1)';strainraw{numcross}.gammaout(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},1)'];
                        evecelm2 = [strainraw{numcross}.exout(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},2)';strainraw{numcross}.eyout(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},2)';strainraw{numcross}.gammaout(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},2)'];
                    end
                    
                    clear evec1 evec2
                    if buckl_loc(i)<=constant.str.xyz((nsec-1)*3+2)
                        evec1 = evecelm1;
%                         evec1 = evecelm2;
                    else
                        evec1(1,:) = interp1(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1,:);evecelm2(1,:)],buckl_loc(i));
                        evec1(2,:) = interp1(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2,:);evecelm2(2,:)],buckl_loc(i));
                        evec1(3,:) = interp1(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3,:);evecelm2(3,:)],buckl_loc(i));
                    end
                    if buckl_loc(i+1)>=constant.str.xyz((nsec-1)*3+5)
                        evec2 = evecelm2;
%                         evec2 = evecelm1;
                    else
                        evec2(1,:) = interp1(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1,:);evecelm2(1,:)],buckl_loc(i+1));
                        evec2(2,:) = interp1(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2,:);evecelm2(2,:)],buckl_loc(i+1));
                        evec2(3,:) = interp1(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3,:);evecelm2(3,:)],buckl_loc(i+1));
                    end
                    
                    N1 = ABD(1:3,1:3)*evec1;
                    N2 = ABD(1:3,1:3)*evec2;
                    
                    clear r1 r2 dr1dDpart dr2dDpart dr1dN dr2dN
                    for l=1:size(N1,2)
                        [r1(:,l),dr1dDpart{l},dr1dN{l}] = quadbuckl(N1(:,l),ABD(4:6,4:6),2,constant.buckl.stiff{i}{ncross}{k+Npan(ncross)});
                        [r2(:,l),dr2dDpart{l},dr2dN{l}] = quadbuckl(N2(:,l),ABD(4:6,4:6),2,constant.buckl.stiff{i}{ncross}{k+Npan(ncross)});
                    end
                    
                    % Identify critical modes
                    if size(r1,2)>1
                        [~,ind1] = sort(r1(1,:),'descend');
                        [~,ind2] = sort(r2(1,:),'descend');
                        
                        r = [r1(:,ind1(1));r1(:,ind1(2));r2(:,ind2(1));r2(:,ind2(2))];
                        
                    else
                        r = [r1;r1;r2;r2];
                    end
                    
                    rcell{i}{ncross}(k+Npan(ncross),:) = r';
                    rvec{ntot}(k+Npan(ncross),:) = r';
                    if ders == 1
                        if size(r1,2)>1
                            if tailflag == 1
                                dN1dA = zeros(6,ndA);
                                dN1dD = zeros(6,ndD);
                                dN1dt = zeros(6,ndt);
                                
                                dN2dA = zeros(6,ndA);
                                dN2dD = zeros(6,ndD);
                                dN2dt = zeros(6,ndt);
                                
                                dr1dA = zeros(4,ndA);
                                dr1dD = zeros(4,ndD);
                                dr1dt = zeros(4,ndt);
                                
                                dr2dA = zeros(4,ndA);
                                dr2dD = zeros(4,ndD);
                                dr2dt = zeros(4,ndt);
                                
                                if gustflag == 1
                                    dN1ddv = zeros(6,nddv);
                                    dN2ddv = zeros(6,nddv);
                                    dr1ddv = zeros(4,nddv);
                                    dr2ddv = zeros(4,nddv);
                                end
                                
                                dABDdA = zeros(9,ndA);
                                for m=1:2
                                    % Loads on first end of the element
                                    if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                        devecelm1dA = [0*strainraw{numcross}.dexoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        devecelm2dA = [0*strainraw{numcross}.dexoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        
                                        devecelm1dD = [0*strainraw{numcross}.dexoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        devecelm2dD = [0*strainraw{numcross}.dexoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        
                                        devecelm1dt = [0*strainraw{numcross}.dexoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        devecelm2dt = [0*strainraw{numcross}.dexoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        
                                        if gustflag == 1
                                            devecelm1ddv = [0*strainraw{numcross}.dexoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                            devecelm2ddv = [0*strainraw{numcross}.dexoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        end
                                    else
                                        devecelm1dA = [strainraw{numcross}.dexoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        devecelm2dA = [strainraw{numcross}.dexoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        
                                        devecelm1dD = [strainraw{numcross}.dexoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        devecelm2dD = [strainraw{numcross}.dexoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        
                                        devecelm1dt = [strainraw{numcross}.dexoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        devecelm2dt = [strainraw{numcross}.dexoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        
                                        if gustflag == 1
                                            devecelm1ddv = [strainraw{numcross}.dexoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                            devecelm2ddv = [strainraw{numcross}.dexoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        end
                                    end
                                    
                                    if buckl_loc(i)<=constant.str.xyz((nsec-1)*3+2)
                                        devec1dA = devecelm1dA;
                                        devec1dD = devecelm1dD;
                                        devec1dt = devecelm1dt;
                                        
                                        if gustflag == 1
                                            devec1ddv = devecelm1ddv;
                                        end
                                    else
                                        [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1,ind1(m));evecelm2(1,ind1(m))],buckl_loc(i),1,2);
                                        [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2,ind1(m));evecelm2(2,ind1(m))],buckl_loc(i),1,2);
                                        [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3,ind1(m));evecelm2(3,ind1(m))],buckl_loc(i),1,2);
                                        
                                        devec1dA(1,:) =  devec1devecelm1*[devecelm1dA(1,:);devecelm2dA(1,:)];
                                        devec1dA(2,:) =  devec1devecelm2*[devecelm1dA(2,:);devecelm2dA(2,:)];
                                        devec1dA(3,:) =  devec1devecelm3*[devecelm1dA(3,:);devecelm2dA(3,:)];
                                        
                                        devec1dD(1,:) =  devec1devecelm1*[devecelm1dD(1,:);devecelm2dD(1,:)];
                                        devec1dD(2,:) =  devec1devecelm2*[devecelm1dD(2,:);devecelm2dD(2,:)];
                                        devec1dD(3,:) =  devec1devecelm3*[devecelm1dD(3,:);devecelm2dD(3,:)];
                                        
                                        devec1dt(1,:) =  devec1devecelm1*[devecelm1dt(1,:);devecelm2dt(1,:)];
                                        devec1dt(2,:) =  devec1devecelm2*[devecelm1dt(2,:);devecelm2dt(2,:)];
                                        devec1dt(3,:) =  devec1devecelm3*[devecelm1dt(3,:);devecelm2dt(3,:)];
                                        
                                        if gustflag == 1
                                            devec1ddv(1,:) =  devec1devecelm1*[devecelm1ddv(1,:);devecelm2ddv(1,:)];
                                            devec1ddv(2,:) =  devec1devecelm2*[devecelm1ddv(2,:);devecelm2ddv(2,:)];
                                            devec1ddv(3,:) =  devec1devecelm3*[devecelm1ddv(3,:);devecelm2ddv(3,:)];
                                        end
                                    end
                                    
                                    
                                    row    = [1,1,1,2,2,3];     % Row number of independent A and D matrices elements
                                    column = [1,2,3,2,3,3];     % Column number of independent A and D matrices elements
                                    
                                    for l = 1:6
                                        
                                        dABDdAloc = zeros(3,3);
                                        dABDdAloc(row(l),column(l)) = 1;
                                        dABDdAloc(column(l),row(l)) = 1;
                                        
                                        dABDdA(:,(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = reshape(dABDdAloc',[],1);
                                        
                                        dN1dA(3*(m-1)+(1:3),(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = dABDdAloc*evec1(:,ind1(m));
                                        
                                        if row(l)==column(l)
                                            dr1dD(2*(m-1)+(1:2),(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = [dr1dDpart{ind1(m)}(row(l),column(l),1);dr1dDpart{ind1(m)}(row(l),column(l),2)];
                                        else
                                            dr1dD(2*(m-1)+(1:2),(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = [dr1dDpart{ind1(m)}(row(l),column(l),1);dr1dDpart{ind1(m)}(row(l),column(l),2)]+[dr1dDpart{ind1(m)}(column(l),row(l),1);dr1dDpart{ind1(m)}(column(l),row(l),2)];
                                        end
                                    end
                                    
                                    dN1dA(3*(m-1)+(1:3),:)=dN1dA(3*(m-1)+(1:3),:)+ABD(1:3,1:3)*devec1dA;
                                    dN1dD(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1dD;
                                    dN1dt(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1dt;
                                    
                                    if gustflag == 1
                                        dN1ddv(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1ddv;
                                    end
                                    
                                    dr1dA(2*(m-1)+(1:2),:)=dr1dA(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dA(3*(m-1)+(1:3),:);
                                    dr1dD(2*(m-1)+(1:2),:)=dr1dD(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dD(3*(m-1)+(1:3),:);
                                    dr1dt(2*(m-1)+(1:2),:)=dr1dt(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dt(3*(m-1)+(1:3),:);
                                    
                                    if gustflag == 1
                                        dr1ddv(2*(m-1)+(1:2),:)=dr1ddv(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1ddv(3*(m-1)+(1:3),:);
                                    end
                                    
                                    % Loads on second end of the element
                                    if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                        devecelm1dA = [0*strainraw{numcross}.dexoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        devecelm2dA = [0*strainraw{numcross}.dexoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        
                                        devecelm1dD = [0*strainraw{numcross}.dexoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        devecelm2dD = [0*strainraw{numcross}.dexoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        
                                        devecelm1dt = [0*strainraw{numcross}.dexoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        devecelm2dt = [0*strainraw{numcross}.dexoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        
                                        if gustflag == 1
                                            devecelm1ddv = [0*strainraw{numcross}.dexoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                            devecelm2ddv = [0*strainraw{numcross}.dexoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        end
                                    else
                                        devecelm1dA = [strainraw{numcross}.dexoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        devecelm2dA = [strainraw{numcross}.dexoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        
                                        devecelm1dD = [strainraw{numcross}.dexoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        devecelm2dD = [strainraw{numcross}.dexoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        
                                        devecelm1dt = [strainraw{numcross}.dexoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        devecelm2dt = [strainraw{numcross}.dexoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        
                                        if gustflag == 1
                                            devecelm1ddv = [strainraw{numcross}.dexoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                            devecelm2ddv = [strainraw{numcross}.dexoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        end
                                    end
                                    
                                    if buckl_loc(i+1)>=constant.str.xyz((nsec-1)*3+5)
                                        devec2dA = devecelm2dA;
                                        devec2dD = devecelm2dD;
                                        devec2dt = devecelm2dt;
                                        
                                        if gustflag == 1
                                            devec2ddv = devecelm2ddv;
                                        end
                                    else
                                        [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1,ind2(m));evecelm2(1,ind2(m))],buckl_loc(i+1),1,2);
                                        [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2,ind2(m));evecelm2(2,ind2(m))],buckl_loc(i+1),1,2);
                                        [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3,ind2(m));evecelm2(3,ind2(m))],buckl_loc(i+1),1,2);
                                        
                                        devec2dA(1,:) =  devec2devecelm1*[devecelm1dA(1,:);devecelm2dA(1,:)];
                                        devec2dA(2,:) =  devec2devecelm2*[devecelm1dA(2,:);devecelm2dA(2,:)];
                                        devec2dA(3,:) =  devec2devecelm3*[devecelm1dA(3,:);devecelm2dA(3,:)];
                                        
                                        devec2dD(1,:) =  devec2devecelm1*[devecelm1dD(1,:);devecelm2dD(1,:)];
                                        devec2dD(2,:) =  devec2devecelm2*[devecelm1dD(2,:);devecelm2dD(2,:)];
                                        devec2dD(3,:) =  devec2devecelm3*[devecelm1dD(3,:);devecelm2dD(3,:)];
                                        
                                        devec2dt(1,:) =  devec2devecelm1*[devecelm1dt(1,:);devecelm2dt(1,:)];
                                        devec2dt(2,:) =  devec2devecelm2*[devecelm1dt(2,:);devecelm2dt(2,:)];
                                        devec2dt(3,:) =  devec2devecelm3*[devecelm1dt(3,:);devecelm2dt(3,:)];
                                        
                                        if gustflag == 1
                                            devec2ddv(1,:) =  devec2devecelm1*[devecelm1ddv(1,:);devecelm2ddv(1,:)];
                                            devec2ddv(2,:) =  devec2devecelm2*[devecelm1ddv(2,:);devecelm2ddv(2,:)];
                                            devec2ddv(3,:) =  devec2devecelm3*[devecelm1ddv(3,:);devecelm2ddv(3,:)];
                                        end
                                    end
                                    row    = [1,1,1,2,2,3];     % Row number of independent A and D matrices elements
                                    column = [1,2,3,2,3,3];     % Column number of independent A and D matrices elements
                                    
                                    for l = 1:6
                                        
                                        dABDdAloc = zeros(3,3);
                                        dABDdAloc(row(l),column(l)) = 1;
                                        dABDdAloc(column(l),row(l)) = 1;
                                        
                                        dABDdA(:,(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = reshape(dABDdAloc',[],1);
                                        
                                        dN2dA(3*(m-1)+(1:3),(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = dABDdAloc*evec2(:,ind2(m));
                                        
                                        if row(l)==column(l)
                                            dr2dD(2*(m-1)+(1:2),(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = [dr2dDpart{ind2(m)}(row(l),column(l),1);dr2dDpart{ind2(m)}(row(l),column(l),2)];
                                        else
                                            dr2dD(2*(m-1)+(1:2),(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = [dr2dDpart{ind2(m)}(row(l),column(l),1);dr2dDpart{ind2(m)}(row(l),column(l),2)]+[dr2dDpart{ind2(m)}(column(l),row(l),1);dr2dDpart{ind2(m)}(column(l),row(l),2)];
                                        end
                                    end
                                    
                                    dN2dA(3*(m-1)+(1:3),:)=dN2dA(3*(m-1)+(1:3),:)+ABD(1:3,1:3)*devec2dA;
                                    dN2dD(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2dD;
                                    dN2dt(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2dt;
                                    
                                    if gustflag == 1
                                        dN2ddv(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2ddv;
                                    end
                                    
                                    dr2dA(2*(m-1)+(1:2),:)=dr2dA(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dA(3*(m-1)+(1:3),:);
                                    dr2dD(2*(m-1)+(1:2),:)=dr2dD(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dD(3*(m-1)+(1:3),:);
                                    dr2dt(2*(m-1)+(1:2),:)=dr2dt(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dt(3*(m-1)+(1:3),:);
                                    
                                    if gustflag == 1
                                        dr2ddv(2*(m-1)+(1:2),:)=dr2ddv(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2ddv(3*(m-1)+(1:3),:);
                                    end
                                end
                                
                                drdA = [dr1dA;dr2dA];
                                drdD = [dr1dD;dr2dD];
                                drdt = [dr1dt;dr2dt];
                                
                                if gustflag == 1
                                    drddv = [dr1ddv;dr2ddv];
                                end
                                
                                drcelldA{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdA;
                                drcelldD{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdD;
                                drcelldt{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdt;
                                
                                if gustflag == 1
                                    drcellddv{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drddv;
                                end
                                
                                drvecdA{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdA;
                                drvecdD{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdD;
                                drvecdt{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdt;
                                
                                if gustflag == 1
                                    drvecddv{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drddv;
                                end
                            end
                            if morphflag == 1
                                if morph.camber
                                    dN1dparam = zeros(6,ndparam);
                                    dN2dparam = zeros(6,ndparam);
                                    dr1dparam = zeros(4,ndparam);
                                    dr2dparam = zeros(4,ndparam);
                                    
                                    for m=1:2
                                        % Loads on first end of the element
                                        if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                            devecelm1dparam = [0*strainraw{numcross}.dexoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                            devecelm2dparam = [0*strainraw{numcross}.dexoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        else
                                            devecelm1dparam = [strainraw{numcross}.dexoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                            devecelm2dparam = [strainraw{numcross}.dexoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        end
                                        
                                        if buckl_loc(i)<=constant.str.xyz((nsec-1)*3+2)
                                            devec1dparam = devecelm1dparam;
                                        else
                                            [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1,ind1(m));evecelm2(1,ind1(m))],buckl_loc(i),1,2);
                                            [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2,ind1(m));evecelm2(2,ind1(m))],buckl_loc(i),1,2);
                                            [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3,ind1(m));evecelm2(3,ind1(m))],buckl_loc(i),1,2);
                                            
                                            devec1dparam(1,:) =  devec1devecelm1*[devecelm1dparam(1,:);devecelm2dparam(1,:)];
                                            devec1dparam(2,:) =  devec1devecelm2*[devecelm1dparam(2,:);devecelm2dparam(2,:)];
                                            devec1dparam(3,:) =  devec1devecelm3*[devecelm1dparam(3,:);devecelm2dparam(3,:)];
                                        end
                                        
                                        dN1dparam(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1dparam;
                                        
                                        dr1dparam(2*(m-1)+(1:2),:)=dr1dparam(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dparam(3*(m-1)+(1:3),:);
                                        
                                        
                                        % Loads on second end of the element
                                        if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                            
                                            devecelm1dparam = [0*strainraw{numcross}.dexoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                            devecelm2dparam = [0*strainraw{numcross}.dexoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        else
                                            
                                            devecelm1dparam = [strainraw{numcross}.dexoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                            devecelm2dparam = [strainraw{numcross}.dexoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        end
                                        
                                        if buckl_loc(i+1)>=constant.str.xyz((nsec-1)*3+5)
                                            
                                            devec2dparam = devecelm2dparam;
                                        else
                                            [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1,ind2(m));evecelm2(1,ind2(m))],buckl_loc(i+1),1,2);
                                            [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2,ind2(m));evecelm2(2,ind2(m))],buckl_loc(i+1),1,2);
                                            [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3,ind2(m));evecelm2(3,ind2(m))],buckl_loc(i+1),1,2);
                                            
                                            devec2dparam(1,:) =  devec2devecelm1*[devecelm1dparam(1,:);devecelm2dparam(1,:)];
                                            devec2dparam(2,:) =  devec2devecelm2*[devecelm1dparam(2,:);devecelm2dparam(2,:)];
                                            devec2dparam(3,:) =  devec2devecelm3*[devecelm1dparam(3,:);devecelm2dparam(3,:)];
                                        end
                                        
                                        dN2dparam(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2dparam;
                                        
                                        dr2dparam(2*(m-1)+(1:2),:)=dr2dparam(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dparam(3*(m-1)+(1:3),:);
                                    end
                                    
                                    drdparam = [dr1dparam;dr2dparam];
                                    
                                    drcelldparam{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdparam;
                                    
                                    drvecdparam{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdparam;
                                end
                                if morph.twist
                                    dN1dphi = zeros(6,ndphi);
                                    dN2dphi = zeros(6,ndphi);
                                    dr1dphi = zeros(4,ndphi);
                                    dr2dphi = zeros(4,ndphi);
                                    
                                    for m=1:2
                                        % Loads on first end of the element
                                        if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                            devecelm1dphi = [0*strainraw{numcross}.dexoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                            devecelm2dphi = [0*strainraw{numcross}.dexoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        else
                                            devecelm1dphi = [strainraw{numcross}.dexoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                            devecelm2dphi = [strainraw{numcross}.dexoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        end
                                        
                                        if buckl_loc(i)<=constant.str.xyz((nsec-1)*3+2)
                                            devec1dphi = devecelm1dphi;
                                        else
                                            [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1,ind1(m));evecelm2(1,ind1(m))],buckl_loc(i),1,2);
                                            [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2,ind1(m));evecelm2(2,ind1(m))],buckl_loc(i),1,2);
                                            [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3,ind1(m));evecelm2(3,ind1(m))],buckl_loc(i),1,2);
                                            
                                            devec1dphi(1,:) =  devec1devecelm1*[devecelm1dphi(1,:);devecelm2dphi(1,:)];
                                            devec1dphi(2,:) =  devec1devecelm2*[devecelm1dphi(2,:);devecelm2dphi(2,:)];
                                            devec1dphi(3,:) =  devec1devecelm3*[devecelm1dphi(3,:);devecelm2dphi(3,:)];
                                        end
                                        
                                        dN1dphi(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1dphi;
                                        
                                        dr1dphi(2*(m-1)+(1:2),:)=dr1dphi(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dphi(3*(m-1)+(1:3),:);
                                        
                                        
                                        % Loads on second end of the element
                                        if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                            
                                            devecelm1dphi = [0*strainraw{numcross}.dexoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                            devecelm2dphi = [0*strainraw{numcross}.dexoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        else
                                            
                                            devecelm1dphi = [strainraw{numcross}.dexoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                            devecelm2dphi = [strainraw{numcross}.dexoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        end
                                        
                                        if buckl_loc(i+1)>=constant.str.xyz((nsec-1)*3+5)
                                            
                                            devec2dphi = devecelm2dphi;
                                        else
                                            [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1,ind2(m));evecelm2(1,ind2(m))],buckl_loc(i+1),1,2);
                                            [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2,ind2(m));evecelm2(2,ind2(m))],buckl_loc(i+1),1,2);
                                            [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3,ind2(m));evecelm2(3,ind2(m))],buckl_loc(i+1),1,2);
                                            
                                            devec2dphi(1,:) =  devec2devecelm1*[devecelm1dphi(1,:);devecelm2dphi(1,:)];
                                            devec2dphi(2,:) =  devec2devecelm2*[devecelm1dphi(2,:);devecelm2dphi(2,:)];
                                            devec2dphi(3,:) =  devec2devecelm3*[devecelm1dphi(3,:);devecelm2dphi(3,:)];
                                        end
                                        
                                        dN2dphi(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2dphi;
                                        
                                        dr2dphi(2*(m-1)+(1:2),:)=dr2dphi(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dphi(3*(m-1)+(1:3),:);
                                    end
                                    
                                    drdphi = [dr1dphi;dr2dphi];
                                    
                                    drcelldphi{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdphi;
                                    
                                    drvecdphi{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdphi;
                                end
                                if morph.fold
                                    dN1dtheta = zeros(6,ndtheta);
                                    dN2dtheta = zeros(6,ndtheta);
                                    dr1dtheta = zeros(4,ndtheta);
                                    dr2dtheta = zeros(4,ndtheta);
                                    
                                    for m=1:2
                                        % Loads on first end of the element
                                        if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                            devecelm1dtheta = [0*strainraw{numcross}.dexoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                            devecelm2dtheta = [0*strainraw{numcross}.dexoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        else
                                            devecelm1dtheta = [strainraw{numcross}.dexoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                            devecelm2dtheta = [strainraw{numcross}.dexoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        end
                                        
                                        if buckl_loc(i)<=constant.str.xyz((nsec-1)*3+2)
                                            devec1dtheta = devecelm1dtheta;
                                        else
                                            [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1,ind1(m));evecelm2(1,ind1(m))],buckl_loc(i),1,2);
                                            [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2,ind1(m));evecelm2(2,ind1(m))],buckl_loc(i),1,2);
                                            [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3,ind1(m));evecelm2(3,ind1(m))],buckl_loc(i),1,2);
                                            
                                            devec1dtheta(1,:) =  devec1devecelm1*[devecelm1dtheta(1,:);devecelm2dtheta(1,:)];
                                            devec1dtheta(2,:) =  devec1devecelm2*[devecelm1dtheta(2,:);devecelm2dtheta(2,:)];
                                            devec1dtheta(3,:) =  devec1devecelm3*[devecelm1dtheta(3,:);devecelm2dtheta(3,:)];
                                        end
                                        
                                        dN1dtheta(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1dtheta;
                                        
                                        dr1dtheta(2*(m-1)+(1:2),:)=dr1dtheta(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dtheta(3*(m-1)+(1:3),:);
                                        
                                        
                                        % Loads on second end of the element
                                        if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                            
                                            devecelm1dtheta = [0*strainraw{numcross}.dexoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                            devecelm2dtheta = [0*strainraw{numcross}.dexoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        else
                                            
                                            devecelm1dtheta = [strainraw{numcross}.dexoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                            devecelm2dtheta = [strainraw{numcross}.dexoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        end
                                        
                                        if buckl_loc(i+1)>=constant.str.xyz((nsec-1)*3+5)
                                            
                                            devec2dtheta = devecelm2dtheta;
                                        else
                                            [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1,ind2(m));evecelm2(1,ind2(m))],buckl_loc(i+1),1,2);
                                            [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2,ind2(m));evecelm2(2,ind2(m))],buckl_loc(i+1),1,2);
                                            [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3,ind2(m));evecelm2(3,ind2(m))],buckl_loc(i+1),1,2);
                                            
                                            devec2dtheta(1,:) =  devec2devecelm1*[devecelm1dtheta(1,:);devecelm2dtheta(1,:)];
                                            devec2dtheta(2,:) =  devec2devecelm2*[devecelm1dtheta(2,:);devecelm2dtheta(2,:)];
                                            devec2dtheta(3,:) =  devec2devecelm3*[devecelm1dtheta(3,:);devecelm2dtheta(3,:)];
                                        end
                                        
                                        dN2dtheta(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2dtheta;
                                        
                                        dr2dtheta(2*(m-1)+(1:2),:)=dr2dtheta(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dtheta(3*(m-1)+(1:3),:);
                                    end
                                    
                                    drdtheta = [dr1dtheta;dr2dtheta];
                                    
                                    drcelldtheta{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdtheta;
                                    
                                    drvecdtheta{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdtheta;
                                end
                                if morph.shear
                                    dN1dpsi = zeros(6,ndpsi);
                                    dN2dpsi = zeros(6,ndpsi);
                                    dr1dpsi = zeros(4,ndpsi);
                                    dr2dpsi = zeros(4,ndpsi);
                                    
                                    for m=1:2
                                        % Loads on first end of the element
                                        if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                            devecelm1dpsi = [0*strainraw{numcross}.dexoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                            devecelm2dpsi = [0*strainraw{numcross}.dexoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        else
                                            devecelm1dpsi = [strainraw{numcross}.dexoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                            devecelm2dpsi = [strainraw{numcross}.dexoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        end
                                        
                                        if buckl_loc(i)<=constant.str.xyz((nsec-1)*3+2)
                                            devec1dpsi = devecelm1dpsi;
                                        else
                                            [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1,ind1(m));evecelm2(1,ind1(m))],buckl_loc(i),1,2);
                                            [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2,ind1(m));evecelm2(2,ind1(m))],buckl_loc(i),1,2);
                                            [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3,ind1(m));evecelm2(3,ind1(m))],buckl_loc(i),1,2);
                                            
                                            devec1dpsi(1,:) =  devec1devecelm1*[devecelm1dpsi(1,:);devecelm2dpsi(1,:)];
                                            devec1dpsi(2,:) =  devec1devecelm2*[devecelm1dpsi(2,:);devecelm2dpsi(2,:)];
                                            devec1dpsi(3,:) =  devec1devecelm3*[devecelm1dpsi(3,:);devecelm2dpsi(3,:)];
                                        end
                                        
                                        dN1dpsi(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1dpsi;
                                        
                                        dr1dpsi(2*(m-1)+(1:2),:)=dr1dpsi(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dpsi(3*(m-1)+(1:3),:);
                                        
                                        
                                        % Loads on second end of the element
                                        if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                            
                                            devecelm1dpsi = [0*strainraw{numcross}.dexoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                            devecelm2dpsi = [0*strainraw{numcross}.dexoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        else
                                            
                                            devecelm1dpsi = [strainraw{numcross}.dexoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                            devecelm2dpsi = [strainraw{numcross}.dexoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        end
                                        
                                        if buckl_loc(i+1)>=constant.str.xyz((nsec-1)*3+5)
                                            
                                            devec2dpsi = devecelm2dpsi;
                                        else
                                            [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1,ind2(m));evecelm2(1,ind2(m))],buckl_loc(i+1),1,2);
                                            [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2,ind2(m));evecelm2(2,ind2(m))],buckl_loc(i+1),1,2);
                                            [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3,ind2(m));evecelm2(3,ind2(m))],buckl_loc(i+1),1,2);
                                            
                                            devec2dpsi(1,:) =  devec2devecelm1*[devecelm1dpsi(1,:);devecelm2dpsi(1,:)];
                                            devec2dpsi(2,:) =  devec2devecelm2*[devecelm1dpsi(2,:);devecelm2dpsi(2,:)];
                                            devec2dpsi(3,:) =  devec2devecelm3*[devecelm1dpsi(3,:);devecelm2dpsi(3,:)];
                                        end
                                        
                                        dN2dpsi(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2dpsi;
                                        
                                        dr2dpsi(2*(m-1)+(1:2),:)=dr2dpsi(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dpsi(3*(m-1)+(1:3),:);
                                    end
                                    
                                    drdpsi = [dr1dpsi;dr2dpsi];
                                    
                                    drcelldpsi{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdpsi;
                                    
                                    drvecdpsi{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdpsi;
                                end
                                if morph.span
                                    dN1dext = zeros(6,ndext);
                                    dN2dext = zeros(6,ndext);
                                    dr1dext = zeros(4,ndext);
                                    dr2dext = zeros(4,ndext);
                                    
                                    for m=1:2
                                        % Loads on first end of the element
                                        if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                            devecelm1dext = [0*strainraw{numcross}.dexoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                            devecelm2dext = [0*strainraw{numcross}.dexoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);0*strainraw{numcross}.deyoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        else
                                            devecelm1dext = [strainraw{numcross}.dexoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                            devecelm2dext = [strainraw{numcross}.dexoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.deyoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:);strainraw{numcross}.dgammaoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind1(m)),:)];
                                        end
                                        
                                        if buckl_loc(i)<=constant.str.xyz((nsec-1)*3+2)
                                            devec1dext = devecelm1dext;
                                        else
                                            [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1,ind1(m));evecelm2(1,ind1(m))],buckl_loc(i),1,2);
                                            [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2,ind1(m));evecelm2(2,ind1(m))],buckl_loc(i),1,2);
                                            [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3,ind1(m));evecelm2(3,ind1(m))],buckl_loc(i),1,2);
                                            
                                            devec1dext(1,:) =  devec1devecelm1*[devecelm1dext(1,:);devecelm2dext(1,:)];
                                            devec1dext(2,:) =  devec1devecelm2*[devecelm1dext(2,:);devecelm2dext(2,:)];
                                            devec1dext(3,:) =  devec1devecelm3*[devecelm1dext(3,:);devecelm2dext(3,:)];
                                        end
                                        
                                        dN1dext(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec1dext;
                                        
                                        dr1dext(2*(m-1)+(1:2),:)=dr1dext(2*(m-1)+(1:2),:)+dr1dN{ind1(m)}'*dN1dext(3*(m-1)+(1:3),:);
                                        
                                        
                                        % Loads on second end of the element
                                        if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                            
                                            devecelm1dext = [0*strainraw{numcross}.dexoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                            devecelm2dext = [0*strainraw{numcross}.dexoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);0*strainraw{numcross}.deyoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        else
                                            
                                            devecelm1dext = [strainraw{numcross}.dexoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                            devecelm2dext = [strainraw{numcross}.dexoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.deyoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:);strainraw{numcross}.dgammaoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)}(ind2(m)),:)];
                                        end
                                        
                                        if buckl_loc(i+1)>=constant.str.xyz((nsec-1)*3+5)
                                            
                                            devec2dext = devecelm2dext;
                                        else
                                            [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1,ind2(m));evecelm2(1,ind2(m))],buckl_loc(i+1),1,2);
                                            [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2,ind2(m));evecelm2(2,ind2(m))],buckl_loc(i+1),1,2);
                                            [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3,ind2(m));evecelm2(3,ind2(m))],buckl_loc(i+1),1,2);
                                            
                                            devec2dext(1,:) =  devec2devecelm1*[devecelm1dext(1,:);devecelm2dext(1,:)];
                                            devec2dext(2,:) =  devec2devecelm2*[devecelm1dext(2,:);devecelm2dext(2,:)];
                                            devec2dext(3,:) =  devec2devecelm3*[devecelm1dext(3,:);devecelm2dext(3,:)];
                                        end
                                        
                                        dN2dext(3*(m-1)+(1:3),:)=ABD(1:3,1:3)*devec2dext;
                                        
                                        dr2dext(2*(m-1)+(1:2),:)=dr2dext(2*(m-1)+(1:2),:)+dr2dN{ind2(m)}'*dN2dext(3*(m-1)+(1:3),:);
                                    end
                                    
                                    drdext = [dr1dext;dr2dext];
                                    
                                    drcelldext{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdext;
                                    
                                    drvecdext{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdext;
                                end
                            end
                        else
                            if tailflag == 1
                                dN1dA = zeros(3,ndA);
                                dN1dD = zeros(3,ndD);
                                dN1dt = zeros(3,ndt);
                                
                                dN2dA = zeros(3,ndA);
                                dN2dD = zeros(3,ndD);
                                dN2dt = zeros(3,ndt);
                                
                                dr1dA = zeros(2,ndA);
                                dr1dD = zeros(2,ndD);
                                dr1dt = zeros(2,ndt);
                                
                                dr2dA = zeros(2,ndA);
                                dr2dD = zeros(2,ndD);
                                dr2dt = zeros(2,ndt);
                                
                                if gustflag == 1
                                    dN1ddv = zeros(3,nddv);
                                    dN2ddv = zeros(3,nddv);
                                    dr1ddv = zeros(2,nddv);
                                    dr2ddv = zeros(2,nddv);
                                end
                                
                                dABDdA = zeros(9,ndA);
                                
                                % Loads on first end of the element
                                if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                    devecelm1dA = [0*strainraw{numcross}.dexoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    devecelm2dA = [0*strainraw{numcross}.dexoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    
                                    devecelm1dD = [0*strainraw{numcross}.dexoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    devecelm2dD = [0*strainraw{numcross}.dexoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    
                                    devecelm1dt = [0*strainraw{numcross}.dexoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    devecelm2dt = [0*strainraw{numcross}.dexoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    
                                    if gustflag == 1
                                        devecelm1ddv = [0*strainraw{numcross}.dexoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                        devecelm2ddv = [0*strainraw{numcross}.dexoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    end
                                else
                                    devecelm1dA = [strainraw{numcross}.dexoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdA{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    devecelm2dA = [strainraw{numcross}.dexoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdA{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    
                                    devecelm1dD = [strainraw{numcross}.dexoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdD{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    devecelm2dD = [strainraw{numcross}.dexoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdD{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    
                                    devecelm1dt = [strainraw{numcross}.dexoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdt{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    devecelm2dt = [strainraw{numcross}.dexoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdt{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    
                                    if gustflag == 1
                                        devecelm1ddv = [strainraw{numcross}.dexoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutddv{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                        devecelm2ddv = [strainraw{numcross}.dexoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutddv{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    end
                                end
                                
                                if buckl_loc(i)<=constant.str.xyz((nsec-1)*3+2)
                                    devec1dA = devecelm1dA;
                                    devec1dD = devecelm1dD;
                                    devec1dt = devecelm1dt;
                                    
                                    if gustflag == 1
                                        devec1ddv = devecelm1ddv;
                                    end
                                else
                                    [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i),1,2);
                                    [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i),1,2);
                                    [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i),1,2);
                                    
                                    devec1dA(1,:) =  devec1devecelm1*[devecelm1dA(1,:);devecelm2dA(1,:)];
                                    devec1dA(2,:) =  devec1devecelm2*[devecelm1dA(2,:);devecelm2dA(2,:)];
                                    devec1dA(3,:) =  devec1devecelm3*[devecelm1dA(3,:);devecelm2dA(3,:)];
                                    
                                    devec1dD(1,:) =  devec1devecelm1*[devecelm1dD(1,:);devecelm2dD(1,:)];
                                    devec1dD(2,:) =  devec1devecelm2*[devecelm1dD(2,:);devecelm2dD(2,:)];
                                    devec1dD(3,:) =  devec1devecelm3*[devecelm1dD(3,:);devecelm2dD(3,:)];
                                    
                                    devec1dt(1,:) =  devec1devecelm1*[devecelm1dt(1,:);devecelm2dt(1,:)];
                                    devec1dt(2,:) =  devec1devecelm2*[devecelm1dt(2,:);devecelm2dt(2,:)];
                                    devec1dt(3,:) =  devec1devecelm3*[devecelm1dt(3,:);devecelm2dt(3,:)];
                                    
                                    if gustflag == 1
                                        devec1ddv(1,:) =  devec1devecelm1*[devecelm1ddv(1,:);devecelm2ddv(1,:)];
                                        devec1ddv(2,:) =  devec1devecelm2*[devecelm1ddv(2,:);devecelm2ddv(2,:)];
                                        devec1ddv(3,:) =  devec1devecelm3*[devecelm1ddv(3,:);devecelm2ddv(3,:)];
                                    end
                                end
                                
                                row    = [1,1,1,2,2,3];     % Row number of independent A and D matrices elements
                                column = [1,2,3,2,3,3];     % Column number of independent A and D matrices elements
                                
                                for l = 1:6
                                    
                                    dABDdAloc = zeros(3,3);
                                    dABDdAloc(row(l),column(l)) = 1;
                                    dABDdAloc(column(l),row(l)) = 1;
                                    
                                    dABDdA(:,(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = reshape(dABDdAloc',[],1);
                                    
                                    dN1dA(:,(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = dABDdAloc*evec1;
                                    
                                    if row(l)==column(l)
                                        dr1dD(:,(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = [dr1dDpart{1}(row(l),column(l),1);dr1dDpart{1}(row(l),column(l),2)];
                                    else
                                        dr1dD(:,(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = [dr1dDpart{1}(row(l),column(l),1);dr1dDpart{1}(row(l),column(l),2)]+[dr1dDpart{1}(column(l),row(l),1);dr1dDpart{1}(column(l),row(l),2)];
                                    end
                                end
                                
                                dN1dA=dN1dA+ABD(1:3,1:3)*devec1dA;
                                dN1dD=ABD(1:3,1:3)*devec1dD;
                                dN1dt=ABD(1:3,1:3)*devec1dt;
                                
                                if gustflag == 1
                                    dN1ddv=ABD(1:3,1:3)*devec1ddv;
                                end
                                
                                dr1dA=dr1dA+dr1dN{1}'*dN1dA;
                                dr1dD=dr1dD+dr1dN{1}'*dN1dD;
                                dr1dt=dr1dt+dr1dN{1}'*dN1dt;
                                
                                if gustflag == 1
                                    dr1ddv=dr1ddv+dr1dN{1}'*dN1ddv;
                                end
                                
                                % Loads on second end of the element
                                if buckl_loc(i+1)>=constant.str.xyz((nsec-1)*3+5)
                                    devec2dA = devecelm2dA;
                                    devec2dD = devecelm2dD;
                                    devec2dt = devecelm2dt;
                                    
                                    if gustflag == 1
                                        devec2ddv = devecelm2ddv;
                                    end
                                else
                                    [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i+1),1,2);
                                    [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i+1),1,2);
                                    [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i+1),1,2);
                                    
                                    devec2dA(1,:) =  devec2devecelm1*[devecelm1dA(1,:);devecelm2dA(1,:)];
                                    devec2dA(2,:) =  devec2devecelm2*[devecelm1dA(2,:);devecelm2dA(2,:)];
                                    devec2dA(3,:) =  devec2devecelm3*[devecelm1dA(3,:);devecelm2dA(3,:)];
                                    
                                    devec2dD(1,:) =  devec2devecelm1*[devecelm1dD(1,:);devecelm2dD(1,:)];
                                    devec2dD(2,:) =  devec2devecelm2*[devecelm1dD(2,:);devecelm2dD(2,:)];
                                    devec2dD(3,:) =  devec2devecelm3*[devecelm1dD(3,:);devecelm2dD(3,:)];
                                    
                                    devec2dt(1,:) =  devec2devecelm1*[devecelm1dt(1,:);devecelm2dt(1,:)];
                                    devec2dt(2,:) =  devec2devecelm2*[devecelm1dt(2,:);devecelm2dt(2,:)];
                                    devec2dt(3,:) =  devec2devecelm3*[devecelm1dt(3,:);devecelm2dt(3,:)];
                                    
                                    if gustflag == 1
                                        devec2ddv(1,:) =  devec2devecelm1*[devecelm1ddv(1,:);devecelm2ddv(1,:)];
                                        devec2ddv(2,:) =  devec2devecelm2*[devecelm1ddv(2,:);devecelm2ddv(2,:)];
                                        devec2ddv(3,:) =  devec2devecelm3*[devecelm1ddv(3,:);devecelm2ddv(3,:)];
                                    end
                                end
                                
                                row    = [1,1,1,2,2,3];     % Row number of independent A and D matrices elements
                                column = [1,2,3,2,3,3];     % Column number of independent A and D matrices elements
                                
                                for l = 1:6
                                    
                                    dABDdAloc = zeros(3,3);
                                    dABDdAloc(row(l),column(l)) = 1;
                                    dABDdAloc(column(l),row(l)) = 1;
                                    
                                    dABDdA(:,(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = reshape(dABDdAloc',[],1);
                                    
                                    dN2dA(:,(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = dABDdAloc*evec2;
                                    
                                    if row(l)==column(l)
                                        dr2dD(:,(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = [dr2dDpart{1}(row(l),column(l),1);dr2dDpart{1}(row(l),column(l),2)];
                                    else
                                        dr2dD(:,(6*(constant.buckl.lam{i}{ncross}(k+Npan(ncross),1)-1)+l)) = [dr2dDpart{1}(row(l),column(l),1);dr2dDpart{1}(row(l),column(l),2)]+[dr2dDpart{1}(column(l),row(l),1);dr2dDpart{1}(column(l),row(l),2)];
                                    end
                                end
                                
                                dN2dA=dN2dA+ABD(1:3,1:3)*devec2dA;
                                dN2dD=ABD(1:3,1:3)*devec2dD;
                                dN2dt=ABD(1:3,1:3)*devec2dt;
                                
                                if gustflag == 1
                                    dN2ddv=ABD(1:3,1:3)*devec2ddv;
                                end
                                
                                dr2dA=dr2dN{1}'*dN2dA;
                                dr2dD=dr2dD+dr2dN{1}'*dN2dD;
                                dr2dt=dr2dN{1}'*dN2dt;
                                
                                if gustflag == 1
                                    dr2ddv=dr2dN{1}'*dN2ddv;
                                end
                                
                                drdA = [dr1dA;dr1dA;dr2dA;dr2dA];
                                drdD = [dr1dD;dr1dD;dr2dD;dr2dD];
                                drdt = [dr1dt;dr1dt;dr2dt;dr2dt];
                                
                                if gustflag == 1
                                    drddv = [dr1ddv;dr1ddv;dr2ddv;dr2ddv];
                                end
                                
                                drcelldA{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdA;
                                drcelldD{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdD;
                                drcelldt{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdt;
                                
                                if gustflag == 1
                                    drcellddv{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drddv;
                                end
                                
                                drvecdA{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdA;
                                drvecdD{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdD;
                                drvecdt{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdt;
                                
                                if gustflag == 1
                                    drvecddv{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drddv;
                                end
                            end
                            if morphflag == 1
                                if morph.camber
                                    dN1dparam = zeros(3,ndparam);
                                    dN2dparam = zeros(3,ndparam);
                                    dr1dparam = zeros(2,ndparam);
                                    dr2dparam = zeros(2,ndparam);
                                    
                                    % Loads on first end of the element
                                    if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                        devecelm1dparam = [0*strainraw{numcross}.dexoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                        devecelm2dparam = [0*strainraw{numcross}.dexoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    else
                                        devecelm1dparam = [strainraw{numcross}.dexoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdparam{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                        devecelm2dparam = [strainraw{numcross}.dexoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdparam{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    end
                                    
                                    if buckl_loc(i)<=constant.str.xyz((nsec-1)*3+2)
                                        devec1dparam = devecelm1dparam;
                                    else
                                        [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i),1,2);
                                        [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i),1,2);
                                        [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i),1,2);
                                        
                                        devec1dparam(1,:) =  devec1devecelm1*[devecelm1dparam(1,:);devecelm2dparam(1,:)];
                                        devec1dparam(2,:) =  devec1devecelm2*[devecelm1dparam(2,:);devecelm2dparam(2,:)];
                                        devec1dparam(3,:) =  devec1devecelm3*[devecelm1dparam(3,:);devecelm2dparam(3,:)];
                                    end
                                    
                                    dN1dparam=ABD(1:3,1:3)*devec1dparam;
                                    
                                    dr1dparam=dr1dparam+dr1dN{1}'*dN1dparam;
                                    
                                    % Loads on second end of the element
                                    if buckl_loc(i+1)>=constant.str.xyz((nsec-1)*3+5)
                                        devec2dparam = devecelm2dparam;
                                    else
                                        [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i+1),1,2);
                                        [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i+1),1,2);
                                        [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i+1),1,2);
                                        
                                        devec2dparam(1,:) =  devec2devecelm1*[devecelm1dparam(1,:);devecelm2dparam(1,:)];
                                        devec2dparam(2,:) =  devec2devecelm2*[devecelm1dparam(2,:);devecelm2dparam(2,:)];
                                        devec2dparam(3,:) =  devec2devecelm3*[devecelm1dparam(3,:);devecelm2dparam(3,:)];
                                    end
                                    
                                    dN2dparam=ABD(1:3,1:3)*devec2dparam;
                                    
                                    dr2dparam=dr2dN{1}'*dN2dparam;
                                    
                                    drdparam = [dr1dparam;dr1dparam;dr2dparam;dr2dparam];
                                    
                                    drcelldparam{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdparam;
                                    
                                    drvecdparam{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdparam;
                                end
                                if morph.twist
                                    dN1dphi = zeros(3,ndphi);
                                    dN2dphi = zeros(3,ndphi);
                                    dr1dphi = zeros(2,ndphi);
                                    dr2dphi = zeros(2,ndphi);
                                    
                                    % Loads on first end of the element
                                    if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                        devecelm1dphi = [0*strainraw{numcross}.dexoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                        devecelm2dphi = [0*strainraw{numcross}.dexoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    else
                                        devecelm1dphi = [strainraw{numcross}.dexoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdphi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                        devecelm2dphi = [strainraw{numcross}.dexoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdphi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    end
                                    
                                    if buckl_loc(i)<=constant.str.xyz((nsec-1)*3+2)
                                        devec1dphi = devecelm1dphi;
                                    else
                                        [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i),1,2);
                                        [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i),1,2);
                                        [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i),1,2);
                                        
                                        devec1dphi(1,:) =  devec1devecelm1*[devecelm1dphi(1,:);devecelm2dphi(1,:)];
                                        devec1dphi(2,:) =  devec1devecelm2*[devecelm1dphi(2,:);devecelm2dphi(2,:)];
                                        devec1dphi(3,:) =  devec1devecelm3*[devecelm1dphi(3,:);devecelm2dphi(3,:)];
                                    end
                                    
                                    dN1dphi=ABD(1:3,1:3)*devec1dphi;
                                    
                                    dr1dphi=dr1dphi+dr1dN{1}'*dN1dphi;
                                    
                                    % Loads on second end of the element
                                    if buckl_loc(i+1)>=constant.str.xyz((nsec-1)*3+5)
                                        devec2dphi = devecelm2dphi;
                                    else
                                        [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i+1),1,2);
                                        [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i+1),1,2);
                                        [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i+1),1,2);
                                        
                                        devec2dphi(1,:) =  devec2devecelm1*[devecelm1dphi(1,:);devecelm2dphi(1,:)];
                                        devec2dphi(2,:) =  devec2devecelm2*[devecelm1dphi(2,:);devecelm2dphi(2,:)];
                                        devec2dphi(3,:) =  devec2devecelm3*[devecelm1dphi(3,:);devecelm2dphi(3,:)];
                                    end
                                    
                                    dN2dphi=ABD(1:3,1:3)*devec2dphi;
                                    
                                    dr2dphi=dr2dN{1}'*dN2dphi;
                                    
                                    drdphi = [dr1dphi;dr1dphi;dr2dphi;dr2dphi];
                                    
                                    drcelldphi{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdphi;
                                    
                                    drvecdphi{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdphi;
                                end
                                if morph.fold
                                    dN1dtheta = zeros(3,ndtheta);
                                    dN2dtheta = zeros(3,ndtheta);
                                    dr1dtheta = zeros(2,ndtheta);
                                    dr2dtheta = zeros(2,ndtheta);
                                    
                                    % Loads on first end of the element
                                    if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                        devecelm1dtheta = [0*strainraw{numcross}.dexoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                        devecelm2dtheta = [0*strainraw{numcross}.dexoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    else
                                        devecelm1dtheta = [strainraw{numcross}.dexoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdtheta{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                        devecelm2dtheta = [strainraw{numcross}.dexoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdtheta{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    end
                                    
                                    if buckl_loc(i)<=constant.str.xyz((nsec-1)*3+2)
                                        devec1dtheta = devecelm1dtheta;
                                    else
                                        [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i),1,2);
                                        [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i),1,2);
                                        [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i),1,2);
                                        
                                        devec1dtheta(1,:) =  devec1devecelm1*[devecelm1dtheta(1,:);devecelm2dtheta(1,:)];
                                        devec1dtheta(2,:) =  devec1devecelm2*[devecelm1dtheta(2,:);devecelm2dtheta(2,:)];
                                        devec1dtheta(3,:) =  devec1devecelm3*[devecelm1dtheta(3,:);devecelm2dtheta(3,:)];
                                    end
                                    
                                    dN1dtheta=ABD(1:3,1:3)*devec1dtheta;
                                    
                                    dr1dtheta=dr1dtheta+dr1dN{1}'*dN1dtheta;
                                    
                                    % Loads on second end of the element
                                    if buckl_loc(i+1)>=constant.str.xyz((nsec-1)*3+5)
                                        devec2dtheta = devecelm2dtheta;
                                    else
                                        [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i+1),1,2);
                                        [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i+1),1,2);
                                        [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i+1),1,2);
                                        
                                        devec2dtheta(1,:) =  devec2devecelm1*[devecelm1dtheta(1,:);devecelm2dtheta(1,:)];
                                        devec2dtheta(2,:) =  devec2devecelm2*[devecelm1dtheta(2,:);devecelm2dtheta(2,:)];
                                        devec2dtheta(3,:) =  devec2devecelm3*[devecelm1dtheta(3,:);devecelm2dtheta(3,:)];
                                    end
                                    
                                    dN2dtheta=ABD(1:3,1:3)*devec2dtheta;
                                    
                                    dr2dtheta=dr2dN{1}'*dN2dtheta;
                                    
                                    drdtheta = [dr1dtheta;dr1dtheta;dr2dtheta;dr2dtheta];
                                    
                                    drcelldtheta{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdtheta;
                                    
                                    drvecdtheta{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdtheta;
                                end
                                if morph.shear
                                    dN1dpsi = zeros(3,ndpsi);
                                    dN2dpsi = zeros(3,ndpsi);
                                    dr1dpsi = zeros(2,ndpsi);
                                    dr2dpsi = zeros(2,ndpsi);
                                    
                                    % Loads on first end of the element
                                    if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                        devecelm1dpsi = [0*strainraw{numcross}.dexoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                        devecelm2dpsi = [0*strainraw{numcross}.dexoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    else
                                        devecelm1dpsi = [strainraw{numcross}.dexoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdpsi{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                        devecelm2dpsi = [strainraw{numcross}.dexoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdpsi{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    end
                                    
                                    if buckl_loc(i)<=constant.str.xyz((nsec-1)*3+2)
                                        devec1dpsi = devecelm1dpsi;
                                    else
                                        [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i),1,2);
                                        [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i),1,2);
                                        [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i),1,2);
                                        
                                        devec1dpsi(1,:) =  devec1devecelm1*[devecelm1dpsi(1,:);devecelm2dpsi(1,:)];
                                        devec1dpsi(2,:) =  devec1devecelm2*[devecelm1dpsi(2,:);devecelm2dpsi(2,:)];
                                        devec1dpsi(3,:) =  devec1devecelm3*[devecelm1dpsi(3,:);devecelm2dpsi(3,:)];
                                    end
                                    
                                    dN1dpsi=ABD(1:3,1:3)*devec1dpsi;
                                    
                                    dr1dpsi=dr1dpsi+dr1dN{1}'*dN1dpsi;
                                    
                                    % Loads on second end of the element
                                    if buckl_loc(i+1)>=constant.str.xyz((nsec-1)*3+5)
                                        devec2dpsi = devecelm2dpsi;
                                    else
                                        [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i+1),1,2);
                                        [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i+1),1,2);
                                        [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i+1),1,2);
                                        
                                        devec2dpsi(1,:) =  devec2devecelm1*[devecelm1dpsi(1,:);devecelm2dpsi(1,:)];
                                        devec2dpsi(2,:) =  devec2devecelm2*[devecelm1dpsi(2,:);devecelm2dpsi(2,:)];
                                        devec2dpsi(3,:) =  devec2devecelm3*[devecelm1dpsi(3,:);devecelm2dpsi(3,:)];
                                    end
                                    
                                    dN2dpsi=ABD(1:3,1:3)*devec2dpsi;
                                    
                                    dr2dpsi=dr2dN{1}'*dN2dpsi;
                                    
                                    drdpsi = [dr1dpsi;dr1dpsi;dr2dpsi;dr2dpsi];
                                    
                                    drcelldpsi{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdpsi;
                                    
                                    drvecdpsi{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdpsi;
                                end
                                if morph.span
                                    dN1dext = zeros(3,ndext);
                                    dN2dext = zeros(3,ndext);
                                    dr1dext = zeros(2,ndext);
                                    dr2dext = zeros(2,ndext);
                                    
                                    % Loads on first end of the element
                                    if constant.buckl.type{i}{ncross}(k+Npan(ncross)) == 3 % Spars, only account for shear load
                                        devecelm1dext = [0*strainraw{numcross}.dexoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                        devecelm2dext = [0*strainraw{numcross}.dexoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);0*strainraw{numcross}.deyoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    else
                                        devecelm1dext = [strainraw{numcross}.dexoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdext{1}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                        devecelm2dext = [strainraw{numcross}.dexoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.deyoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:);strainraw{numcross}.dgammaoutdext{2}(constant.buckl.crosselm{i}{ncross}{k+Npan(ncross)},:)];
                                    end
                                    
                                    if buckl_loc(i)<=constant.str.xyz((nsec-1)*3+2)
                                        devec1dext = devecelm1dext;
                                    else
                                        [~,devec1devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i),1,2);
                                        [~,devec1devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i),1,2);
                                        [~,devec1devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i),1,2);
                                        
                                        devec1dext(1,:) =  devec1devecelm1*[devecelm1dext(1,:);devecelm2dext(1,:)];
                                        devec1dext(2,:) =  devec1devecelm2*[devecelm1dext(2,:);devecelm2dext(2,:)];
                                        devec1dext(3,:) =  devec1devecelm3*[devecelm1dext(3,:);devecelm2dext(3,:)];
                                    end
                                    
                                    dN1dext=ABD(1:3,1:3)*devec1dext;
                                    
                                    dr1dext=dr1dext+dr1dN{1}'*dN1dext;
                                    
                                    % Loads on second end of the element
                                    if buckl_loc(i+1)>=constant.str.xyz((nsec-1)*3+5)
                                        devec2dext = devecelm2dext;
                                    else
                                        [~,devec2devecelm1,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(1);evecelm2(1)],buckl_loc(i+1),1,2);
                                        [~,devec2devecelm2,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(2);evecelm2(2)],buckl_loc(i+1),1,2);
                                        [~,devec2devecelm3,~] = dinterpolation(constant.str.xyz((nsec-1)*3+(2:3:5)),[evecelm1(3);evecelm2(3)],buckl_loc(i+1),1,2);
                                        
                                        devec2dext(1,:) =  devec2devecelm1*[devecelm1dext(1,:);devecelm2dext(1,:)];
                                        devec2dext(2,:) =  devec2devecelm2*[devecelm1dext(2,:);devecelm2dext(2,:)];
                                        devec2dext(3,:) =  devec2devecelm3*[devecelm1dext(3,:);devecelm2dext(3,:)];
                                    end
                                    
                                    dN2dext=ABD(1:3,1:3)*devec2dext;
                                    
                                    dr2dext=dr2dN{1}'*dN2dext;
                                    
                                    drdext = [dr1dext;dr1dext;dr2dext;dr2dext];
                                    
                                    drcelldext{i}{ncross}(8*(k+Npan(ncross)-1)+(1:8),:) = drdext;
                                    
                                    drvecdext{ntot}(8*(k+Npan(ncross)-1)+(1:8),:) = drdext;
                                end
                            end
                        end
                    end
                end
                
                Npan(ncross) = Npan(ncross) + constant.buckl.Npan{i}(j-nsec1+1,ncross);
                panelm(i,ncross) = Npan(ncross);
            end
        end
    end
end

buckl.rfull = rcell;

% Identify the critical buckling panels per laminate
rmat = cell2mat(rvec');
if ders == 1
    if tailflag == 1
        drmatdA = cell2mat(drvecdA');
        drmatdD = cell2mat(drvecdD');
        drmatdt = cell2mat(drvecdt');
        
        if gustflag == 1
            drmatddv = cell2mat(drvecddv');
        end
    end
    if morphflag == 1
        if morph.camber
            drmatdparam = cell2mat(drvecdparam');
        end
        if morph.twist
            drmatdphi = cell2mat(drvecdphi');
        end
        if morph.fold
            drmatdtheta = cell2mat(drvecdtheta');
        end
        if morph.shear
            drmatdpsi = cell2mat(drvecdpsi');
        end
        if morph.span
            drmatdext = cell2mat(drvecdext');
        end
    end
end

BPL = constant.opt.constraint.BucklPerLam;
for i=1:length(constant.lam.ID)
    
    rind = find(lammat==i);
    
    if length(rind)>=BPL
        rlam = rmat(rind,:);
        rlam_crit = max(rlam,[],2);
        
        [~,crit_panel] = sort(rlam_crit,'descend');
        
        crit_panel = sort(crit_panel(1:BPL),'ascend');
        
        buckl.r(8*BPL*(i-1)+(1:8*BPL),1) = reshape((rlam(crit_panel(1:BPL),:))',[],1);
        
        buckl.pan{i}(:,1) = rind(crit_panel);
        buckl.sec{i}(:,1) = bucklsec(rind(crit_panel));
        buckl.cross{i}(:,1) = bucklcross(rind(crit_panel));
        
        if ders == 1
            indvec = reshape((8*(rind(:,ones(8,1))-1)+repmat(1:8,length(rind),1))',[],1);
            crit_panelvec = reshape((8*(crit_panel(1:BPL,ones(8,1))-1)+repmat(1:8,BPL,1))',[],1);
            if tailflag == 1
                drlamdA = drmatdA(indvec,:);
                drlamdD = drmatdD(indvec,:);
                drlamdt = drmatdt(indvec,:);
                
                if gustflag == 1
                    drlamddv = drmatddv(indvec,:);
                end
                
                buckl.drdA(8*BPL*(i-1)+(1:8*BPL),:) = drlamdA(crit_panelvec,:);
                buckl.drdD(8*BPL*(i-1)+(1:8*BPL),:) = drlamdD(crit_panelvec,:);
                buckl.drdt(8*BPL*(i-1)+(1:8*BPL),:) = drlamdt(crit_panelvec,:);
                
                if gustflag == 1
                    buckl.drddv(8*BPL*(i-1)+(1:8*BPL),:) = drlamddv(crit_panelvec,:);
                end
            end
            if morphflag == 1
                if morph.camber
                    drlamdparam = drmatdparam(indvec,:);
                    buckl.drdparam(8*BPL*(i-1)+(1:8*BPL),:) = drlamdparam(crit_panelvec,:);
                end
                if morph.twist
                    drlamdphi = drmatdphi(indvec,:);
                    buckl.drdphi(8*BPL*(i-1)+(1:8*BPL),:) = drlamdphi(crit_panelvec,:);
                end
                if morph.fold
                    drlamdtheta = drmatdtheta(indvec,:);
                    buckl.drdtheta(8*BPL*(i-1)+(1:8*BPL),:) = drlamdtheta(crit_panelvec,:);
                end
                if morph.shear
                    drlamdpsi = drmatdpsi(indvec,:);
                    buckl.drdpsi(8*BPL*(i-1)+(1:8*BPL),:) = drlamdpsi(crit_panelvec,:);
                end
                if morph.span
                    drlamdext = drmatdext(indvec,:);
                    buckl.drdext(8*BPL*(i-1)+(1:8*BPL),:) = drlamdext(crit_panelvec,:);
                end
            end
        end
        
    else
        rlam = rmat(rind,:);
        rlam_crit = max(rlam,[],2);
        
        [~,crit_panel] = sort(rlam_crit,'descend');
        
        crit_panel = sort(crit_panel,'ascend');
        for j=1:ceil(BPL/length(rind))
            crit_panel(length(rind)*(j-1)+(1:length(rind)),1) = crit_panel(1:length(rind));
        end
        
        crit_panel = crit_panel(1:BPL);
        
        buckl.r(8*BPL*(i-1)+(1:8*BPL),1) = reshape((rlam(crit_panel,:))',[],1);
        buckl.pan{i}(:,1) = rind(crit_panel);
        buckl.sec{i}(:,1) = bucklsec(rind(crit_panel));
        buckl.cross{i}(:,1) = bucklcross(rind(crit_panel));
        
        if ders == 1
            indvec = reshape((8*(rind(:,ones(8,1))-1)+repmat(1:8,length(rind),1))',[],1);
            crit_panelvec = reshape((8*(crit_panel(:,ones(8,1))-1)+repmat(1:8,BPL,1))',[],1);
            if tailflag == 1
                drlamdA = drmatdA(indvec,:);
                drlamdD = drmatdD(indvec,:);
                drlamdt = drmatdt(indvec,:);
                
                if gustflag == 1
                    drlamddv = drmatddv(indvec,:);
                end
                
                buckl.drdA(8*BPL*(i-1)+(1:8*BPL),:) = drlamdA(crit_panelvec,:);
                buckl.drdD(8*BPL*(i-1)+(1:8*BPL),:) = drlamdD(crit_panelvec,:);
                buckl.drdt(8*BPL*(i-1)+(1:8*BPL),:) = drlamdt(crit_panelvec,:);
                
                if gustflag == 1
                    buckl.drddv(8*BPL*(i-1)+(1:8*BPL),:) = drlamddv(crit_panelvec,:);
                end
            end
            if morphflag == 1
                if morph.camber
                    drlamdparam = drmatdparam(indvec,:);
                    buckl.drdparam(8*BPL*(i-1)+(1:8*BPL),:) = drlamdparam(crit_panelvec,:);
                end
                if morph.twist
                    drlamdphi = drmatdphi(indvec,:);
                    buckl.drdphi(8*BPL*(i-1)+(1:8*BPL),:) = drlamdphi(crit_panelvec,:);
                end
                if morph.fold
                    drlamdtheta = drmatdtheta(indvec,:);
                    buckl.drdtheta(8*BPL*(i-1)+(1:8*BPL),:) = drlamdtheta(crit_panelvec,:);
                end
                if morph.shear
                    drlamdpsi = drmatdpsi(indvec,:);
                    buckl.drdpsi(8*BPL*(i-1)+(1:8*BPL),:) = drlamdpsi(crit_panelvec,:);
                end
                if morph.span
                    drlamdext = drmatdext(indvec,:);
                    buckl.drdext(8*BPL*(i-1)+(1:8*BPL),:) = drlamdext(crit_panelvec,:);
                end
            end
        end
    end
end

