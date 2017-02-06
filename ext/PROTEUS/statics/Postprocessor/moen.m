function [statics] = moen(constant,statics,ders,tailflag)

nstep = length(statics.morph.sc);

if isfield(statics.morph,'twist')
    if ders == 1
        statics.morph.twist.dEposdphiini = sparse(constant.str.Ns,constant.str.Ns);
        statics.morph.twist.dEposdphifin = sparse(constant.str.Ns,constant.str.Ns);
        statics.morph.twist.dEnegdphiini = sparse(constant.str.Ns,constant.str.Ns);
        statics.morph.twist.dEnegdphifin = sparse(constant.str.Ns,constant.str.Ns);
        
        if isfield(statics.morph,'shear')
            statics.morph.twist.dEposdpsiini = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.twist.dEposdpsifin = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.twist.dEnegdpsiini = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.twist.dEnegdpsifin = sparse(constant.str.Ns,constant.str.Ns);
        end
        
        
        if isfield(statics.morph,'fold')
            statics.morph.twist.dEposdthetaini = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.twist.dEposdthetafin = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.twist.dEnegdthetaini = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.twist.dEnegdthetafin = sparse(constant.str.Ns,constant.str.Ns);
        end
        
        if isfield(statics.morph,'camber')
            statics.morph.twist.dEposdparamini = sparse(constant.str.Ns,length(constant.morph.camber.loc));
            statics.morph.twist.dEposdparamfin = sparse(constant.str.Ns,length(constant.morph.camber.loc));
            statics.morph.twist.dEnegdparamini = sparse(constant.str.Ns,length(constant.morph.camber.loc));
            statics.morph.twist.dEnegdparamfin = sparse(constant.str.Ns,length(constant.morph.camber.loc));
        end
        
        if isfield(statics.morph,'span')
            statics.morph.twist.dEposdextini = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.twist.dEposdextfin = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.twist.dEnegdextini = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.twist.dEnegdextfin = sparse(constant.str.Ns,constant.str.Ns);
        end
        
        if tailflag == 1
            statics.morph.twist.dEposdC = sparse(constant.str.Ns,numel(statics.str.C));
            statics.morph.twist.dEposdt = sparse(constant.str.Ns,length(constant.lam.ID));
            statics.morph.twist.dEnegdC = sparse(constant.str.Ns,numel(statics.str.C));
            statics.morph.twist.dEnegdt = sparse(constant.str.Ns,length(constant.lam.ID));
        end
    end
    
    
    for i = 1:constant.str.Ns
        Edat = (statics.morph.twist.angle(i,2:end)-statics.morph.twist.angle(i,1:end-1)).*(statics.morph.twist.Mphi(i,2:end)+statics.morph.twist.Mphi(i,1:end-1))/2;
        
%         keyboard
        statics.morph.twist.E(i,1) = sum((statics.morph.twist.angle(i,2:end)-statics.morph.twist.angle(i,1:end-1)).*(statics.morph.twist.Mphi(i,2:end)+statics.morph.twist.Mphi(i,1:end-1))/2);
        statics.morph.twist.Epos(i,1) = sum(max(Edat,0));
        statics.morph.twist.Eneg(i,1) = sum(min(Edat,0));
        indpos = find(max(Edat,0));
        indneg = find(min(Edat,0));
%         if isempty(indpos) && isempty(indneg)
%            indpos = 1:nstep-1;
%            indneg = 1:nstep-1;
%         end
        if ders == 1
            statics.morph.twist.dEdphiini(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMphidphiini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMphidphiini(nstep*(i-1)+(1:nstep-1),:))/2);
            statics.morph.twist.dEdphifin(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMphidphifin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMphidphifin(nstep*(i-1)+(1:nstep-1),:))/2);
            
            statics.morph.twist.dEposdphiini(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMphidphiini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMphidphiini(nstep*(i-1)+(indpos),:))/2);
            statics.morph.twist.dEposdphifin(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMphidphifin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMphidphifin(nstep*(i-1)+(indpos),:))/2);
            
            statics.morph.twist.dEnegdphiini(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMphidphiini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMphidphiini(nstep*(i-1)+(indneg),:))/2);
            statics.morph.twist.dEnegdphifin(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMphidphifin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMphidphifin(nstep*(i-1)+(indneg),:))/2);
            
            
            for j=1:nstep-1
                statics.morph.twist.dEdphiini(i,i) = statics.morph.twist.dEdphiini(i,i)+(statics.morph.twist.dangledini{j+1}(i,i)-statics.morph.twist.dangledini{j}(i,i))'.*(statics.morph.twist.Mphi(i,j+1)+statics.morph.twist.Mphi(i,j))/2;
                statics.morph.twist.dEdphifin(i,i) = statics.morph.twist.dEdphifin(i,i)+(statics.morph.twist.dangledfin{j+1}(i,i)-statics.morph.twist.dangledfin{j}(i,i))'.*(statics.morph.twist.Mphi(i,j+1)+statics.morph.twist.Mphi(i,j))/2;
                if ~isempty(find(j==indpos,1))
                    statics.morph.twist.dEposdphiini(i,i) = statics.morph.twist.dEposdphiini(i,i)+(statics.morph.twist.dangledini{j+1}(i,i)-statics.morph.twist.dangledini{j}(i,i))'.*(statics.morph.twist.Mphi(i,j+1)+statics.morph.twist.Mphi(i,j))/2;
                    statics.morph.twist.dEposdphifin(i,i) = statics.morph.twist.dEposdphifin(i,i)+(statics.morph.twist.dangledfin{j+1}(i,i)-statics.morph.twist.dangledfin{j}(i,i))'.*(statics.morph.twist.Mphi(i,j+1)+statics.morph.twist.Mphi(i,j))/2;
                end
                if ~isempty(find(j==indneg,1))
                    statics.morph.twist.dEnegdphiini(i,i) = statics.morph.twist.dEnegdphiini(i,i)+(statics.morph.twist.dangledini{j+1}(i,i)-statics.morph.twist.dangledini{j}(i,i))'.*(statics.morph.twist.Mphi(i,j+1)+statics.morph.twist.Mphi(i,j))/2;
                    statics.morph.twist.dEnegdphifin(i,i) = statics.morph.twist.dEnegdphifin(i,i)+(statics.morph.twist.dangledfin{j+1}(i,i)-statics.morph.twist.dangledfin{j}(i,i))'.*(statics.morph.twist.Mphi(i,j+1)+statics.morph.twist.Mphi(i,j))/2;
                end
            end
            if isfield(statics.morph,'shear')
                statics.morph.twist.dEdpsiini(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMphidpsiini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMphidpsiini(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.twist.dEdpsifin(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMphidpsifin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMphidpsifin(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.twist.dEposdpsiini(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMphidpsiini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMphidpsiini(nstep*(i-1)+(indpos),:))/2);
                statics.morph.twist.dEposdpsifin(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMphidpsifin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMphidpsifin(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.twist.dEnegdpsiini(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMphidpsiini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMphidpsiini(nstep*(i-1)+(indneg),:))/2);
                statics.morph.twist.dEnegdpsifin(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMphidpsifin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMphidpsifin(nstep*(i-1)+(indneg),:))/2);
            end
            if isfield(statics.morph,'fold')
                statics.morph.twist.dEdthetaini(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMphidthetaini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMphidthetaini(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.twist.dEdthetafin(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMphidthetafin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMphidthetafin(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.twist.dEposdthetaini(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMphidthetaini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMphidthetaini(nstep*(i-1)+(indpos),:))/2);
                statics.morph.twist.dEposdthetafin(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMphidthetafin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMphidthetafin(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.twist.dEnegdthetaini(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMphidthetaini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMphidthetaini(nstep*(i-1)+(indneg),:))/2);
                statics.morph.twist.dEnegdthetafin(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMphidthetafin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMphidthetafin(nstep*(i-1)+(indneg),:))/2);
            
            end
            if isfield(statics.morph,'camber')
                statics.morph.twist.dEdparamini(i,:) = sum((statics.morph.twist.angle(i*ones(length(constant.morph.camber.loc),1),2:end)'-statics.morph.twist.angle(i*ones(length(constant.morph.camber.loc),1),1:end-1)').*(statics.sens.dMphidparamini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMphidparamini(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.twist.dEdparamfin(i,:) = sum((statics.morph.twist.angle(i*ones(length(constant.morph.camber.loc),1),2:end)'-statics.morph.twist.angle(i*ones(length(constant.morph.camber.loc),1),1:end-1)').*(statics.sens.dMphidparamfin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMphidparamfin(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.twist.dEposdparamini(i,:) = sum((statics.morph.twist.angle(i*ones(length(constant.morph.camber.loc),1),indpos+1)'-statics.morph.twist.angle(i*ones(length(constant.morph.camber.loc),1),indpos)').*(statics.sens.dMphidparamini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMphidparamini(nstep*(i-1)+(indpos),:))/2);
                statics.morph.twist.dEposdparamfin(i,:) = sum((statics.morph.twist.angle(i*ones(length(constant.morph.camber.loc),1),indpos+1)'-statics.morph.twist.angle(i*ones(length(constant.morph.camber.loc),1),indpos)').*(statics.sens.dMphidparamfin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMphidparamfin(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.twist.dEnegdparamini(i,:) = sum((statics.morph.twist.angle(i*ones(length(constant.morph.camber.loc),1),indneg+1)'-statics.morph.twist.angle(i*ones(length(constant.morph.camber.loc),1),indneg)').*(statics.sens.dMphidparamini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMphidparamini(nstep*(i-1)+(indneg),:))/2);
                statics.morph.twist.dEnegdparamfin(i,:) = sum((statics.morph.twist.angle(i*ones(length(constant.morph.camber.loc),1),indneg+1)'-statics.morph.twist.angle(i*ones(length(constant.morph.camber.loc),1),indneg)').*(statics.sens.dMphidparamfin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMphidparamfin(nstep*(i-1)+(indneg),:))/2);
            
            end
            if isfield(statics.morph,'span')
                statics.morph.twist.dEdextini(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMphidextini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMphidextini(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.twist.dEdextfin(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMphidextfin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMphidextfin(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.twist.dEposdextini(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMphidextini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMphidextini(nstep*(i-1)+(indpos),:))/2);
                statics.morph.twist.dEposdextfin(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMphidextfin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMphidextfin(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.twist.dEnegdextini(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMphidextini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMphidextini(nstep*(i-1)+(indneg),:))/2);
                statics.morph.twist.dEnegdextfin(i,:) = sum((statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.twist.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMphidextfin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMphidextfin(nstep*(i-1)+(indneg),:))/2);
            
            end
            if tailflag == 1
                statics.morph.twist.dEdC(i,:) = sum((statics.morph.twist.angle(i*ones(numel(statics.str.C),1),2:end)'-statics.morph.twist.angle(i*ones(numel(statics.str.C),1),1:end-1)').*(statics.sens.dMphidC(nstep*(i-1)+(2:nstep),:)+statics.sens.dMphidC(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.twist.dEdt(i,:) = sum((statics.morph.twist.angle(i*ones(length(constant.lam.ID),1),2:end)'-statics.morph.twist.angle(i*ones(length(constant.lam.ID),1),1:end-1)').*(statics.sens.dMphidt(nstep*(i-1)+(2:nstep),:)+statics.sens.dMphidt(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.twist.dEposdC(i,:) = sum((statics.morph.twist.angle(i*ones(numel(statics.str.C),1),indpos+1)'-statics.morph.twist.angle(i*ones(numel(statics.str.C),1),indpos)').*(statics.sens.dMphidC(nstep*(i-1)+(indpos+1),:)+statics.sens.dMphidC(nstep*(i-1)+(indpos),:))/2);
                statics.morph.twist.dEposdt(i,:) = sum((statics.morph.twist.angle(i*ones(length(constant.lam.ID),1),indpos+1)'-statics.morph.twist.angle(i*ones(length(constant.lam.ID),1),indpos)').*(statics.sens.dMphidt(nstep*(i-1)+(indpos+1),:)+statics.sens.dMphidt(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.twist.dEnegdC(i,:) = sum((statics.morph.twist.angle(i*ones(numel(statics.str.C),1),indneg+1)'-statics.morph.twist.angle(i*ones(numel(statics.str.C),1),indneg)').*(statics.sens.dMphidC(nstep*(i-1)+(indneg+1),:)+statics.sens.dMphidC(nstep*(i-1)+(indneg),:))/2);
                statics.morph.twist.dEnegdt(i,:) = sum((statics.morph.twist.angle(i*ones(length(constant.lam.ID),1),indneg+1)'-statics.morph.twist.angle(i*ones(length(constant.lam.ID),1),indneg)').*(statics.sens.dMphidt(nstep*(i-1)+(indneg+1),:)+statics.sens.dMphidt(nstep*(i-1)+(indneg),:))/2);
            
            end
        end
    end
end
if isfield(statics.morph,'shear')
    
    if ders == 1
        statics.morph.shear.dEposdpsiini = sparse(constant.str.Ns,constant.str.Ns);
        statics.morph.shear.dEposdpsifin = sparse(constant.str.Ns,constant.str.Ns);
        statics.morph.shear.dEnegdpsiini = sparse(constant.str.Ns,constant.str.Ns);
        statics.morph.shear.dEnegdpsifin = sparse(constant.str.Ns,constant.str.Ns);
        
        if isfield(statics.morph,'twist')
            statics.morph.shear.dEposdphiini = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.shear.dEposdphifin = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.shear.dEnegdphiini = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.shear.dEnegdphifin = sparse(constant.str.Ns,constant.str.Ns);
        end
        
        
        if isfield(statics.morph,'fold')
            statics.morph.shear.dEposdthetaini = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.shear.dEposdthetafin = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.shear.dEnegdthetaini = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.shear.dEnegdthetafin = sparse(constant.str.Ns,constant.str.Ns);
        end
        
        if isfield(statics.morph,'camber')
            statics.morph.shear.dEposdparamini = sparse(constant.str.Ns,length(constant.morph.camber.loc));
            statics.morph.shear.dEposdparamfin = sparse(constant.str.Ns,length(constant.morph.camber.loc));
            statics.morph.shear.dEnegdparamini = sparse(constant.str.Ns,length(constant.morph.camber.loc));
            statics.morph.shear.dEnegdparamfin = sparse(constant.str.Ns,length(constant.morph.camber.loc));
        end
        
        if isfield(statics.morph,'span')
            statics.morph.shear.dEposdextini = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.shear.dEposdextfin = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.shear.dEnegdextini = sparse(constant.str.Ns,constant.str.Ns);
            statics.morph.shear.dEnegdextfin = sparse(constant.str.Ns,constant.str.Ns);
        end
        
        if tailflag == 1
            statics.morph.shear.dEposdC = sparse(constant.str.Ns,numel(statics.str.C));
            statics.morph.shear.dEposdt = sparse(constant.str.Ns,length(constant.lam.ID));
            statics.morph.shear.dEnegdC = sparse(constant.str.Ns,numel(statics.str.C));
            statics.morph.shear.dEnegdt = sparse(constant.str.Ns,length(constant.lam.ID));
        end
    end
    
    for i = 1:constant.str.Ns
        
        Edat = (statics.morph.shear.angle(i,2:end)-statics.morph.shear.angle(i,1:end-1)).*(statics.morph.shear.Mpsi(i,2:end)+statics.morph.shear.Mpsi(i,1:end-1))/2;
        
        statics.morph.shear.E(i,1) = sum((statics.morph.shear.angle(i,2:end)-statics.morph.shear.angle(i,1:end-1)).*(statics.morph.shear.Mpsi(i,2:end)+statics.morph.shear.Mpsi(i,1:end-1))/2);
        statics.morph.shear.Epos(i,1) = sum(max(Edat,0));
        statics.morph.shear.Eneg(i,1) = sum(min(Edat,0));
        indpos = find(max(Edat,0));
        indneg = find(min(Edat,0));
        
        if ders == 1
            statics.morph.shear.dEdpsiini(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMpsidpsiini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMpsidpsiini(nstep*(i-1)+(1:nstep-1),:))/2);
            statics.morph.shear.dEdpsifin(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMpsidpsifin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMpsidpsifin(nstep*(i-1)+(1:nstep-1),:))/2);
            
            statics.morph.shear.dEposdpsiini(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMpsidpsiini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMpsidpsiini(nstep*(i-1)+(indpos),:))/2);
            statics.morph.shear.dEposdpsifin(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMpsidpsifin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMpsidpsifin(nstep*(i-1)+(indpos),:))/2);
            
            statics.morph.shear.dEnegdpsiini(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMpsidpsiini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMpsidpsiini(nstep*(i-1)+(indneg),:))/2);
            statics.morph.shear.dEnegdpsifin(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMpsidpsifin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMpsidpsifin(nstep*(i-1)+(indneg),:))/2);
            
            
            for j=1:nstep-1
                statics.morph.shear.dEdpsiini(i,i) = statics.morph.shear.dEdpsiini(i,i)+(statics.morph.shear.dangledini{j+1}(i,i)-statics.morph.shear.dangledini{j}(i,i))'.*(statics.morph.shear.Mpsi(i,j+1)+statics.morph.shear.Mpsi(i,j))/2;
                statics.morph.shear.dEdpsifin(i,i) = statics.morph.shear.dEdpsifin(i,i)+(statics.morph.shear.dangledfin{j+1}(i,i)-statics.morph.shear.dangledfin{j}(i,i))'.*(statics.morph.shear.Mpsi(i,j+1)+statics.morph.shear.Mpsi(i,j))/2;
                if ~isempty(find(j==indpos,1))
                    statics.morph.shear.dEposdpsiini(i,i) = statics.morph.shear.dEposdpsiini(i,i)+(statics.morph.shear.dangledini{j+1}(i,i)-statics.morph.shear.dangledini{j}(i,i))'.*(statics.morph.shear.Mpsi(i,j+1)+statics.morph.shear.Mpsi(i,j))/2;
                    statics.morph.shear.dEposdpsifin(i,i) = statics.morph.shear.dEposdpsifin(i,i)+(statics.morph.shear.dangledfin{j+1}(i,i)-statics.morph.shear.dangledfin{j}(i,i))'.*(statics.morph.shear.Mpsi(i,j+1)+statics.morph.shear.Mpsi(i,j))/2;
                end
                if ~isempty(find(j==indneg,1))
                    statics.morph.shear.dEnegdpsiini(i,i) = statics.morph.shear.dEnegdpsiini(i,i)+(statics.morph.shear.dangledini{j+1}(i,i)-statics.morph.shear.dangledini{j}(i,i))'.*(statics.morph.shear.Mpsi(i,j+1)+statics.morph.shear.Mpsi(i,j))/2;
                    statics.morph.shear.dEnegdpsifin(i,i) = statics.morph.shear.dEnegdpsifin(i,i)+(statics.morph.shear.dangledfin{j+1}(i,i)-statics.morph.shear.dangledfin{j}(i,i))'.*(statics.morph.shear.Mpsi(i,j+1)+statics.morph.shear.Mpsi(i,j))/2;
                end
            end
            if isfield(statics.morph,'twist')
                statics.morph.shear.dEdphiini(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMpsidphiini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMpsidphiini(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.shear.dEdphifin(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMpsidphifin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMpsidphifin(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.shear.dEposdphiini(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMpsidphiini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMpsidphiini(nstep*(i-1)+(indpos),:))/2);
                statics.morph.shear.dEposdphifin(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMpsidphifin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMpsidphifin(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.shear.dEnegdphiini(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMpsidphiini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMpsidphiini(nstep*(i-1)+(indneg),:))/2);
                statics.morph.shear.dEnegdphifin(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMpsidphifin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMpsidphifin(nstep*(i-1)+(indneg),:))/2);
            
            end
            if isfield(statics.morph,'fold')
                statics.morph.shear.dEdthetaini(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMpsidthetaini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMpsidthetaini(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.shear.dEdthetafin(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMpsidthetafin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMpsidthetafin(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.shear.dEposdthetaini(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMpsidthetaini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMpsidthetaini(nstep*(i-1)+(indpos),:))/2);
                statics.morph.shear.dEposdthetafin(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMpsidthetafin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMpsidthetafin(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.shear.dEnegdthetaini(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMpsidthetaini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMpsidthetaini(nstep*(i-1)+(indneg),:))/2);
                statics.morph.shear.dEnegdthetafin(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMpsidthetafin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMpsidthetafin(nstep*(i-1)+(indneg),:))/2);
            
            end
            if isfield(statics.morph,'camber')
                statics.morph.shear.dEdparamini(i,:) = sum((statics.morph.shear.angle(i*ones(length(constant.morph.camber.loc),1),2:end)'-statics.morph.shear.angle(i*ones(length(constant.morph.camber.loc),1),1:end-1)').*(statics.sens.dMpsidparamini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMpsidparamini(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.shear.dEdparamfin(i,:) = sum((statics.morph.shear.angle(i*ones(length(constant.morph.camber.loc),1),2:end)'-statics.morph.shear.angle(i*ones(length(constant.morph.camber.loc),1),1:end-1)').*(statics.sens.dMpsidparamfin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMpsidparamfin(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.shear.dEposdparamini(i,:) = sum((statics.morph.shear.angle(i*ones(length(constant.morph.camber.loc),1),indpos+1)'-statics.morph.shear.angle(i*ones(length(constant.morph.camber.loc),1),indpos)').*(statics.sens.dMpsidparamini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMpsidparamini(nstep*(i-1)+(indpos),:))/2);
                statics.morph.shear.dEposdparamfin(i,:) = sum((statics.morph.shear.angle(i*ones(length(constant.morph.camber.loc),1),indpos+1)'-statics.morph.shear.angle(i*ones(length(constant.morph.camber.loc),1),indpos)').*(statics.sens.dMpsidparamfin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMpsidparamfin(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.shear.dEnegdparamini(i,:) = sum((statics.morph.shear.angle(i*ones(length(constant.morph.camber.loc),1),indneg+1)'-statics.morph.shear.angle(i*ones(length(constant.morph.camber.loc),1),indneg)').*(statics.sens.dMpsidparamini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMpsidparamini(nstep*(i-1)+(indneg),:))/2);
                statics.morph.shear.dEnegdparamfin(i,:) = sum((statics.morph.shear.angle(i*ones(length(constant.morph.camber.loc),1),indneg+1)'-statics.morph.shear.angle(i*ones(length(constant.morph.camber.loc),1),indneg)').*(statics.sens.dMpsidparamfin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMpsidparamfin(nstep*(i-1)+(indneg),:))/2);
            
            end
            if isfield(statics.morph,'span')
                statics.morph.shear.dEdextini(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMpsidextini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMpsidextini(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.shear.dEdextfin(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMpsidextfin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMpsidextfin(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.shear.dEposdextini(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMpsidextini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMpsidextini(nstep*(i-1)+(indpos),:))/2);
                statics.morph.shear.dEposdextfin(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMpsidextfin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMpsidextfin(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.shear.dEnegdextini(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMpsidextini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMpsidextini(nstep*(i-1)+(indneg),:))/2);
                statics.morph.shear.dEnegdextfin(i,:) = sum((statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.shear.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMpsidextfin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMpsidextfin(nstep*(i-1)+(indneg),:))/2);
            
            end
            if tailflag == 1
                statics.morph.shear.dEdC(i,:) = sum((statics.morph.shear.angle(i*ones(numel(statics.str.C),1),2:end)'-statics.morph.shear.angle(i*ones(numel(statics.str.C),1),1:end-1)').*(statics.sens.dMpsidC(nstep*(i-1)+(2:nstep),:)+statics.sens.dMpsidC(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.shear.dEdt(i,:) = sum((statics.morph.shear.angle(i*ones(length(constant.lam.ID),1),2:end)'-statics.morph.shear.angle(i*ones(length(constant.lam.ID),1),1:end-1)').*(statics.sens.dMpsidt(nstep*(i-1)+(2:nstep),:)+statics.sens.dMpsidt(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.shear.dEposdC(i,:) = sum((statics.morph.shear.angle(i*ones(numel(statics.str.C),1),indpos+1)'-statics.morph.shear.angle(i*ones(numel(statics.str.C),1),indpos)').*(statics.sens.dMpsidC(nstep*(i-1)+(indpos+1),:)+statics.sens.dMpsidC(nstep*(i-1)+(indpos),:))/2);
                statics.morph.shear.dEposdt(i,:) = sum((statics.morph.shear.angle(i*ones(length(constant.lam.ID),1),indpos+1)'-statics.morph.shear.angle(i*ones(length(constant.lam.ID),1),indpos)').*(statics.sens.dMpsidt(nstep*(i-1)+(indpos+1),:)+statics.sens.dMpsidt(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.shear.dEnegdC(i,:) = sum((statics.morph.shear.angle(i*ones(numel(statics.str.C),1),indneg+1)'-statics.morph.shear.angle(i*ones(numel(statics.str.C),1),indneg)').*(statics.sens.dMpsidC(nstep*(i-1)+(indneg+1),:)+statics.sens.dMpsidC(nstep*(i-1)+(indneg),:))/2);
                statics.morph.shear.dEnegdt(i,:) = sum((statics.morph.shear.angle(i*ones(length(constant.lam.ID),1),indneg+1)'-statics.morph.shear.angle(i*ones(length(constant.lam.ID),1),indneg)').*(statics.sens.dMpsidt(nstep*(i-1)+(indneg+1),:)+statics.sens.dMpsidt(nstep*(i-1)+(indneg),:))/2);
            
            end
        end
    end
end
if isfield(statics.morph,'fold')
    for i = 1:constant.str.Ns
        Edat = (statics.morph.fold.angle(i,2:end)-statics.morph.fold.angle(i,1:end-1)).*(statics.morph.fold.Mtheta(i,2:end)+statics.morph.fold.Mtheta(i,1:end-1))/2;
        
        statics.morph.fold.E(i,1) = sum((statics.morph.fold.angle(i,2:end)-statics.morph.fold.angle(i,1:end-1)).*(statics.morph.fold.Mtheta(i,2:end)+statics.morph.fold.Mtheta(i,1:end-1))/2);
        
        statics.morph.fold.Epos(i,1) = sum(max(Edat,0));
        statics.morph.fold.Eneg(i,1) = sum(min(Edat,0));
        indpos = find(max(Edat,0));
        indneg = find(min(Edat,0));
        
        
        if ders == 1
            statics.morph.fold.dEdthetaini(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMthetadthetaini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMthetadthetaini(nstep*(i-1)+(1:nstep-1),:))/2);
            statics.morph.fold.dEdthetafin(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMthetadthetafin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMthetadthetafin(nstep*(i-1)+(1:nstep-1),:))/2);
            
            statics.morph.fold.dEposdthetaini(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMthetadthetaini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMthetadthetaini(nstep*(i-1)+(indpos),:))/2);
            statics.morph.fold.dEposdthetafin(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMthetadthetafin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMthetadthetafin(nstep*(i-1)+(indpos),:))/2);
            
            statics.morph.fold.dEnegdthetaini(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMthetadthetaini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMthetadthetaini(nstep*(i-1)+(indneg),:))/2);
            statics.morph.fold.dEnegdthetafin(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMthetadthetafin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMthetadthetafin(nstep*(i-1)+(indneg),:))/2);
            
            for j=1:nstep-1
                statics.morph.fold.dEdthetaini(i,i) = statics.morph.fold.dEdthetaini(i,i)+(statics.morph.fold.dangledini{j+1}(i,i)-statics.morph.fold.dangledini{j}(i,i))'.*(statics.morph.fold.Mtheta(i,j+1)+statics.morph.fold.Mtheta(i,j))/2;
                statics.morph.fold.dEdthetafin(i,i) = statics.morph.fold.dEdthetafin(i,i)+(statics.morph.fold.dangledfin{j+1}(i,i)-statics.morph.fold.dangledfin{j}(i,i))'.*(statics.morph.fold.Mtheta(i,j+1)+statics.morph.fold.Mtheta(i,j))/2;
                
                if ~isempty(find(j==indpos,1))
                    statics.morph.fold.dEposdthetaini(i,i) = statics.morph.fold.dEposdthetaini(i,i)+(statics.morph.fold.dangledini{j+1}(i,i)-statics.morph.fold.dangledini{j}(i,i))'.*(statics.morph.fold.Mtheta(i,j+1)+statics.morph.fold.Mtheta(i,j))/2;
                    statics.morph.fold.dEposdthetafin(i,i) = statics.morph.fold.dEposdthetafin(i,i)+(statics.morph.fold.dangledfin{j+1}(i,i)-statics.morph.fold.dangledfin{j}(i,i))'.*(statics.morph.fold.Mtheta(i,j+1)+statics.morph.fold.Mtheta(i,j))/2;
                end
                if ~isempty(find(j==indneg,1))
                    statics.morph.fold.dEnegdthetaini(i,i) = statics.morph.fold.dEnegdthetaini(i,i)+(statics.morph.fold.dangledini{j+1}(i,i)-statics.morph.fold.dangledini{j}(i,i))'.*(statics.morph.fold.Mtheta(i,j+1)+statics.morph.fold.Mtheta(i,j))/2;
                    statics.morph.fold.dEnegdthetafin(i,i) = statics.morph.fold.dEnegdthetafin(i,i)+(statics.morph.fold.dangledfin{j+1}(i,i)-statics.morph.fold.dangledfin{j}(i,i))'.*(statics.morph.fold.Mtheta(i,j+1)+statics.morph.fold.Mtheta(i,j))/2;
                end
            end
            if isfield(statics.morph,'twist')
                statics.morph.fold.dEdphiini(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMthetadphiini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMthetadphiini(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.fold.dEdphifin(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMthetadphifin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMthetadphifin(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.fold.dEposdphiini(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMthetadphiini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMthetadphiini(nstep*(i-1)+(indpos),:))/2);
                statics.morph.fold.dEposdphifin(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMthetadphifin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMthetadphifin(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.fold.dEnegdphiini(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMthetadphiini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMthetadphiini(nstep*(i-1)+(indneg),:))/2);
                statics.morph.fold.dEnegdphifin(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMthetadphifin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMthetadphifin(nstep*(i-1)+(indneg),:))/2);
            end
            if isfield(statics.morph,'shear')
                statics.morph.fold.dEdpsiini(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMthetadpsiini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMthetadpsiini(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.fold.dEdpsifin(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMthetadpsifin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMthetadpsifin(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.fold.dEposdpsiini(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMthetadpsiini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMthetadpsiini(nstep*(i-1)+(indpos),:))/2);
                statics.morph.fold.dEposdpsifin(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMthetadpsifin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMthetadpsifin(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.fold.dEnegdpsiini(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMthetadpsiini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMthetadpsiini(nstep*(i-1)+(indneg),:))/2);
                statics.morph.fold.dEnegdpsifin(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMthetadpsifin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMthetadpsifin(nstep*(i-1)+(indneg),:))/2);
            end
            if isfield(statics.morph,'camber')
                statics.morph.fold.dEdparamini(i,:) = sum((statics.morph.fold.angle(i*ones(length(constant.morph.camber.loc),1),2:end)'-statics.morph.fold.angle(i*ones(length(constant.morph.camber.loc),1),1:end-1)').*(statics.sens.dMthetadparamini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMthetadparamini(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.fold.dEdparamfin(i,:) = sum((statics.morph.fold.angle(i*ones(length(constant.morph.camber.loc),1),2:end)'-statics.morph.fold.angle(i*ones(length(constant.morph.camber.loc),1),1:end-1)').*(statics.sens.dMthetadparamfin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMthetadparamfin(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.fold.dEposdparamini(i,:) = sum((statics.morph.fold.angle(i*ones(length(constant.morph.camber.loc),1),indpos+1)'-statics.morph.fold.angle(i*ones(length(constant.morph.camber.loc),1),indpos)').*(statics.sens.dMthetadparamini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMthetadparamini(nstep*(i-1)+(indpos),:))/2);
                statics.morph.fold.dEposdparamfin(i,:) = sum((statics.morph.fold.angle(i*ones(length(constant.morph.camber.loc),1),indpos+1)'-statics.morph.fold.angle(i*ones(length(constant.morph.camber.loc),1),indpos)').*(statics.sens.dMthetadparamfin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMthetadparamfin(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.fold.dEnegdparamini(i,:) = sum((statics.morph.fold.angle(i*ones(length(constant.morph.camber.loc),1),indneg+1)'-statics.morph.fold.angle(i*ones(length(constant.morph.camber.loc),1),indneg)').*(statics.sens.dMthetadparamini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMthetadparamini(nstep*(i-1)+(indneg),:))/2);
                statics.morph.fold.dEnegdparamfin(i,:) = sum((statics.morph.fold.angle(i*ones(length(constant.morph.camber.loc),1),indneg+1)'-statics.morph.fold.angle(i*ones(length(constant.morph.camber.loc),1),indneg)').*(statics.sens.dMthetadparamfin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMthetadparamfin(nstep*(i-1)+(indneg),:))/2);
            end
            if isfield(statics.morph,'span')
                statics.morph.fold.dEdextini(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMthetadextini(nstep*(i-1)+(2:nstep),:)+statics.sens.dMthetadextini(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.fold.dEdextfin(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),2:end)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),1:end-1)').*(statics.sens.dMthetadextfin(nstep*(i-1)+(2:nstep),:)+statics.sens.dMthetadextfin(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.fold.dEposdextini(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMthetadextini(nstep*(i-1)+(indpos+1),:)+statics.sens.dMthetadextini(nstep*(i-1)+(indpos),:))/2);
                statics.morph.fold.dEposdextfin(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indpos)').*(statics.sens.dMthetadextfin(nstep*(i-1)+(indpos+1),:)+statics.sens.dMthetadextfin(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.fold.dEnegdextini(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMthetadextini(nstep*(i-1)+(indneg+1),:)+statics.sens.dMthetadextini(nstep*(i-1)+(indneg),:))/2);
                statics.morph.fold.dEnegdextfin(i,:) = sum((statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg+1)'-statics.morph.fold.angle(i*ones(constant.str.Ns,1),indneg)').*(statics.sens.dMthetadextfin(nstep*(i-1)+(indneg+1),:)+statics.sens.dMthetadextfin(nstep*(i-1)+(indneg),:))/2);
            end
            if tailflag == 1
                statics.morph.fold.dEdC(i,:) = sum((statics.morph.fold.angle(i*ones(numel(statics.str.C),1),2:end)'-statics.morph.fold.angle(i*ones(numel(statics.str.C),1),1:end-1)').*(statics.sens.dMthetadC(nstep*(i-1)+(2:nstep),:)+statics.sens.dMthetadC(nstep*(i-1)+(1:nstep-1),:))/2);
                statics.morph.fold.dEdt(i,:) = sum((statics.morph.fold.angle(i*ones(length(constant.lam.ID),1),2:end)'-statics.morph.fold.angle(i*ones(length(constant.lam.ID),1),1:end-1)').*(statics.sens.dMthetadt(nstep*(i-1)+(2:nstep),:)+statics.sens.dMthetadt(nstep*(i-1)+(1:nstep-1),:))/2);
                
                statics.morph.fold.dEposdC(i,:) = sum((statics.morph.fold.angle(i*ones(numel(statics.str.C),1),indpos+1)'-statics.morph.fold.angle(i*ones(numel(statics.str.C),1),indpos)').*(statics.sens.dMthetadC(nstep*(i-1)+(indpos+1),:)+statics.sens.dMthetadC(nstep*(i-1)+(indpos),:))/2);
                statics.morph.fold.dEposdt(i,:) = sum((statics.morph.fold.angle(i*ones(length(constant.lam.ID),1),indpos+1)'-statics.morph.fold.angle(i*ones(length(constant.lam.ID),1),indpos)').*(statics.sens.dMthetadt(nstep*(i-1)+(indpos+1),:)+statics.sens.dMthetadt(nstep*(i-1)+(indpos),:))/2);
                
                statics.morph.fold.dEnegdC(i,:) = sum((statics.morph.fold.angle(i*ones(numel(statics.str.C),1),indneg+1)'-statics.morph.fold.angle(i*ones(numel(statics.str.C),1),indneg)').*(statics.sens.dMthetadC(nstep*(i-1)+(indneg+1),:)+statics.sens.dMthetadC(nstep*(i-1)+(indneg),:))/2);
                statics.morph.fold.dEnegdt(i,:) = sum((statics.morph.fold.angle(i*ones(length(constant.lam.ID),1),indneg+1)'-statics.morph.fold.angle(i*ones(length(constant.lam.ID),1),indneg)').*(statics.sens.dMthetadt(nstep*(i-1)+(indneg+1),:)+statics.sens.dMthetadt(nstep*(i-1)+(indneg),:))/2);
            end
        end
    end
end


if isfield(statics.morph,'camber')
    statics.morph.camber.E = zeros(constant.aero.nbs,1);
    statics.morph.camber.Epos = zeros(constant.aero.nbs,1);
    statics.morph.camber.Eneg = zeros(constant.aero.nbs,1);
    if ders == 1
        statics.morph.camber.dEdparamini = sparse(constant.aero.nbs,length(constant.morph.camber.loc));
        statics.morph.camber.dEdparamfin = sparse(constant.aero.nbs,length(constant.morph.camber.loc));
        
        statics.morph.camber.dEposdparamini = sparse(constant.aero.nbs,length(constant.morph.camber.loc));
        statics.morph.camber.dEposdparamfin = sparse(constant.aero.nbs,length(constant.morph.camber.loc));
        
        statics.morph.camber.dEnegdparamini = sparse(constant.aero.nbs,length(constant.morph.camber.loc));
        statics.morph.camber.dEnegdparamfin = sparse(constant.aero.nbs,length(constant.morph.camber.loc));
        
        
        if isfield(statics.morph,'shear')
            statics.morph.camber.dEdpsiini = sparse(constant.aero.nbs,constant.str.Ns);
            statics.morph.camber.dEdpsifin = sparse(constant.aero.nbs,constant.str.Ns);
            
            statics.morph.camber.dEposdpsiini = sparse(constant.aero.nbs,constant.str.Ns);
            statics.morph.camber.dEposdpsifin = sparse(constant.aero.nbs,constant.str.Ns);
            
            statics.morph.camber.dEnegdpsiini = sparse(constant.aero.nbs,constant.str.Ns);
            statics.morph.camber.dEnegdpsifin = sparse(constant.aero.nbs,constant.str.Ns);
        end
        if isfield(statics.morph,'fold')
            statics.morph.camber.dEdthetaini = sparse(constant.aero.nbs,constant.str.Ns);
            statics.morph.camber.dEdthetafin = sparse(constant.aero.nbs,constant.str.Ns);
            
            statics.morph.camber.dEposdthetaini = sparse(constant.aero.nbs,constant.str.Ns);
            statics.morph.camber.dEposdthetafin = sparse(constant.aero.nbs,constant.str.Ns);
            
            statics.morph.camber.dEnegdthetaini = sparse(constant.aero.nbs,constant.str.Ns);
            statics.morph.camber.dEnegdthetafin = sparse(constant.aero.nbs,constant.str.Ns);
        end
        if isfield(statics.morph,'twist')
            statics.morph.camber.dEdphiini = sparse(constant.aero.nbs,constant.str.Ns);
            statics.morph.camber.dEdphifin = sparse(constant.aero.nbs,constant.str.Ns);
            
            statics.morph.camber.dEposdphiini = sparse(constant.aero.nbs,constant.str.Ns);
            statics.morph.camber.dEposdphifin = sparse(constant.aero.nbs,constant.str.Ns);
            
            statics.morph.camber.dEnegdphiini = sparse(constant.aero.nbs,constant.str.Ns);
            statics.morph.camber.dEnegdphifin = sparse(constant.aero.nbs,constant.str.Ns);
        end
        if isfield(statics.morph,'span')
            statics.morph.camber.dEdextini = sparse(constant.aero.nbs,constant.str.Ns);
            statics.morph.camber.dEdextfin = sparse(constant.aero.nbs,constant.str.Ns);
            
            statics.morph.camber.dEposdextini = sparse(constant.aero.nbs,constant.str.Ns);
            statics.morph.camber.dEposdextfin = sparse(constant.aero.nbs,constant.str.Ns);
            
            statics.morph.camber.dEnegdextini = sparse(constant.aero.nbs,constant.str.Ns);
            statics.morph.camber.dEnegdextfin = sparse(constant.aero.nbs,constant.str.Ns);
        end
        if tailflag == 1
            statics.morph.camber.dEdC = sparse(constant.aero.nbs,numel(statics.str.C));
            statics.morph.camber.dEdt = sparse(constant.aero.nbs,length(constant.lam.ID));
            
            statics.morph.camber.dEposdC = sparse(constant.aero.nbs,numel(statics.str.C));
            statics.morph.camber.dEposdt = sparse(constant.aero.nbs,length(constant.lam.ID));
            
            statics.morph.camber.dEnegdC = sparse(constant.aero.nbs,numel(statics.str.C));
            statics.morph.camber.dEnegdt = sparse(constant.aero.nbs,length(constant.lam.ID));
        end
    end
    
    
    nsteps = size(statics.morph.camber.Faloc,2)-1;
    for i = 1:size(statics.morph.camber.Fadist,1)/3
        % Find the section of the aerodynamic force
        Fsec = rem((i-1),constant.aero.nbs)+1;
        
        % Find the corresponding row of the aerodynamic force
        Frow = ceil(i/constant.aero.nbs);
        
        % Integrate to find dr
        rowstart = find(statics.morph.camber.xloccam{1}(:,1) == statics.morph.camber.xlocstat{1}(Frow,1));
        rowend = find(statics.morph.camber.xloccam{1}(:,1) == statics.morph.camber.xlocstat{1}(Frow+1,1));
        
%         drx = trapz(statics.morph.camber.xloccam{1}(rowstart:rowend,Fsec),...
%         (statics.morph.camber.Faloc(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+1:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+1,2:end)-...
%         statics.morph.camber.Faloc(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+1:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+1,1:end-1)))/...
%         (statics.morph.camber.xloccam{1}(rowend,Fsec)-statics.morph.camber.xloccam{1}(rowstart,Fsec));
        
        dFalocmatx = (statics.morph.camber.Faloc(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+1:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+1,2:end)-...
        statics.morph.camber.Faloc(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+1:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+1,1:end-1));

        drx = sum((statics.morph.camber.xloccam{1}(rowstart+1:rowend,ones(nsteps,1))-statics.morph.camber.xloccam{1}(rowstart:rowend-1,ones(nsteps,1))).*...
            (dFalocmatx(1:end-1,:)+dFalocmatx(2:end,:))/2)/...
        (statics.morph.camber.xloccam{1}(rowend,1)-statics.morph.camber.xloccam{1}(rowstart,1));
        
%         drz = trapz(statics.morph.camber.xloccam{1}(rowstart:rowend,Fsec),...
%         (statics.morph.camber.Faloc(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+3:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+3,2:end)-...
%         statics.morph.camber.Faloc(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+3:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+3,1:end-1)))/...
%         (statics.morph.camber.xloccam{1}(rowend,Fsec)-statics.morph.camber.xloccam{1}(rowstart,Fsec));
        
        dFalocmatz = (statics.morph.camber.Faloc(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+3:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+3,2:end)-...
        statics.morph.camber.Faloc(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+3:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+3,1:end-1));
        
        drz = sum((statics.morph.camber.xloccam{1}(rowstart+1:rowend,ones(nsteps,1))-statics.morph.camber.xloccam{1}(rowstart:rowend-1,ones(nsteps,1))).*...
            (dFalocmatz(1:end-1,:)+dFalocmatz(2:end,:))/2)/...
        (statics.morph.camber.xloccam{1}(rowend,1)-statics.morph.camber.xloccam{1}(rowstart,1));
    
    
    
%         % Define the undeformed camber motion
%         dr = [statics.morph.camber.Faloc(3*(i-1)+1,2:end)-statics.morph.camber.Faloc(3*(i-1)+1,1:end-1);
%             zeros(1,size(statics.morph.camber.Faloc,2)-1);
%             statics.morph.camber.Faloc(3*(i-1)+3,2:end)-statics.morph.camber.Faloc(3*(i-1)+3,1:end-1)];
        
        % Define the undeformed camber motion
        dr = [drx;
            zeros(1,nsteps);
            drz];
        
        for j = 1:nstep
            % Rotate dr with the section angle
            Rdir = statics.morph.camber.Farot{j}(3*(Fsec-1)+(1:3),:);
            Fa(:,j) = Rdir'*statics.morph.camber.Fadist(3*(i-1)+(1:3),j);
        end
        
        statics.morph.camber.forces.Fdist{Fsec}(3*(Frow-1)+(1:3),:) = Fa;
        statics.morph.camber.forces.dr{Fsec}(3*(Frow-1)+(1:3),:) = dr;
        statics.morph.camber.forces.Floc{Fsec}(3*(Frow-1)+(1:3),:) = statics.morph.camber.Falocstat(3*(i-1)+(1:3),:);
        
        % Average the load of a morphing step 
        Fa = (Fa(:,2:end)+Fa(:,1:end-1))/2;
        % take care that the work that needs to be delivered is the
        % negative of the aerodynamic work, since sum(W) = 0
%         statics.morph.camber.E(Fsec,1) = -sum(sum(Fa.*dr));

        Edat = -sum(Fa.*dr);

        statics.morph.camber.E(Fsec,1) = statics.morph.camber.E(Fsec,1)-sum(sum(Fa.*dr));
        
        statics.morph.camber.Epos(Fsec,1) = statics.morph.camber.Epos(Fsec,1)+sum(max(Edat,0));
        statics.morph.camber.Eneg(Fsec,1) = statics.morph.camber.Eneg(Fsec,1)+sum(min(Edat,0));
        indpos = find(max(Edat,0));
        indneg = find(min(Edat,0));
        
         
%         statics.morph.camber.E(Fsec,1) = -sum(Fa.*dr,1)';
%         if Fsec<=6
%             for plotind=1:length(Fa)
%                 figure(Fsec)
%                 subplot(2,1,1)
%                 hold all
%                 plot((statics.morph.camber.xlocstat{1}(Frow,1)+statics.morph.camber.xlocstat{1}(Frow+1,1))/2*ones(2,1)+[0;drx(plotind)],[0;drz(plotind)],'x')
%                 subplot(2,1,2)
%                 hold all
%                 plot((statics.morph.camber.xlocstat{1}(Frow,1)+statics.morph.camber.xlocstat{1}(Frow+1,1))/2,Fa(3,plotind)/(statics.morph.camber.xlocstat{1}(Frow+1,1)-statics.morph.camber.xlocstat{1}(Frow,1)),'x')
%             end
%         end
        
        if ders == 1
            for k=1:length(constant.morph.camber.loc)
                dFalocdparamini = reshape(statics.sens.dFalocdparamini(:,k),size(statics.morph.camber.Faloc,2),size(statics.morph.camber.Faloc,1))';
                dFadistdparamini = reshape(statics.sens.dFadistdparamini(:,k),size(statics.morph.camber.Fadist,2),size(statics.morph.camber.Fadist,1))';
                
%                 % Define the undeformed camber motion
%                 drdparamini = [dFalocdparamini(3*(i-1)+1,2:end)-dFalocdparamini(3*(i-1)+1,1:end-1);
%                     zeros(1,size(statics.morph.camber.Faloc,2)-1);
%                     dFalocdparamini(3*(i-1)+3,2:end)-dFalocdparamini(3*(i-1)+3,1:end-1)];
                
                dFalocmatxdparamini = (dFalocdparamini(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+1:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+1,2:end)-...
                    dFalocdparamini(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+1:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+1,1:end-1));
                
                drxdparamini = sum((statics.morph.camber.xloccam{1}(rowstart+1:rowend,ones(nsteps,1))-statics.morph.camber.xloccam{1}(rowstart:rowend-1,ones(nsteps,1))).*...
                    (dFalocmatxdparamini(1:end-1,:)+dFalocmatxdparamini(2:end,:))/2)/...
                    (statics.morph.camber.xloccam{1}(rowend,1)-statics.morph.camber.xloccam{1}(rowstart,1));
                
                dFalocmatzdparamini = (dFalocdparamini(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+3:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+3,2:end)-...
                    dFalocdparamini(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+3:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+3,1:end-1));
                
                drzdparamini = sum((statics.morph.camber.xloccam{1}(rowstart+1:rowend,ones(nsteps,1))-statics.morph.camber.xloccam{1}(rowstart:rowend-1,ones(nsteps,1))).*...
                    (dFalocmatzdparamini(1:end-1,:)+dFalocmatzdparamini(2:end,:))/2)/...
                    (statics.morph.camber.xloccam{1}(rowend,1)-statics.morph.camber.xloccam{1}(rowstart,1));
                
                
                % Define the undeformed camber motion
                drdparamini = [drxdparamini;
                    zeros(1,nsteps);
                    drzdparamini];
                
                for j = 1:nstep
                    dFarotdparamini = reshape(statics.sens.dFarotdparamini(j:nstep:end,k),size(statics.morph.camber.Farot{j},2),size(statics.morph.camber.Farot{j},1))';
                    
                    % Rotate dr with the section angle
                    Rdir = statics.morph.camber.Farot{j}(3*(Fsec-1)+(1:3),:);
                    dRdirdparamini = dFarotdparamini(3*(Fsec-1)+(1:3),:);
                    dFadparamini(:,j) = dRdirdparamini'*statics.morph.camber.Fadist(3*(i-1)+(1:3),j)+Rdir'*dFadistdparamini(3*(i-1)+(1:3),j);
                end
                % Average the load of a morphing step
                dFadparamini = (dFadparamini(:,2:end)+dFadparamini(:,1:end-1))/2;
                % take care that the work that needs to be delivered is the
                % negative of the aerodynamic work, since sum(W) = 0
                dEdat = -sum(Fa.*drdparamini+dFadparamini.*dr);
                statics.morph.camber.dEdparamini(Fsec,k) = statics.morph.camber.dEdparamini(Fsec,k)-sum(sum(Fa.*drdparamini+dFadparamini.*dr));
                statics.morph.camber.dEposdparamini(Fsec,k) = statics.morph.camber.dEposdparamini(Fsec,k)+sum(dEdat(indpos));
                statics.morph.camber.dEnegdparamini(Fsec,k) = statics.morph.camber.dEnegdparamini(Fsec,k)+sum(dEdat(indneg));
                
                dFalocdparamfin = reshape(statics.sens.dFalocdparamfin(:,k),size(statics.morph.camber.Faloc,2),size(statics.morph.camber.Faloc,1))';
                dFadistdparamfin = reshape(statics.sens.dFadistdparamfin(:,k),size(statics.morph.camber.Fadist,2),size(statics.morph.camber.Fadist,1))';
                
                % Define the undeformed camber motion
%                 drdparamfin = [dFalocdparamfin(3*(i-1)+1,2:end)-dFalocdparamfin(3*(i-1)+1,1:end-1);
%                     zeros(1,size(statics.morph.camber.Faloc,2)-1);
%                     dFalocdparamfin(3*(i-1)+3,2:end)-dFalocdparamfin(3*(i-1)+3,1:end-1)];
                
                dFalocmatxdparamfin = (dFalocdparamfin(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+1:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+1,2:end)-...
                    dFalocdparamfin(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+1:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+1,1:end-1));
                
                drxdparamfin = sum((statics.morph.camber.xloccam{1}(rowstart+1:rowend,ones(nsteps,1))-statics.morph.camber.xloccam{1}(rowstart:rowend-1,ones(nsteps,1))).*...
                    (dFalocmatxdparamfin(1:end-1,:)+dFalocmatxdparamfin(2:end,:))/2)/...
                    (statics.morph.camber.xloccam{1}(rowend,1)-statics.morph.camber.xloccam{1}(rowstart,1));
                
                dFalocmatzdparamfin = (dFalocdparamfin(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+3:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+3,2:end)-...
                    dFalocdparamfin(3*((rowstart-1)*constant.aero.nbs+(Fsec-1))+3:3*constant.aero.nbs:3*((rowend-1)*constant.aero.nbs+(Fsec-1))+3,1:end-1));
                
                drzdparamfin = sum((statics.morph.camber.xloccam{1}(rowstart+1:rowend,ones(nsteps,1))-statics.morph.camber.xloccam{1}(rowstart:rowend-1,ones(nsteps,1))).*...
                    (dFalocmatzdparamfin(1:end-1,:)+dFalocmatzdparamfin(2:end,:))/2)/...
                    (statics.morph.camber.xloccam{1}(rowend,1)-statics.morph.camber.xloccam{1}(rowstart,1));
                
                
                % Define the undeformed camber motion
                drdparamfin = [drxdparamfin;
                    zeros(1,nsteps);
                    drzdparamfin];
                
                for j = 1:nstep
                    dFarotdparamfin = reshape(statics.sens.dFarotdparamfin(j:nstep:end,k),size(statics.morph.camber.Farot{j},2),size(statics.morph.camber.Farot{j},1))';
                    
                    % Rotate dr with the section angle
                    Rdir = statics.morph.camber.Farot{j}(3*(Fsec-1)+(1:3),:);
                    dRdirdparamfin = dFarotdparamfin(3*(Fsec-1)+(1:3),:);
                    dFadparamfin(:,j) = dRdirdparamfin'*statics.morph.camber.Fadist(3*(i-1)+(1:3),j)+Rdir'*dFadistdparamfin(3*(i-1)+(1:3),j);
                end
                % Average the load of a morphing step
                dFadparamfin = (dFadparamfin(:,2:end)+dFadparamfin(:,1:end-1))/2;
                % take care that the work that needs to be delivered is the
                % negative of the aerodynamic work, since sum(W) = 0
                dEdat = -sum(Fa.*drdparamfin+dFadparamfin.*dr);
                statics.morph.camber.dEdparamfin(Fsec,k) = statics.morph.camber.dEdparamfin(Fsec,k)-sum(sum(Fa.*drdparamfin+dFadparamfin.*dr));
                
                statics.morph.camber.dEposdparamfin(Fsec,k) = statics.morph.camber.dEposdparamfin(Fsec,k)+sum(dEdat(indpos));
                statics.morph.camber.dEnegdparamfin(Fsec,k) = statics.morph.camber.dEnegdparamfin(Fsec,k)+sum(dEdat(indneg));
            end
            if isfield(statics.morph,'shear')
                for k=1:constant.str.Ns
                    dFadistdpsiini = reshape(statics.sens.dFadistdpsiini(:,k),size(statics.morph.camber.Fadist,2),size(statics.morph.camber.Fadist,1))';
                    
                    for j = 1:nstep
                        dFarotdpsiini = reshape(statics.sens.dFarotdpsiini(j:nstep:end,k),size(statics.morph.camber.Farot{j},2),size(statics.morph.camber.Farot{j},1))';
                        
                        % Rotate dr with the section angle
                        Rdir = statics.morph.camber.Farot{j}(3*(Fsec-1)+(1:3),:);
                        dRdirdpsiini = dFarotdpsiini(3*(Fsec-1)+(1:3),:);
                        dFadpsiini(:,j) = dRdirdpsiini'*statics.morph.camber.Fadist(3*(i-1)+(1:3),j)+Rdir'*dFadistdpsiini(3*(i-1)+(1:3),j);
                    end
                    % Average the load of a morphing step
                    dFadpsiini = (dFadpsiini(:,2:end)+dFadpsiini(:,1:end-1))/2;
                    % take care that the work that needs to be delivered is the
                    % negative of the aerodynamic work, since sum(W) = 0
                    dEdat = -sum(dFadpsiini.*dr);
                    statics.morph.camber.dEdpsiini(Fsec,k) = statics.morph.camber.dEdpsiini(Fsec,k)-sum(sum(dFadpsiini.*dr));
                    statics.morph.camber.dEposdpsiini(Fsec,k) = statics.morph.camber.dEposdpsiini(Fsec,k)+sum(dEdat(indpos));
                    statics.morph.camber.dEnegdpsiini(Fsec,k) = statics.morph.camber.dEnegdpsiini(Fsec,k)+sum(dEdat(indneg));
                    
                    dFadistdpsifin = reshape(statics.sens.dFadistdpsifin(:,k),size(statics.morph.camber.Fadist,2),size(statics.morph.camber.Fadist,1))';
                    
                    for j = 1:nstep
                        dFarotdpsifin = reshape(statics.sens.dFarotdpsifin(j:nstep:end,k),size(statics.morph.camber.Farot{j},2),size(statics.morph.camber.Farot{j},1))';
                        
                        % Rotate dr with the section angle
                        Rdir = statics.morph.camber.Farot{j}(3*(Fsec-1)+(1:3),:);
                        dRdirdpsifin = dFarotdpsifin(3*(Fsec-1)+(1:3),:);
                        dFadpsifin(:,j) = dRdirdpsifin'*statics.morph.camber.Fadist(3*(i-1)+(1:3),j)+Rdir'*dFadistdpsifin(3*(i-1)+(1:3),j);
                    end
                    % Average the load of a morphing step
                    dFadpsifin = (dFadpsifin(:,2:end)+dFadpsifin(:,1:end-1))/2;
                    % take care that the work that needs to be delivered is the
                    % negative of the aerodynamic work, since sum(W) = 0
                    dEdat = -sum(dFadpsifin.*dr);
                    statics.morph.camber.dEdpsifin(Fsec,k) = statics.morph.camber.dEdpsifin(Fsec,k)-sum(sum(dFadpsifin.*dr));
                    statics.morph.camber.dEposdpsifin(Fsec,k) = statics.morph.camber.dEposdpsifin(Fsec,k)+sum(dEdat(indpos));
                    statics.morph.camber.dEnegdpsifin(Fsec,k) = statics.morph.camber.dEnegdpsifin(Fsec,k)+sum(dEdat(indneg));
                end
            end
            if isfield(statics.morph,'fold')
                for k=1:constant.str.Ns
                    dFadistdthetaini = reshape(statics.sens.dFadistdthetaini(:,k),size(statics.morph.camber.Fadist,2),size(statics.morph.camber.Fadist,1))';
                    
                    for j = 1:nstep
                        dFarotdthetaini = reshape(statics.sens.dFarotdthetaini(j:nstep:end,k),size(statics.morph.camber.Farot{j},2),size(statics.morph.camber.Farot{j},1))';
                        
                        % Rotate dr with the section angle
                        Rdir = statics.morph.camber.Farot{j}(3*(Fsec-1)+(1:3),:);
                        dRdirdthetaini = dFarotdthetaini(3*(Fsec-1)+(1:3),:);
                        dFadthetaini(:,j) = dRdirdthetaini'*statics.morph.camber.Fadist(3*(i-1)+(1:3),j)+Rdir'*dFadistdthetaini(3*(i-1)+(1:3),j);
                    end
                    % Average the load of a morphing step
                    dFadthetaini = (dFadthetaini(:,2:end)+dFadthetaini(:,1:end-1))/2;
                    % take care that the work that needs to be delivered is the
                    % negative of the aerodynamic work, since sum(W) = 0
                    dEdat = -sum(dFadthetaini.*dr);
                    statics.morph.camber.dEdthetaini(Fsec,k) = statics.morph.camber.dEdthetaini(Fsec,k)-sum(sum(dFadthetaini.*dr));
                    statics.morph.camber.dEposdthetaini(Fsec,k) = statics.morph.camber.dEposdthetaini(Fsec,k)+sum(dEdat(indpos));
                    statics.morph.camber.dEnegdthetaini(Fsec,k) = statics.morph.camber.dEnegdthetaini(Fsec,k)+sum(dEdat(indneg));
                    
                    dFadistdthetafin = reshape(statics.sens.dFadistdthetafin(:,k),size(statics.morph.camber.Fadist,2),size(statics.morph.camber.Fadist,1))';
                    
                    for j = 1:nstep
                        dFarotdthetafin = reshape(statics.sens.dFarotdthetafin(j:nstep:end,k),size(statics.morph.camber.Farot{j},2),size(statics.morph.camber.Farot{j},1))';
                        
                        % Rotate dr with the section angle
                        Rdir = statics.morph.camber.Farot{j}(3*(Fsec-1)+(1:3),:);
                        dRdirdthetafin = dFarotdthetafin(3*(Fsec-1)+(1:3),:);
                        dFadthetafin(:,j) = dRdirdthetafin'*statics.morph.camber.Fadist(3*(i-1)+(1:3),j)+Rdir'*dFadistdthetafin(3*(i-1)+(1:3),j);
                    end
                    % Average the load of a morphing step
                    dFadthetafin = (dFadthetafin(:,2:end)+dFadthetafin(:,1:end-1))/2;
                    % take care that the work that needs to be delivered is the
                    % negative of the aerodynamic work, since sum(W) = 0
                    dEdat = -sum(dFadthetafin.*dr);
                    statics.morph.camber.dEdthetafin(Fsec,k) = statics.morph.camber.dEdthetafin(Fsec,k)-sum(sum(dFadthetafin.*dr));
                    statics.morph.camber.dEposdthetafin(Fsec,k) = statics.morph.camber.dEposdthetafin(Fsec,k)+sum(dEdat(indpos));
                    statics.morph.camber.dEnegdthetafin(Fsec,k) = statics.morph.camber.dEnegdthetafin(Fsec,k)+sum(dEdat(indneg));
                    
                end
            end
            if isfield(statics.morph,'span')
                for k=1:constant.str.Ns
                    dFadistdextini = reshape(statics.sens.dFadistdextini(:,k),size(statics.morph.camber.Fadist,2),size(statics.morph.camber.Fadist,1))';
                    
                    for j = 1:nstep
                        dFarotdextini = reshape(statics.sens.dFarotdextini(j:nstep:end,k),size(statics.morph.camber.Farot{j},2),size(statics.morph.camber.Farot{j},1))';
                        
                        % Rotate dr with the section angle
                        Rdir = statics.morph.camber.Farot{j}(3*(Fsec-1)+(1:3),:);
                        dRdirdextini = dFarotdextini(3*(Fsec-1)+(1:3),:);
                        dFadextini(:,j) = dRdirdextini'*statics.morph.camber.Fadist(3*(i-1)+(1:3),j)+Rdir'*dFadistdextini(3*(i-1)+(1:3),j);
                    end
                    % Average the load of a morphing step
                    dFadextini = (dFadextini(:,2:end)+dFadextini(:,1:end-1))/2;
                    % take care that the work that needs to be delivered is the
                    % negative of the aerodynamic work, since sum(W) = 0
                    dEdat = -sum(dFadextini.*dr);
                    statics.morph.camber.dEdextini(Fsec,k) = statics.morph.camber.dEdextini(Fsec,k)-sum(sum(dFadextini.*dr));
                    statics.morph.camber.dEposdextini(Fsec,k) = statics.morph.camber.dEposdextini(Fsec,k)+sum(dEdat(indpos));
                    statics.morph.camber.dEnegdextini(Fsec,k) = statics.morph.camber.dEnegdextini(Fsec,k)+sum(dEdat(indneg));
                    
                    dFadistdextfin = reshape(statics.sens.dFadistdextfin(:,k),size(statics.morph.camber.Fadist,2),size(statics.morph.camber.Fadist,1))';
                    
                    for j = 1:nstep
                        dFarotdextfin = reshape(statics.sens.dFarotdextfin(j:nstep:end,k),size(statics.morph.camber.Farot{j},2),size(statics.morph.camber.Farot{j},1))';
                        
                        % Rotate dr with the section angle
                        Rdir = statics.morph.camber.Farot{j}(3*(Fsec-1)+(1:3),:);
                        dRdirdextfin = dFarotdextfin(3*(Fsec-1)+(1:3),:);
                        dFadextfin(:,j) = dRdirdextfin'*statics.morph.camber.Fadist(3*(i-1)+(1:3),j)+Rdir'*dFadistdextfin(3*(i-1)+(1:3),j);
                    end
                    % Average the load of a morphing step
                    dFadextfin = (dFadextfin(:,2:end)+dFadextfin(:,1:end-1))/2;
                    % take care that the work that needs to be delivered is the
                    % negative of the aerodynamic work, since sum(W) = 0
                    dEdat = -sum(dFadextfin.*dr);
                    statics.morph.camber.dEdextfin(Fsec,k) = statics.morph.camber.dEdextfin(Fsec,k)-sum(sum(dFadextfin.*dr));
                    statics.morph.camber.dEposdextfin(Fsec,k) = statics.morph.camber.dEposdextfin(Fsec,k)+sum(dEdat(indpos));
                    statics.morph.camber.dEnegdextfin(Fsec,k) = statics.morph.camber.dEnegdextfin(Fsec,k)+sum(dEdat(indneg));
                    
                end
            end
            if isfield(statics.morph,'twist')
                for k=1:constant.str.Ns
                    dFadistdphiini = reshape(statics.sens.dFadistdphiini(:,k),size(statics.morph.camber.Fadist,2),size(statics.morph.camber.Fadist,1))';
                    
                    for j = 1:nstep
                        dFarotdphiini = reshape(statics.sens.dFarotdphiini(j:nstep:end,k),size(statics.morph.camber.Farot{j},2),size(statics.morph.camber.Farot{j},1))';
                        
                        % Rotate dr with the section angle
                        Rdir = statics.morph.camber.Farot{j}(3*(Fsec-1)+(1:3),:);
                        dRdirdphiini = dFarotdphiini(3*(Fsec-1)+(1:3),:);
                        dFadphiini(:,j) = dRdirdphiini'*statics.morph.camber.Fadist(3*(i-1)+(1:3),j)+Rdir'*dFadistdphiini(3*(i-1)+(1:3),j);
                    end
                    % Average the load of a morphing step
                    dFadphiini = (dFadphiini(:,2:end)+dFadphiini(:,1:end-1))/2;
                    % take care that the work that needs to be delivered is the
                    % negative of the aerodynamic work, since sum(W) = 0
                    dEdat = -sum(dFadphiini.*dr);
                    statics.morph.camber.dEdphiini(Fsec,k) = statics.morph.camber.dEdphiini(Fsec,k)-sum(sum(dFadphiini.*dr));
                    statics.morph.camber.dEposdphiini(Fsec,k) = statics.morph.camber.dEposdphiini(Fsec,k)+sum(dEdat(indpos));
                    statics.morph.camber.dEnegdphiini(Fsec,k) = statics.morph.camber.dEnegdphiini(Fsec,k)+sum(dEdat(indneg));
                    
                    dFadistdphifin = reshape(statics.sens.dFadistdphifin(:,k),size(statics.morph.camber.Fadist,2),size(statics.morph.camber.Fadist,1))';
                    
                    for j = 1:nstep
                        dFarotdphifin = reshape(statics.sens.dFarotdphifin(j:nstep:end,k),size(statics.morph.camber.Farot{j},2),size(statics.morph.camber.Farot{j},1))';
                        
                        % Rotate dr with the section angle
                        Rdir = statics.morph.camber.Farot{j}(3*(Fsec-1)+(1:3),:);
                        dRdirdphifin = dFarotdphifin(3*(Fsec-1)+(1:3),:);
                        dFadphifin(:,j) = dRdirdphifin'*statics.morph.camber.Fadist(3*(i-1)+(1:3),j)+Rdir'*dFadistdphifin(3*(i-1)+(1:3),j);
                    end
                    % Average the load of a morphing step
                    dFadphifin = (dFadphifin(:,2:end)+dFadphifin(:,1:end-1))/2;
                    % take care that the work that needs to be delivered is the
                    % negative of the aerodynamic work, since sum(W) = 0
                    dEdat = -sum(dFadphifin.*dr);
                    statics.morph.camber.dEdphifin(Fsec,k) = statics.morph.camber.dEdphifin(Fsec,k)-sum(sum(dFadphifin.*dr));
                    statics.morph.camber.dEposdphifin(Fsec,k) = statics.morph.camber.dEposdphifin(Fsec,k)+sum(dEdat(indpos));
                    statics.morph.camber.dEnegdphifin(Fsec,k) = statics.morph.camber.dEnegdphifin(Fsec,k)+sum(dEdat(indneg));
                    
                end
            end
            if tailflag == 1
                for k=1:numel(statics.str.C)
                    dFadistdC = reshape(statics.sens.dFadistdC(:,k),size(statics.morph.camber.Fadist,2),size(statics.morph.camber.Fadist,1))';
                    
                    for j = 1:nstep
                        dFarotdC = reshape(statics.sens.dFarotdC(j:nstep:end,k),size(statics.morph.camber.Farot{j},2),size(statics.morph.camber.Farot{j},1))';
                        
                        % Rotate dr with the section angle
                        Rdir = statics.morph.camber.Farot{j}(3*(Fsec-1)+(1:3),:);
                        dRdirdC = dFarotdC(3*(Fsec-1)+(1:3),:);
                        dFadC(:,j) = dRdirdC'*statics.morph.camber.Fadist(3*(i-1)+(1:3),j)+Rdir'*dFadistdC(3*(i-1)+(1:3),j);
                    end
                    % Average the load of a morphing step
                    dFadC = (dFadC(:,2:end)+dFadC(:,1:end-1))/2;
                    % take care that the work that needs to be delivered is the
                    % negative of the aerodynamic work, since sum(W) = 0
                    dEdat = -sum(dFadC.*dr);
                    statics.morph.camber.dEdC(Fsec,k) = statics.morph.camber.dEdC(Fsec,k)-sum(sum(dFadC.*dr));
                    statics.morph.camber.dEposdC(Fsec,k) = statics.morph.camber.dEdC(Fsec,k)+sum(dEdat(indpos));
                    statics.morph.camber.dEdC(Fsec,k) = statics.morph.camber.dEdC(Fsec,k)+sum(dEdat(indneg));
                end
                for k=1:length(constant.lam.ID)
                    dFadistdt = reshape(statics.sens.dFadistdt(:,k),size(statics.morph.camber.Fadist,2),size(statics.morph.camber.Fadist,1))';
                    
                    for j = 1:nstep
                        dFarotdt = reshape(statics.sens.dFarotdt(j:nstep:end,k),size(statics.morph.camber.Farot{j},2),size(statics.morph.camber.Farot{j},1))';
                        
                        % Rotate dr with the section angle
                        Rdir = statics.morph.camber.Farot{j}(3*(Fsec-1)+(1:3),:);
                        dRdirdt = dFarotdt(3*(Fsec-1)+(1:3),:);
                        dFadt(:,j) = dRdirdt'*statics.morph.camber.Fadist(3*(i-1)+(1:3),j)+Rdir'*dFadistdt(3*(i-1)+(1:3),j);
                    end
                    % Average the load of a morphing step
                    dFadt = (dFadt(:,2:end)+dFadt(:,1:end-1))/2;
                    % take care that the work that needs to be delivered is the
                    % negative of the aerodynamic work, since sum(W) = 0
                    dEdat = -sum(dFadt.*dr);
                    statics.morph.camber.dEdt(Fsec,k) = statics.morph.camber.dEdt(Fsec,k)-sum(sum(dFadt.*dr));
                    statics.morph.camber.dEposdt(Fsec,k) = statics.morph.camber.dEposdt(Fsec,k)+sum(dEdat(indpos));
                    statics.morph.camber.dEnegdt(Fsec,k) = statics.morph.camber.dEnegdt(Fsec,k)+sum(dEdat(indneg));
                end
            end
        end
%         % Define the undeformed camber motion
%         dr = [statics.morph.camber.Faloc(3*(i-1)+1,2:end)-statics.morph.camber.Faloc(3*(i-1)+1,1:end-1);
%             zeros(1,size(statics.morph.camber.Faloc,2)-1);
%             statics.morph.camber.Faloc(3*(i-1)+3,2:end)-statics.morph.camber.Faloc(3*(i-1)+3,1:end-1);];
%         % Find the section of the aerodynamic force
%         Fsec = rem((i-1),constant.aero.nbs)+1;
%         for j = 1:nstep-1
%             % Rotate dr with the section angle
%             Rdir = statics.morph.camber.Farot{j}(3*(Fsec-1)+(1:3),:);
%             Rdir2 = statics.morph.camber.Farot{j+1}(3*(Fsec-1)+(1:3),:);
%             dr1(:,j) = Rdir*dr(:,j);
%             dr2(:,j) = Rdir2*dr(:,j);
%         end
%         
%         dr = (dr1+dr2)/2;
%         % Average the load of a morphing step 
%         Fa = (statics.morph.camber.Fadist(3*(i-1)+(1:3),2:end)+statics.morph.camber.Fadist(3*(i-1)+(1:3),1:end-1))/2
%         % take care that the work that needs to be delivered is the
%         % negative of the aerodynamic work, since sum(W) = 0
%         statics.morph.camber.E2(Fsec,1) = -sum(sum(Fa.*dr));
    end
end
%}