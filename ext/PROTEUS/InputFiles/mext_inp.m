
function lumped = mext_inp (str,lumped)

if isfield(lumped,'mass')
    
    nr_massfields=length(lumped.mass);
    
    for i=1:nr_massfields;
        
        nr_masses= length(lumped.mass{i});
        
        for j=1:nr_masses;
            
            location=lumped.location{i}(j,:);
            xyz=reshape(str.xyz,3,str.Ns+1);
            
            
            for k=1:str.Ns;
                d(k)=norm(location-xyz(:,k)');
            end
            
            [c index]=min(abs(d));
            
            if (index==str.Ns+1)
                element=index-1;
            else
                element=index;
            end
            
            
            ecc_local=projection_lumped(str,element,location);
            if (ecc_local(1)<0)
                element=element-1;
            end
            
            if (element==0)
                element=1;
            end
            
            ecc_local=projection_lumped(str,element,location);
            
            
            
            lumped.ecc{i}(j,1:3)=ecc_local';
            lumped.element{i}(j)=element;
            
            
            
        end
        
    end
end

%

function [ecc_local] = projection_lumped(str,element,location)

v1=location'-str.xyz(3*(element-1)+(1:3));
R0 = str.R0((element-1)*3+(1:3),:);
ecc_local=R0'*v1;
L_el=norm(str.xyz(3*(element)+(1:3))-str.xyz(3*(element-1)+(1:3)));
ecc_local(1)=ecc_local(1)/L_el;



