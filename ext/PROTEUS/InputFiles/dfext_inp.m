function fext = dfext_inp(str,fext,ders)

nr_fextfields=length(fext.magnitude);


for i=1:nr_fextfields
    
 
   
  nr_fext=size(fext.magnitude{i},1);
  
  for j=1:nr_fext
      
       fextloc=fext.location{i}(j,:);
       xyz=reshape(str.xyz,3,str.Ns+1);
      

%determine the element which is closest to the external force 
               for k=1:str.Ns+1;
                    d(k)=norm(fextloc-xyz(:,k)');
               end
                
               
                [c index]=min(abs(d));


                    if (index==str.Ns+1)
                        element=index-1;
                    else 
                        element=index;
                    end

                if (projection(fextloc,str,element)<0)
                    element=element-1;
                end

                if(element==0)
                    element=1;
                end


                xi=projection(fextloc,str,element);

                node1=element;
                node2=element+1;

                dof1=str.EFT(element,1:6);
                dof2=str.EFT(element,7:12);

                xa_undeformed=(1-xi)*(str.xyz((node1-1)*3+(1:3)))+xi*(str.xyz((node2-1)*3+(1:3)));


                %rigid link, vector connecting position with artificial undeformed node;
                v0=fextloc-xa_undeformed';
                if ders == 1 && isempty(fext.dlocationdt{i}) == 0
                    fext.dv0dt{i}(3*(j-1)+(1:3),:) = fext.dlocationdt{i}(3*(j-1)+(1:3),:);
                end
                fext.dof1{i}(j,:)=dof1;
                fext.dof2{i}(j,:)=dof2;
                fext.xi{i}(j,:)=xi;
                fext.v0{i}(j,:)=v0;
                fext.element{i}(j,:)=element;
  end
                fext.nr_fextfields=nr_fextfields;
                fext.nr_fext{i}=nr_fext;
end



function [xi] = projection(fextloc,str,element)

v1=fextloc'-str.xyz(3*(element-1)+(1:3));
Ro = str.R0((element-1)*3+(1:3),:);
ecc_local=Ro'*v1;
L_el=norm(str.xyz(3*(element)+(1:3))-str.xyz(3*(element-1)+(1:3)));
xi=ecc_local(1)/L_el;
