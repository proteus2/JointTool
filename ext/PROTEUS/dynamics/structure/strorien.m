function [str] = strorien(str,sens,lin)
strxyz = reshape(str.xyz0,3,[])';
e1elm = reshape(str.e1elm,3,[]);
e2elm = reshape(str.e2elm,3,[]);
e3elm = reshape(str.e3elm,3,[]);


for i=1:str.Nel
   str.elm.ell(i)=norm(strxyz(i+1,:)-strxyz(i,:));
   if sens == 1
       if lin ~= 1
           dstrxyz2 = sparse(3,numel(str.xyz));
           dstrxyz2(:,1+i*3:3+i*3) = eye(3);
           dstrxyz1 = sparse(3,numel(str.xyz));
           dstrxyz1(:,1+(i-1)*3:3+(i-1)*3) = eye(3);
           str.elm.delldxyz0(i,:) = dnorm((strxyz(i+1,:)-strxyz(i,:))',(dstrxyz2-dstrxyz1));
       else
           str.elm.delldxyz0(i,:) = sparse(1,numel(str.xyz));
       end
   end
end

Relm = permute(cat(3, e1elm, e2elm, e3elm),[2 1 3]);

% convert to local element coordinates
str.elm.R = zeros(str.Nel,12,12);
str.elm.R(:,1:3,1:3)     = Relm;
str.elm.R(:,4:6,4:6)     = Relm;
str.elm.R(:,7:9,7:9)     = Relm;
str.elm.R(:,10:12,10:12) = Relm;

if sens == 1
    dRde1loc = zeros(9,3);
    dRde1loc(1,1) = 1;
    dRde1loc(4,2) = 1;
    dRde1loc(7,3) = 1;
    dRde1 = sparse(9*str.Nel,3*str.Nel);
    dRde2loc = zeros(9,3);
    dRde2loc(2,1) = 1;
    dRde2loc(5,2) = 1;
    dRde2loc(8,3) = 1;
    dRde3 = sparse(9*str.Nel,3*str.Nel);
    dRde3loc = zeros(9,3);
    dRde3loc(3,1) = 1;
    dRde3loc(6,2) = 1;
    dRde3loc(9,3) = 1;
    dRde1 = sparse(9*str.Nel,3*str.Nel);
    for i=1:str.Nel
        dRde1((i-1)*9+(1:9),(i-1)*3+(1:3)) = dRde1loc;
        dRde2((i-1)*9+(1:9),(i-1)*3+(1:3)) = dRde2loc;
        dRde3((i-1)*9+(1:9),(i-1)*3+(1:3)) = dRde3loc;
    end
    
    str.dRde1 = dRde1;
    str.dRde2 = dRde2;
    str.dRde3 = dRde3;
end

function [dnorma] = dnorm(a,da)

dnorma = 1/2/norm(a)*(sum(2*a(:,ones(1,size(da,2))).*da));





