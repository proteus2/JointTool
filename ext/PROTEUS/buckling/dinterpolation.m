function [dyadx,dyady,dyadxa] = dinterpolation(x,y,xa,loc1,loc2)

dyady1 = 1-(xa-x(loc1))./(x(loc2)-x(loc1));
dyady2 = (xa-x(loc1))./(x(loc2)-x(loc1));

dyadxa = (y(loc2)-y(loc1))./(x(loc2)-x(loc1));
dyadxa = sparse(diag(dyadxa));

dyadx1 = (y(loc2)-y(loc1)).*(xa-x(loc2))./(x(loc2)-x(loc1)).^2;
dyadx2 = (y(loc2)-y(loc1)).*(x(loc1)-xa)./(x(loc2)-x(loc1)).^2;

% dyadxa = sparse(length(xa),length(xa));

dyady = sparse(length(xa),length(y));
dyadx = sparse(length(xa),length(x));
for i=1:length(xa)
          dyady(i,loc1(i))=dyady1(i);
          dyadx(i,loc1(i))=dyadx1(i);
          dyady(i,loc2(i))=dyady2(i);
          dyadx(i,loc2(i))=dyadx2(i);
end