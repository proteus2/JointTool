function [xy_discr,elmloc_discr,lamflag,elmconv,type_discr] = discret(xy,elmloc,numel,laminates,type)

circumf = sum(  sqrt(sum((xy(elmloc(:,1),:)-xy(elmloc(:,2),:)).^2,2)) );

elml = circumf/numel;

for i=1:size(elmloc,1)
    if i==1
        numelloc = ceil(norm(xy(elmloc(i,1),:)-xy(elmloc(i,2),:))/elml);
        xy_discr_nu(1:numelloc+1,:)   = [linspace(xy(elmloc(i,1),1),xy(elmloc(i,2),1),numelloc+1)',linspace(xy(elmloc(i,1),2),xy(elmloc(i,2),2),numelloc+1)'];
        elmloc_discr_nu(1:numelloc,:) = [1:numelloc;2:numelloc+1]';
        lamflag(1:numelloc,1)         = laminates(i);
        elmconv(1:numelloc,1) = i;
        type_discr(1:numelloc,1)         = type(i);
    else
        numelloc = ceil(norm(xy(elmloc(i,1),:)-xy(elmloc(i,2),:))/elml);
        elmloc_discr_nu(end+1:end+numelloc,:) = [size(xy_discr_nu,1)+(1:numelloc);size(xy_discr_nu,1)+(2:numelloc+1)]';
        xy_discr_nu(end+1:end+numelloc+1,:) = [linspace(xy(elmloc(i,1),1),xy(elmloc(i,2),1),numelloc+1)',linspace(xy(elmloc(i,1),2),xy(elmloc(i,2),2),numelloc+1)'];
        lamflag(end+1:end+numelloc,1) = laminates(i);
        elmconv(end+1:end+numelloc,1) = i;
        type_discr(end+1:end+numelloc,1) = type(i);
    end
end

[b, m, n] = unique(xy_discr_nu,'rows','first');
msort = sort(m);
xy_discr = xy_discr_nu(msort,:);
for i=1:size(elmloc_discr_nu,1)
   elmxy1  = xy_discr_nu(elmloc_discr_nu(i,1),:);
   elmxy2  = xy_discr_nu(elmloc_discr_nu(i,2),:);
   elmloc1 = find(xy_discr(:,1) == elmxy1(1) & xy_discr(:,2) == elmxy1(2));
   elmloc2 = find(xy_discr(:,1) == elmxy2(1) & xy_discr(:,2) == elmxy2(2));
   elmloc_discr(i,:) = [elmloc1 elmloc2];
end

% figure
% hold all
% plot(xy(:,1),xy(:,2),'+')
% plot(xy_discr(:,1),xy_discr(:,2),'s red')
% keyboard