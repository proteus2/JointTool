function [cross,elmconv] = prepr_cross(cross)

Nel = length(cross.yzlocal);

for j=1:Nel
    for k=1:length(cross.yzlocal{j})
        [yz_discr,elmloc_discr,lamflag,elmconv{j}{k},type] = discret(cross.yzlocal{j}{k},cross.elmloc{j}{k},cross.numel,cross.lam{j}{k},cross.type{j}{k}); % Discretize the cross-section into at least numel elements
        cross.yzlocal{j}{k} = yz_discr;
        cross.elmloc{j}{k} = elmloc_discr;
        cross.lam{j}{k} = lamflag;
        cross.elmID{j}{k} = (1:size(elmloc_discr,1))';
        cross.type{j}{k} = type;
    end
end