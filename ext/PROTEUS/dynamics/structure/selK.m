function [dK,dKrs] = selK(dKs,constant)

Nel = constant.str.Ns;
Nnode = constant.str.Ndof/6;

vec1 = [0:6*(Nel+1):6*(Nel+1)*(6*(Nel+1)-1)];
% vec2 = [7:6*(Nel+1)]';
vec2 = [constant.str.frdof]';
% Filter the fixed columns
dK = dKs(reshape(bsxfun(@plus,vec2,vec1),[],1),:);
dKrs = dKs(reshape(bsxfun(@plus,vec2,vec1),[],1),:);

% Filter the fixed rows
vec3 = (constant.str.frdof-1)*6*Nel;
vec4 = (1:6*Nel)';

dK = dK(reshape(bsxfun(@plus,vec4,vec3),[],1),:);

% Select the fixed rows
vec3 = (constant.str.fxdof-1)*6*Nel;
vec4 = (1:6*Nel)';

dKrs = dKrs(reshape(bsxfun(@plus,vec4,vec3),[],1),:);


% dumind = reshape(1:(6*Nel)^2,6*Nel,6*Nel)';
% 
% ind = reshape(dumind',[],1);
% 
% dK = dK(ind,:);
% keyboard