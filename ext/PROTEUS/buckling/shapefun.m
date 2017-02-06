function [L,L1,L2,w,x] = shapefun(n)

if n < 1
    L = []; L1 = []; L2 = []; w = []; x = [];
    return
end


% Gauss weights and abscissas
T = diag([n:-1:1]./sqrt(2*[n+1:-1:2]-1)./sqrt(2*[n:-1:1]-1),-1); T = T + T';
[V,xg] = eig(T);
w = V(end,:).^2;
x = diag(xg)';
w = 2*w/sum(w);
% Lobatto shape functions and first and second derivatives
L = zeros(n,length(x)); L1 = L; L2 = L;

L2(1,:) = 1;
L1(1,:) = x;
L(1,:) = (x.*L1(1,:)-1)/2;

if n == 1, return; end

L2(2,:) = 3*x;
L1(2,:) = 3*L(1,:)+1;
L(2,:) = (x.*L1(2,:)-L1(1,:))/3;

for k = 3:n
    L1(k,:) = ((2*k-1)*x.*L1(k-1,:)-(k-1)*L1(k-2,:))/k;
    L(k,:) = (x.*L1(k,:)-L1(k-1,:))/(k+1);
    L2(k,:) = k*L1(k-1,:)+x.*L2(k-1,:);
end

end
