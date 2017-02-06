function [r,Psi,s,a] = quadbuckl(N,D,nmodes,stiff)

K11 = stiff.K11;
K12 = stiff.K12;
K16 = stiff.K16;
K22 = stiff.K22;
K26 = stiff.K26;
K66 = stiff.K66;
Kxx = stiff.Kxx;
Kxy = stiff.Kxy;
Kyy = stiff.Kyy;

s = zeros(3,nmodes);
% Phi = zeros(3,3,nmodes);
Psi = zeros(3,3,nmodes);
Kg = N(1)*Kxx+N(2)*Kyy+N(3)*Kxy;
K = D(1,1)*K11+D(1,2)*K12+D(1,3)*K16+D(2,2)*K22+D(2,3)*K26+D(3,3)*K66;

% R = chol(K);

% A = (R'\Kg)/R;

% [a,r] = poseig(A,nmodes);
% for j = 1:nmodes
%     %     ae = a(:,j);
%     ae = R\a(:,j);
%     se = [ae'*Kxx*ae;ae'*Kyy*ae;ae'*Kxy*ae];
%     Psie = -r(j)*[ ...
%         ae'*K11*ae/1,ae'*K12*ae/2,ae'*K16*ae/2; ...
%         ae'*K12*ae/2,ae'*K22*ae/1,ae'*K26*ae/2; ...
%         ae'*K16*ae/2,ae'*K26*ae/2,ae'*K66*ae/1];
%     s(:,j)= se;
% %     Phi(:,:,j) = -D*Psie*D;
%     Psi(:,:,j) = Psie;
% end

[a,r] = poseig(K,Kg,nmodes);
for j = 1:nmodes
    %     ae = a(:,j);
    ae = a(:,j);
    se = [ae'*Kxx*ae;ae'*Kyy*ae;ae'*Kxy*ae];
    Psie = -r(j)*[ ...
        ae'*K11*ae/1,ae'*K12*ae/2,ae'*K16*ae/2; ...
        ae'*K12*ae/2,ae'*K22*ae/1,ae'*K26*ae/2; ...
        ae'*K16*ae/2,ae'*K26*ae/2,ae'*K66*ae/1];
    s(:,j)= se;
%     Phi(:,:,j) = -D*Psie*D;
    Psi(:,:,j) = Psie;
end
end

% function [pVk,pek] = poseig(A,k)
% 
% n = size(A,1);
% [V,D] = eig(A);
% e = diag(D); % eigenvalues
% e(imag(e)~=0) = 0;
% pei = find(e>0); % indices of positive eigenvalues
% pn = length(pei); % number of positive eigenvalues
% pe = e(pei);
% [pe,pes] = sort(-pe); % sorting
% pei = pei(pes);
% 
% if(k<=pn)
%     pek = -pe(1:k);
%     pVk = V(:,pei(1:k));
% else
%     pek(1:pn) = -pe(1:pn);
%     pek(pn+1:k) = 0;
%     pVk = V(:,pei(1:pn));
%     pVk(1:n,pn+1:k) = 0;
% end
% end

function [pVk,pek] = poseig(K,Kg,k)

n = size(K,1);
[V,D] = eig(Kg,K);
e = diag(D); % eigenvalues
pei = find(e>0); % indices of positive eigenvalues
pn = length(pei); % number of positive eigenvalues
pe = e(pei);
[pe,pes] = sort(-pe); % sorting
pei = pei(pes);

if(k<=pn)
    pek = -pe(1:k);
    pVk = V(:,pei(1:k));
else
    pek(1:pn) = -pe(1:pn);
    pek(pn+1:k) = 0;
    pVk = V(:,pei(1:pn));
    pVk(1:n,pn+1:k) = 0;
end
end
