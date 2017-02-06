function [lamparvec] = lam2lampar(lam)

lam = [lam; flipud(lam)];

lamparvec = zeros(8,1);
lampar_B  = zeros(4,1);

for i=1:length(lam)
    z1           = -1/2 + (i-1)/length(lam);
    z2           = -1/2 + i/length(lam);
    lamparvec(1) = lamparvec(1) + cosd(2*lam(i))*(z2-z1);
    lamparvec(2) = lamparvec(2) + sind(2*lam(i))*(z2-z1);
    lamparvec(3) = lamparvec(3) + cosd(4*lam(i))*(z2-z1);
    lamparvec(4) = lamparvec(4) + sind(4*lam(i))*(z2-z1);
    lamparvec(5) = lamparvec(5) + 4*cosd(2*lam(i))*(z2^3-z1^3);
    lamparvec(6) = lamparvec(6) + 4*sind(2*lam(i))*(z2^3-z1^3);
    lamparvec(7) = lamparvec(7) + 4*cosd(4*lam(i))*(z2^3-z1^3);
    lamparvec(8) = lamparvec(8) + 4*sind(4*lam(i))*(z2^3-z1^3);
    
    % ---
    if 1    % not used in the code (only for check)
        lampar_B(1)  = lampar_B(1)  + 2*cosd(2*lam(i))*(z2^2-z1^2);
        lampar_B(2)  = lampar_B(2)  + 2*sind(2*lam(i))*(z2^2-z1^2);
        lampar_B(3)  = lampar_B(3)  + 2*cosd(4*lam(i))*(z2^2-z1^2);
        lampar_B(4)  = lampar_B(4)  + 2*sind(4*lam(i))*(z2^2-z1^2);
    end
    % ---
end

if ~isempty(find(abs(lampar_B)>1e-6,1)), display(lampar_B); error('non-symmetric laminate'); end % non-symmetric