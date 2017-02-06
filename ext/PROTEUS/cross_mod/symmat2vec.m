function [out] = symmat2vec(in,typ,d)
%this function converts symmetric matrix into a vector and vise versa
%input: 
%   in: could be a syymetric matrix or vector of a symmetric matrix
%   typ: type of output wanted: 1 if a vector is wanted, -1 if a matrix is  wanted
%   d: size of the symmetric matrix
%output
%   out: a vector or matrix depending on the value for typ

switch typ
    case 1 % matrix to vector
        nd = size(in,1);
        N = 1/2*nd*(nd-1); %number of elements
        out = zeros(N,1);     
        en = 0;
        for i=1:size(in,1)
            vc = diag(in,i-1);
            st = en+1;
            en = en+length(vc);
            out(st:en) = vc;
        end
    case -1 %vector to matrix
        out = zeros(d,d);
        en = 0;
        for i=1:d
            st = en+1;
            en = en+d-(i-1);
            dvec = in(st:en);
            out = out + diag(dvec,i-1);
        end
        out = (out'+out)-diag(in(1:d));
end

end


% % switch typ
% %     case 1 % matrix to vector
% %         out = in(find(triu(in))); %#ok<FNDSB> extract the unique elements of the symetric matrix
% %     case -1 %vector to matrix
% %         sm = 1/2*(sqrt(8*length(in)+1)-1); %size of the matrix
% %         out = zeros(sm); %initialization
% %         out(find(triu(ones(sm)))) = in; %#ok<FNDSB> fill in the upper diagonal
% %         out = out+triu(out,1)'; %fill in the lower diagonal        
% % end