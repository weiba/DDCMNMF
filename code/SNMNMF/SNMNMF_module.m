function module = SNMNMF_module(H,t,isdouble)
%
% Compute the mean and standard deviation of each rows in H to determine 
% the members of module. 
%
% INPUT:
%
% H : A (K x n) matrix, we identify module members based on this matrix.
% t : Threshold.
% isdouble : A binary variable, if isdouble = 1, the matrix H is
%           'double' of the original matrix, when we identify the module
%            member, we should restore that.
%
% OUTPUT:
%
% module : (K x 1) cell array. module{i} contains the feature indexes in 
%          module i. 
% 

[K,m] = size(H);
MH = zeros(K,1); % Records the mean of each row of H.
VH = zeros(K,1); % Records the standard deviation of each row of H.

% Compute MH,VH
for i = 1:K
    h = H(i,:);
    MH(i) = mean(h);
    VH(i) = std(h,0);
end

% Module content.
% Records the indexes of one type of feature in each module.
module = cell(K,1); 

if isdouble
    for i = 1:K
        c = find(H(i,:) >= MH(i) + t*VH(i));
        c(c > m/2) = c(c > m/2) - m/2;  
        % Transform the double (gene) indexes into origin indexes
        c = unique(c);
        module{i,1} = c'; 
    end
else
    for i = 1:K
        c = find(H(i,:) >= MH(i) + t*VH(i));
        module{i,1} = c';
    end
end

clear MH VH H
end
