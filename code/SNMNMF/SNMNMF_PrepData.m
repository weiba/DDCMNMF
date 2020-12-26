function [newInput,isdouble] = SNMNMF_PrepData(Input)
% To ensure the non-negativity of input matrices, we 'double' the original 
% matrix if the matrix has negative elements.
%
% INPUT:
%
% Input: A structure variable, contains 
% Input.data: It should be a (m x n) matrix. Each row represents a sample,
%              and each column is a feature. Input.data = [X1,X2];
% Input.XBlockInd: It should be a (2 x 2) matrix. Each row represents
%                   the column index (start index, end index) of each X
%                   block in the Input.data.
%
% Input.netAdj : It is the adjacent matrix for the network which represents
%                the relationship between all the features in the input 
%                data. A11(n1 x n1),A22(n1 x n2),A12(n1 x n2) are the  
%                adjacency matrices for G1,G2,G3,respectively.
%                Thus, Input.netAdj = [A11 A12; A12' A22].If some network
%                is not available, replace it by zero matrix.
%
% OUTPUT:
%
% newInput: A new structure variable after processing, the components of
%           which are the same as Input.
% isdouble: A binary variable. 0 for no change, 1 for 'double' the original
%           matrix.
%

X = Input.data;
XInd = Input.XBlockInd; 

X1 = X(:,XInd(1,1):XInd(1,2));
X2 = X(:,XInd(2,1):XInd(2,2));
N1 = size(X1,2);
N2 = size(X2,2);
netAdj = Input.netAdj;
A11 = netAdj(1:N1,1:N1);
A12 = netAdj(1:N1,(N1+1):(N1+N2));

if (~isempty(find(X < 0, 1)))
    XX1 = [max(X1,0), max(-X1,0)];
    XX2 = [max(X2,0), max(-X2,0)];
    
    AA11 = [A11, A11;
        A11, A11];
    AA12 = [A12, A12;
        A12, A12];
    isdouble = 1;
    
    newblockInd = [1,2*XInd(1,2);(1+2*XInd(1,2)),2*XInd(2,2)];
else
    XX1 = X1;
    XX2 = X2;
    
    AA11 = A11;
    AA12 = A12;
    isdouble = 0;
    
    newblockInd = XInd;
end

clear X1 X2 A12 A11 Input

newInput.data = [XX1,XX2];
newInput.XBlockInd = newblockInd;
newInput.netAdj = sparse([AA11 AA12; AA12' zeros(size(XX2,2))]);
clear AA12 AA11 XX1 XX2
end

