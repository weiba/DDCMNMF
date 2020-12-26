function [W,H1,H2,TerminalObj,iter] = SNMNMF_algorithm(X1,X2,A11,A12,A22,params)
%
% INPUT:
%
% X1 : (m x n1) non negative input matrix.
% X2 : (m x n2) non negative input matrix.
% A11 : (n1 x n1) adjacent matrix about the features in X1.
% A12 : (n1 x n2) adjacent matrix about the features in X1 and X2.
% A22 : (n2 x n2) adjacent matrix about the features in X2.
% params: a structure variable, containing all the parameters in this
%         algorithm, including the below variables.
% params.thrNet11 : A parameter about the constraint for network A11.
% params.thrNet12 : about the constraint for network A12.
% params.thrNet22 : about the constraint for network A22.
% params.thrXr : threshold to limit the growth of W.
% params.thrXc : threshold to limit the growth of H.      
% params.K : the number of modules it will identify.
% params.maxiter : the number of maximal iteration step for algorithm
% params.tol : the precision for convergence of algorithm 
%
% OUTPUT:
%
% W : (m x K) basis matrix.
% H1 : (K x n1) weight matrix.
% H2 : (K x n2) weight matrix.
% TerminalObj: (iter x 8) matrix to record each term of the objective  
%              function and the last one is the value of objective function.
% iter: The number of iterations when it stops.
%
if nargin == 0
    help SNMNMF_algorithm
    return
end

%%%%%%%%%%%%%%%%%%%%%%%% Test for non-negative %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for negative values in X1 and X2
if (min(min(X1)) < 0) || (min(min(X2)) < 0)
    error('Input matrix elements can not be negative');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for same rows in X1 X2
[rX1,cX1] = size(X1);
[rX2,cX2] = size(X2);
if (rX1 ~= rX2)
    error('Input matrices should have the same rows');
end

rX = rX1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = params.K;
maxiter = params.maxiter;
tol = params.tol;

thrNet11 = params.thrNet11; 
thrNet12 = params.thrNet12; 
thrNet22 = params.thrNet22; 
thrXr = params.thrXr;
thrXc = params.thrXc; 

% Initialize W, H1 and H2 with random values.
W = rand(rX,K);
H1 = rand(K,cX1);
H2 = rand(K,cX2);
obj = zeros(maxiter,8);

old_obj = Inf;

ResultsFile = 'SNMNMF_Results';
%file_path1='W/'
for iter = 1:maxiter
    % Multiplicative update method
    W = W.*([H1 H2]*[X1 X2]')'./(W*([H1 H2]*[H1 H2]'+ thrXr*eye(K)) + eps);
    %str=[file_path1,num2str(iter),'.txt']
    %if(iter<=100)
    %save(str,'W','-ascii')
    %end
    HH1 = H1.*(W'*X1 + thrNet11*H1*A11 + 0.5*thrNet12*H2*A12')./((W'*W + thrXc*ones(K))*H1 + eps);
    H2 = H2.*(W'*X2 + thrNet22*H2*A22 + 0.5*thrNet12*H1*A12)./((W'*W + thrXc*ones(K))*H2 + eps);
    H1 = HH1;
    clear HH1
    res1 = X1 - W*H1;
    res2 = X2 - W*H2;
    obj(iter,1:2) = [sum(sum((res1).^2)),sum(sum((res2).^2))];
    obj(iter,3:5) = [-thrNet11*trace(H1*A11*H1'),-thrNet12*trace(H1*A12*H2'),-thrNet22*trace(H2*A22*H2')]; 
    obj(iter,6:7) = [thrXr*sum(sum(W.*W)),thrXc*(sum(sum(H1.*H1))+sum(sum(H2.*H2)))];
    obj(iter,8) = sum(obj(iter,1:7));

    if(isnan(obj(iter,1)) || isnan(obj(iter,2)))
        TerminalObj = obj(1:iter,:);
        save([ResultsFile '/' ResultsFile '.mat'],'W','H1','H2','TerminalObj')
        error('Please check the results, and then alter your parameters! ')
    end
    new_obj = obj(iter,8);
    abs(new_obj - old_obj)/new_obj
    if(abs(new_obj - old_obj) < tol)
        TerminalObj = obj(1:iter,:);
        return;
    else
        old_obj = new_obj;
    end
    
end
TerminalObj = obj(1:iter,:);
end

