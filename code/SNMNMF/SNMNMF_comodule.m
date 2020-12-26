function [W,H1,H2,Comodule,params] = SNMNMF_comodule(Input,params)
%
% INPUT:
%
% Input: A structure variable for input, including
% Input.data: A (m x n) non-negative matrix, combining two data blocks.
% Input.XBlockInd: A (2 x 2) matrix, the column index for each X block.
% Input.netAdj: A (n x n) adjacency matrix for the relationships between
%               the features in Input.data.
%
% params: A structure variables, including
% params.thrNet11: A parameter about the constraint for network A11.
% params.thrNet12: about the constraint for network A12.
% params.thrNet22: about the constraint for network A22.
% params.thrXr: threshold to limit the growth of W.
% params.thrXc: threshold to limit the growth of H.
% params.isdouble: a number to indicate whether the input matrices are
%                  'double'.0 for original and 1 for double.
% params.thrd_module: three thresholds to select elements in two X blocks 
% to construct md-modules. The first threshold is for selecting samples.
% params.K: A pre-defined number of identifiey md-modules.
% params.nloop: The repeating times of SNMNMF algorithm.
% params.maxiter: The maximal number of iterations for this algorithm.
% params.tol: The precision for convergence of algorithm.
%
% OUTPUT:
%
% W : (m x K) basis matrix.
% H1 : (K x n1) weight matrix.
% H2 : (K x n2) weight matrix.
% Comodule : (K x 3) cell, the comodules identified by the SNMNMF algorithm
%            (Comodule{i,j} records the j_th type of variable indexes in 
%            the i_th identified comodule).
% params : adding something new comparing to the input variable --- params.
% (params.thrNet11, params.thrNet12, params.thrXr, params.thrXc,
% params.thrd_module, params.K, params.isdouble,
% params.nloop, params.maxiter, params.tol)
% params.records : (nloop x 1) cell array, each one is a (iter x 8) vector,
%                  recording the values of each term in objective function
%                  and the sum of them.                 
% params.iterNumList : (nloop x 1) vector, each one is the iteration number
%                      in each loop.
%

ResultsFile = 'SNMNMF_Results';
% To make a new folder to store all the results.
if ~isdir(ResultsFile)
    mkdir(ResultsFile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = Input.data;
XInd = Input.XBlockInd;
nbX = size(XInd,1);
%
if (length(params.thrd_module) ~= (1+nbX))
    error(['Error in variable ''thrd_module'': please input ' ...
        num2str(1+nbX) ' values as thresholds for selecting ' ...
        num2str(1+nbX) ' components.'])
end

X1 = X(:,XInd(1,1):XInd(1,2));
X2 = X(:,XInd(2,1):XInd(2,2));
netAdj = Input.netAdj;
A11 = netAdj(XInd(1,1):XInd(1,2),XInd(1,1):XInd(1,2));
A12 = netAdj(XInd(1,1):XInd(1,2),XInd(2,1):XInd(2,2));
A22 = netAdj(XInd(2,1):XInd(2,2),XInd(2,1):XInd(2,2));
K = params.K;

%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nloop = params.nloop;
[rX1,cX1] = size(X1);
[rX2,cX2] = size(X2);

bestW = zeros(rX1,K);
bestH1 = zeros(K,cX1);
bestH2 = zeros(K,cX2);
Comodule = cell(K,1+nbX);
records = cell(nloop,1);
iterNumList = zeros(nloop,1);

% Records some variables' values during running
fid = fopen([ResultsFile '/SNMNMF_RunRecords.txt'],'w+');
if(fid == -1)
    error('Fail in opening a text file.')
end

fprintf(fid,'%s\n',['K = ' int2str(K) ', thrNet11 = ' num2str(params.thrNet11) ...
    ', thrNet12 = ' num2str(params.thrNet12) ', thrNet22 = ' num2str(params.thrNet22) ...
    ', thrXr = ' num2str(params.thrXr) ', thrXc = ' num2str(params.thrXc) ...
    ', nloop = ' num2str(params.nloop) ', maxiter = ' num2str(params.maxiter) '.']);

bestObj = Inf(1,8);
file_path1='W/'
file_path2='H2/'
for iloop = 1:nloop
    fprintf(1,' Iteration %d\n',iloop);
    % Algorithm.
    [W,H1,H2,Obj,IterStepNum] = SNMNMF_algorithm(X1, X2, A11, A12, A22, params);
    str=[file_path1,num2str(iloop),'.txt']
    save(str,'W','-ascii')
    str1=[file_path2,num2str(iloop),'.txt']
   save(str1,'H2','-ascii')
    % Records the objective values in each step.
    records{iloop,1} = Obj;
    iterNumList(iloop) = IterStepNum;
    
    % If the objective values have improved, then update the variables W, H.
    % and records the new objective values.
    if(Obj(IterStepNum,8) < bestObj(8))
        bestObj = Obj(IterStepNum,:);
        bestW = W;
        bestH1 = H1;
        bestH2 = H2;
        fprintf(fid,'%s\n',['iloop = ' int2str(iloop) ...
            ', bestObj = [' num2str(bestObj(1:7)) ...
            '], sum_Obj = ' num2str(bestObj(8)) '.']); 
    end
end
fclose(fid);

% Compute the modules according to bestW, bestH1 and bestH2
W = bestW;
H1 = bestH1;
H2 = bestH2;
clear bestW bestH1 bestH2;
params.records = records;
params.iterNumList =  iterNumList;

% Output rule
thrd_module = params.thrd_module;
isdouble = params.isdouble;
% Get comodules(index): based on zscore
Comodule(:,1) = SNMNMF_module(W',thrd_module(1),0);

module1 = SNMNMF_module(H1,thrd_module(2),isdouble);
module2 = SNMNMF_module(H2,thrd_module(3),isdouble);
Comodule(:,2:3) = [module1,module2];
clear module1 module2

end

