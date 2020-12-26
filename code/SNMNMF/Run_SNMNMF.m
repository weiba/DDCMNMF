function Run_SNMNMF(Input)
%
% This algorithm is only used for two input matrices, X1, X2 with the same
% rows, and if the interaction networks G1 referring to the interactions
% between the features of X1 and G2 for the interactions between the
% features of X2, G3 for the interactions between the features of X1 and X2
% are available, user could incorporate them into the integrative analysis.
%
% For example,
% >> clear all
% >> Run_SNMNMF(Input);
%
% INPUT:
%
% Input: A structure variable, including a number of components:
% Input.data: It should be a (m x n) matrix. Each row represents a sample,
%              and each column is a feature. Input.data = [X1,X2];
% Input.XBlockInd: It should be a (2 x 2) matrix. Each row represents
%                   the column index (start index, end index) of each X
%                   block in the Input.data.
%
% Input.netAdj : It is the adjacency matrix for the network which represents
%                the relationship between all the features in the input 
%                data. A11(n1 x n1),A22(n1 x n2),A12(n1 x n2) are the  
%                adjacency matrices for G1,G2,G3,respectively.
%                Thus, Input.netAdj = [A11 A12; A12' A22].If some network
%                is not available, replace it by zero matrix.
%
% Input.FeatureLabel: (n x 1) cell array, the labels of features in 
%                      Input.data.
% Input.SampleLabel: (m x 1) cell array, the labels of samples in 
%                     Input.data.
% Input.FeatureType: (1 x 2) cell array, the type names for features in
%                    Input.data.
% Input.params: A structure variable, including the parameters used in
%                this algorithm, including 
% Input.params.thrd_module: 
% three thresholds to select elements in two X blocks to construct 
% md-modules, and the first threshold is for selecting samples.
% Input.params.NClsuter: A pre-defined number of identifiey md-modules.
% Input.params.nloop: The repeating times of SNMNMF algorithm.
% Input.params.maxiter: The maximal number of iterations for this algorithm.
% Input.params.tol: The precision for convergence of algorithm.
% Input.params.thrNet11, Input.params.thrNet12, Input.params.thrNet22:
% three parameters for the network constraints.
% Input.params.thrXr, Input.params.thrXc: two parameters for the sparsity
% of factorized matrices W,H.
%
if nargin == 0
    help Run_SNMNMF
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Load data for input %%%%%%%%%%%%%%%%%%%%%%%%%%%
X = Input.data;
XInd = Input.XBlockInd; 
netAdj = Input.netAdj; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test for data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbX = size(XInd,1);
if(nbX ~= 2)
    error('Your input data X should contain two and only two blocks.')
end

[rNet,cNet] = size(netAdj);
[rX,cX] = size(X);
if(cX ~= cNet || rNet ~= cX)
    error('The dimension of adjacent network does not match with the input matrix!')
end
clear cX rNet cNet netAdj

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ResultsFile = 'SNMNMF_Results';
% To make a new folder to store all the results.
if ~isdir(ResultsFile)
    mkdir(ResultsFile);
end

% Plot the heatmap of original data.
colormap_type = 'blue-white-red';
SNMNMF_plot_X(X,XInd,100,'Original data',colormap_type);

saveas(gcf,[ResultsFile '/Original_data.fig'])
close(figure(100))

%%%%%%%%%%%%%%%%%% All the parameters used in this method %%%%%%%%%%%%%%%%%

params.thrNet11 = Input.params.thrNet11; % constraints for network A11
params.thrNet12 = Input.params.thrNet12; % constraints for network A12
params.thrNet22 = Input.params.thrNet22; % constraint for network A22
params.thrXr = Input.params.thrXr; % limit the growth of W
params.thrXc = Input.params.thrXc; % limit the growth of H
params.thrd_module = Input.params.thrd_module; % thresholds for selecting features
params.K = Input.params.NCluster; % the number of identified modules.

% Some parameters about iteration
params.nloop = Input.params.nloop; % 50
params.maxiter = Input.params.maxiter; % 500
params.tol = Input.params.tol; % 10^(-6)

%%%%%%%%%%%%%%%%%%%%%%%%% Test for parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if((nbX+1) ~= length(params.thrd_module)),
    error(['Error in variable ''thrd_module'': please input ' ...
        num2str(1+nbX) ' values as thresholds for selecting ' ...
        num2str(1+nbX) ' components.'])
end

XBlockSize = zeros(nbX,1);
for aa = 1:nbX,
    XBlockSize(aa) = XInd(aa,2) - XInd(aa,1) + 1;
end

if(params.K > min(min(XBlockSize),rX))
    error('Please input the NCluster such that NCluster < minimum of the number of features in all the X blocks')
end

%%%%%%%%%%%%%%%%%%%%%% Pre-process the real input data %%%%%%%%%%%%%%%%%%%%
% To ensure the non negativity of input matrices.
[newInput,isdouble] = SNMNMF_PrepData(Input);
params.isdouble = isdouble;
clear isdouble

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run SNMNMF algorithm, and output the comodules.
[W,H1,H2,Comodule,params] = SNMNMF_comodule(newInput,params);

save([ResultsFile '/' ResultsFile '.mat'],'Comodule','W','H1','H2','params')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the heatmap of a certain module
% User can set an arbitrary value in [1:1:params.K]
SelectIDForExample = 1; 

H = [H1,H2];
temp = Comodule(SelectIDForExample,:);
while((isempty(temp{1}) || isempty(temp{2}) || isempty(temp{3})) && (SelectIDForExample < params.K))
    SelectIDForExample = SelectIDForExample + 1;
    temp = Comodule(SelectIDForExample,:);
end
clear temp
if(SelectIDForExample >= params.K)
    error('Please set the parameter (thrd_module) smaller.')
else
    vectorForRank.comodule = Comodule(SelectIDForExample,:);
    vectorForRank.w = W(:,SelectIDForExample);
    vectorForRank.h = H(SelectIDForExample,:)';
    vectorForRank.hInd = newInput.XBlockInd;
    
    SNMNMF_plot_results(X,XInd,Input.FeatureType,200 + SelectIDForExample,...
        ['Identified comodule ' num2str(SelectIDForExample)],...
        colormap_type,vectorForRank)
    
    saveas(gcf,[ResultsFile '/Identified_comodule_' num2str(SelectIDForExample) '.fig'])
    close(figure(200+SelectIDForExample))
    
    clear vectorForRank
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the correlation between original data X and reconstructed data
% newX = W*H.

newX = [W*H1,W*H2];
newXInd = newInput.XBlockInd;

figure_title = 'Correlations between original space and reconstructed space';
corrMat = SNMNMF_plot_correlation(X,newX,XInd,newXInd,300,figure_title);

figure_title = 'Correlations_between_original_space_and_reconstructed_space';
saveas(gcf,[ResultsFile '/' figure_title '.fig'])
close(figure(300))

clear newX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output the co-modules (feature indexes) to a text file.
OutputModule2TXT(Comodule,Input.FeatureType,[ResultsFile '/' ResultsFile]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output the co-modules (feature labels) to an array of text files.

% For samples
SampleIndex = Comodule(:,1);
SampleLabel = Input.SampleLabel;
TypeName = 'Sample';
Index2LabelForModuleContent(SampleIndex,SampleLabel,TypeName,ResultsFile)

% For the first type of feature
GeneIndex = Comodule(:,2);
GeneLabel = Input.FeatureLabel(XInd(1,1):XInd(1,2));
TypeName = Input.FeatureType{1};
Index2LabelForModuleContent(GeneIndex,GeneLabel,TypeName,ResultsFile)

% For the second type of feature
miRNAIndex = Comodule(:,3);
miRNALabel = Input.FeatureLabel(XInd(2,1):XInd(2,2));
TypeName = Input.FeatureType{2};
Index2LabelForModuleContent(miRNAIndex,miRNALabel,TypeName,ResultsFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Module size distribution
nSamples = size(X,1);
figure_title = 'Module size distributions';
SNMNMF_plot_distribution(nSamples,XInd,Comodule,Input.FeatureType,400,figure_title);

figure_title = 'Module_size_distributions';
saveas(gcf,[ResultsFile '/' figure_title '.fig'])
close(figure(400))

end