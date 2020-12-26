function SNMNMF_plot_results(X,XInd,FeatureType,fig,figure_title,colormap_type,vectorForRank)
%
% INPUT:
%
% X: (m x n) matrix.
% XInd: (2 x 2) matrix, the column index for each X block.
% FeatureType: (1 x 2) cell array, the type names for two types of features 
%              in matrix X.
% fig: Figure No.
% fig_title: The title of this figure.
% colormap_type: A parameter for heatmaps of these matrices.
%                Options for it are: 'blue-yellow'; 'green-red';
%                'yellow'; 'blue-white-red'; 'default'.
% vectorForRank: Structure variable, containing 
%                vectorForRank.w,vectorForRank.h,vectorForRank.comodule,
%                vectorForRank.hInd.
%

comodule = vectorForRank.comodule;
w = vectorForRank.w;
h = vectorForRank.h;
hInd = vectorForRank.hInd;
nbX = size(hInd,1);
rbX = size(X,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reorder the selected samples and randomly select the same number of
% samples.
selected_samples = comodule{1};
selected_samples = sort_vector(w(selected_samples),selected_samples);
if (isempty(selected_samples)),
    disp('There is no samples in the identified comodule.');
    return;
end

popu = setdiff(1:rbX,selected_samples);
nSample = length(selected_samples);
if(length(popu) >= nSample)
    rand_ri = randsample(popu,nSample);
else
    rand_ri = randsample(popu,nSample,true);
end
clear popu nSample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(fig);
set(gcf,'name',figure_title);
C = exp_colormap(colormap_type,64);

% For each X block, reorder the selected features and then select the same
% number of other features.
for bXi = 1:nbX
    coli = XInd(bXi,1):XInd(bXi,2);

    pX = X(:,coli);
    obj_hi = comodule{bXi + 1};
    sta_bh = hInd(bXi,1);
    end_bh = hInd(bXi,2);
    temp = end_bh - sta_bh + 1;
    cbX = length(coli);
    
    if(cbX == temp)
        obj_h = h(sta_bh:end_bh);
        obj_h_sign = ones(cbX,1);
    elseif(2*cbX == temp)
        pos_h = h(sta_bh:(sta_bh+cbX-1));
        neg_h = h((sta_bh+cbX):end_bh);
        [obj_h,obj_h_sign] = max([pos_h,neg_h],[],2);
    else
        error('');
    end
    
    obj_hi = sort_vector(obj_h(obj_hi),obj_hi);
    obj_hi_sign = obj_h_sign(obj_hi);
    selected_features = [obj_hi(obj_hi_sign == 1);obj_hi(obj_hi_sign == 2)];
    objMat = pX(selected_samples,selected_features);
    
    popu = setdiff(1:cbX,selected_features);
    nFeature = length(selected_features);
    if(length(popu) >= nFeature)
        rand_ci = randsample(popu,nFeature);
    else
        rand_ci = randsample(popu,nFeature,true);
    end
    clear popu nFeature
    
    tempMat21 = pX(rand_ri,selected_features);
    tempMat12 = pX(selected_samples,rand_ci);
    tempMat22 = pX(rand_ri,rand_ci);
    tempMat = [objMat tempMat12;tempMat21 tempMat22];
    clear tempMat21 tempMat12 tempMat22
    
    subplot(1,nbX,bXi);
    
    image(rescale(tempMat)*64);
    colormap(gca,C);
    set(gca,'FontSize',14,'XTick',zeros(1,0),'XTickLabel',{},...
        'YTick',zeros(1,0),'YTickLabel',{});
    % Using a rectangle to box the module.
    rectangle('Position',[0.6,0.6,length(selected_features)-0.1,...
        length(selected_samples)-0.1],'LineWidth',4,'EdgeColor',[1 1 0])
    
    title(['X' num2str(bXi)]);
    xlabel({[FeatureType{bXi} ' subset '],['(2 x ' num2str(length(obj_hi)) ')']});
    if(bXi == 1)
        ylabel(['Sample (2 x ' num2str(length(selected_samples)) ')']);
    end
    clear tempMat
end
end

function newLabels = sort_vector(vec,labels)
[tmp,sorti] = sort(vec,'descend');
newLabels = labels(sorti);
end

function newX = rescale(X)
minX = min(min(X));
maxX = max(max(X));
if (minX == maxX),
    newX = X;
else
    newX = (X-minX)/(maxX-minX);
end
end
