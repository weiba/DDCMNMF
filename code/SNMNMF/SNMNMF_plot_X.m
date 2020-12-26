function SNMNMF_plot_X(X,XInd,fig,figure_title,colormap_type)
%
% INPUT:
%
% X: (m x n) matrix.
% XInd: (2 x 2) matrix, the column index for each X block.
% fig: Figure No.
% figure_title: The title of figure.
% colormap_type: A parameter for heatmaps of these nbX matrices. 
%                Options for it are: 'blue-yellow'; 'green-red';
%                'yellow'; 'blue-white-red'; 'default'.
%

K = size(X,1); % #samples
nbX = size(XInd,1); % #blocks in X

sampleSort = [1:K]'; 

nTopRank = size(X,2);

listTopRankSample = 1:min(K,nTopRank);

figure(fig);
set(gcf,'name',figure_title);
C = exp_colormap(colormap_type,64);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% X %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for aa = 1:nbX,
    subplot(1,nbX,aa);
    coli = XInd(aa,1):XInd(aa,2);
    pX = X(:,coli);
    h = image(rescale(pX(sampleSort(listTopRankSample),:))*64);
    colormap(gca,C);

    set(gca,'FontSize',14);
    set(gca,'XTick',zeros(0,1),'XTickLabel',{},...
        'YTick',zeros(0,1),'YTickLabel',{})
    title(['X' num2str(aa)]);
end
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


