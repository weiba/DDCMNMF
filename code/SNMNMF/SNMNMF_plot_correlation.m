function corrMat = SNMNMF_plot_correlation(X,newX,XInd,newXInd,fig,figure_title)
%
% This function is used for the box-plot of sample-wise correlations of
% original data (X) and reconstructed data (newX).
%
% INPUT:
%
% X: (m x n) matrix.
% newX: (m x n) matrix, newX = W*H.
% XInd: (2 x 2) matrix, the column index for each X block.
% newXInd: Similar to XInd, it is for newX.
% fig: Figure No.
% fig_title: The title of this figure.
%
% OUTPUT:
%
% corrMat: (m x 2) matrix, recording the correlations between original 
%          data (X) and reconstructed data (newX).
%

nbX = size(XInd,1);
rX = size(X,1);

corrMat = zeros(rX,nbX);

for bXi = 1:nbX
    coli = XInd(bXi,1):XInd(bXi,2);
    pX = X(:,coli);
    
    sta_newX = newXInd(bXi,1);
    end_newX = newXInd(bXi,2);
    temp = end_newX - sta_newX + 1;
    cbX = length(coli);
    
    if(cbX == temp)
        pnewX = newX(:,sta_newX:end_newX);
    elseif(2*cbX == temp)
        pos_newX = newX(:,sta_newX:(sta_newX+cbX-1));
        neg_newX = newX(:,(sta_newX+cbX):end_newX);
        pnewX = pos_newX - neg_newX;
    else
        error('');
    end
    
    for i = 1:size(X,1)
        corrMat(i,bXi) = corr(pX(i,:)',pnewX(i,:)');
    end
    
end

figure(fig);
set(gcf,'name',figure_title);

label = cell(1,nbX);

for i = 1:nbX
    label{i} = ['X' num2str(i)];
end

h = boxplot(corrMat,'colors',[0.5,0.5,0.5]);
set(h(7,:),'Visible','off')
h = findobj(gca,'Tag','Box');
for j = 1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),'b','FaceAlpha',.2,'FaceColor','blue');
   line(get(h(j),'XData'),get(h(j),'YData'),'Color',[0.5 0.5 0.5]);
end
clear h

set(gca,'FontSize',14);
set(gca,'XTick',1:nbX,'XTickLabel',label);
ylabel('Correlations')

end
