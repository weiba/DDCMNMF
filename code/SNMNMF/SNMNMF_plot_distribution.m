function SNMNMF_plot_distribution(nSample,XInd,Comodule,FeatureType,fig,figure_title)
%
% Module size distribution
%
% INPUT :
%
% nSample : #samples
% XInd: (2 x 2) matrix. the column index for each X block.
% Comodule : (K x 3) cell array, each row represents a identified
%             md-module.
% fig: Figure No.
% figure_title: The title of this figure.
% FeatureType: (1 x 2) cell array, the type names for two types of features 
%              in matrix X.
%
[K,nblock] = size(Comodule);
nbX = size(XInd,1);

figure(fig);
set(gcf,'name',figure_title);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The size distribution for samples in all the md-modules.
moduleMat = zeros(nSample,K);
for imodule = 1:K
    ele = Comodule{imodule,1};
    moduleMat(ele,imodule) = 1;
end
sta1 = sum(moduleMat,1);

subplot(1,nblock,1);
[nelements,centers] = hist(sta1,10);
bar(centers,nelements,'FaceColor',[1 0.84 0],...
    'EdgeColor',[0.5,0.5,0.5],'BarWidth',0.8)
title('Samples','FontSize',14)
ylabel('Frequency','FontSize',14)

maxVal = max(nelements);
ylim([0,maxVal+1])
yInterVal = max(1,fix(maxVal/2));

set(gca,'YTick',[0,yInterVal:yInterVal:maxVal],...
    'YTickLabel',[0,yInterVal:yInterVal:maxVal])

clear sta1 moduleMat maxVal yInterVal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The size distributions for each type of features in all the md-modules.
for aa = 1:nbX
    coli = XInd(aa,1):XInd(aa,2);
    rbX = length(coli);
    moduleMat = zeros(rbX,K);
    for imodule = 1:K
        ele = Comodule{imodule,aa+1};
        moduleMat(ele,imodule) = 1;
    end
    sta1 = sum(moduleMat,1);
    
    subplot(1,nblock,1+aa);
    [nelements,centers] = hist(sta1,10);
    bar(centers,nelements,'FaceColor',[1 0.84 0],...
        'EdgeColor',[0.5 0.5 0.5],'BarWidth',0.8)
    title(FeatureType{aa},'FontSize',14)
    
    maxVal = max(nelements);
    ylim([0,maxVal+1])
    yInterVal = max(1,fix(maxVal/2));
    
    set(gca,'YTick',[0,yInterVal:yInterVal:maxVal],...
        'YTickLabel',[0,yInterVal:yInterVal:maxVal])
    
    clear sta1 moduleMat
    
end
end