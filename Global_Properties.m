clear all
clc

cd 'C:\Users\User\OneDrive - Universidade Estadual de Campinas\Doutorado_2024\Stenose-Chapter\Processing_Data'
load Matrix_with_Global_Effect.mat
load StenoseClassification.mat

g12 = [group_12;group_21];
g13 = [group_13;group_31];
g23 = [group_23;group_32];
g33 = group_33;
g34 = [group_43;group_34];
gControl = [51:70];

clear group_12 group_21 group_13 group_31 group_23 group_32 
clear group_33 group_43 group_34 group_24 group_41 group_42

%% Median subject and group
clc

% Compute median subject across all runs
clear median_Subj
for subj = 1:length(medianPearsonCorrROI)
    ax_median = NaN(12,12,4);
    for run = 1:length(medianPearsonCorrROI{subj})
        % % Pearson
%         ax_median(:,:,run) = medianPearsonCorrROI{subj}{run}(:,:,3);
        
        % % Cross correlation
%         ax_median(:,:,run) = medianCrossCorrROI{subj}{run}(:,:,3);
        
        % % Lag mtx
        ax_median(:,:,run) = medianLagROI{subj}{run}(:,:,3);
    end
    median_Subj{subj} = nanmedian(ax_median,3); % *Best results for PLOTS
    %medianLag_run{subj} = nanmean(ax_medianLag,3);
end
clear subj run ax_median

% Compute median group across all subject of specifics groups
ax_lagControl = NaN(12,12,length(gControl));
ax_lag_g12 = NaN(12,12,length(g12));
ax_lag_g13 = NaN(12,12,length(g13));
ax_lag_g23 = NaN(12,12,length(g23));
ax_lag_g33 = NaN(12,12,length(g33));
ax_lag_g34 = NaN(12,12,length(g34));

for subj = 1:length(gControl)
    ax_lagControl(:,:,subj) = median_Subj{gControl(subj)};
end
for subj = 1:length(g12)
    ax_lag_g12(:,:,subj) = median_Subj{g12(subj)};
end
for subj = 1:length(g13)
    ax_lag_g13(:,:,subj) = median_Subj{g13(subj)};
end
for subj = 1:length(g23)
    ax_lag_g23(:,:,subj) = median_Subj{g23(subj)};
end
for subj = 1:length(g33)
    ax_lag_g33(:,:,subj) = median_Subj{g33(subj)};
end
for subj = 1:length(g34)
    ax_lag_g34(:,:,subj) = median_Subj{g34(subj)};
end

% % For Pearson and cross-correlation
% Median_Group{1} = tanh(nanmedian(ax_lagControl,3));
% Median_Group{2} = tanh(nanmedian(ax_lag_g12,3));
% Median_Group{3} = tanh(nanmedian(ax_lag_g13,3));
% Median_Group{4} = tanh(nanmedian(ax_lag_g23,3));
% Median_Group{5} = tanh(nanmedian(ax_lag_g33,3));
% Median_Group{6} = tanh(nanmedian(ax_lag_g34,3));

% % For lag
Median_Group{1} = nanmedian(ax_lagControl,3);
Median_Group{2} = nanmedian(ax_lag_g12,3);
Median_Group{3} = nanmedian(ax_lag_g13,3);
Median_Group{4} = nanmedian(ax_lag_g23,3);
Median_Group{5} = nanmedian(ax_lag_g33,3);
Median_Group{6} = nanmedian(ax_lag_g34,3);

clear ax_lagControl ax_lag_g12 ax_lag_g13 ax_lag_g23 ax_lag_g33 ax_lag_g34 subj
% clear median_Subj

%% Create adj mtx only for Pearson and Cross correlation

thr = 0; clear axMtx_Pearson axMtx_CrossCorr
for subj = 1:6
    % For Pearson
    for line = 1:size(Median_Group_Pearson{subj},1)
        for column = 1:size(Median_Group_Pearson{subj},2)
            if Median_Group_Pearson{subj}(line,column) > thr
                axMtx_Pearson{subj}(line,column) = 1;
            else
                axMtx_Pearson{subj}(line,column) = 0;
            end
        end
    end
    clear line column

    % For CrossCorr
    for line = 1:size(Median_Group_CrossCorr{subj},1)
        for column = 1:size(Median_Group_CrossCorr{subj},2)
            if Median_Group_CrossCorr{subj}(line,column) > thr
                axMtx_CrossCorr{subj}(line,column) = 1;
            else
                axMtx_CrossCorr{subj}(line,column) = 0;
            end
        end
    end
    clear line column
end
clear thr subj line column

%% PLOT
for k = 1:18
    ax(k) = subplot(3,6,k);
end

for subj = 1:6
    subplot(ax(subj));
    imagesc(axMtx_Pearson{(subj)},[0 1]); colormap('jet'); 
    
    subplot(ax(subj+6));
    imagesc(axMtx_CrossCorr{(subj)},[0 1]); colormap('jet'); 
    
    subplot(ax(subj+12));
    imagesc(Median_Group_Lag{(subj)},[0 1]); colormap('jet'); 
end
clear k ax subj 
% set(gcf,'Position',[10 800 1900 1800])

%%
thr = .5; clear adjMtx_PC adjMtx_CC
for subj = 1:70
    % % for Pearson correlation
    for line = 1:size(median_Subj_PC{subj},1)
        for column = 1:size(median_Subj_PC{subj},2)
            if median_Subj_PC{subj}(line,column) >= thr 
                adjMtx_PC{subj}(line,column) = 1;
            else
                adjMtx_PC{subj}(line,column) = 0;
            end
        end
    end
    clear line column
    
    % % for Cross correlation
    for line = 1:size(median_Subj_CC{subj},1)
        for column = 1:size(median_Subj_CC{subj},2)
            if median_Subj_CC{subj}(line,column) >= thr 
                adjMtx_CC{subj}(line,column) = 1;
            else
                adjMtx_CC{subj}(line,column) = 0;
            end
        end
    end
    clear line column
end
clear thr subj

% % For Control
% for k = 1:40
%     ax(k) = subplot(2,20,k);
% end
% 
% for subj = 1:length(gControl)
%     subplot(ax(subj));
%     imagesc(adjMtx_PC{gControl(subj)},[-1 1]); colormap('jet'); 
%     
%     subplot(ax(subj+20));
%     imagesc(adjMtx_CC{gControl(subj)},[-1 1]); colormap('jet'); 
% end
% clear k subj
% set(gcf,'Position',[10 800 1900 180])

% % For g12
% for k = 1:20
%     ax(k) = subplot(2,10,k);
% end
% 
% for subj = 1:length(g12)
%     subplot(ax(subj));
%     imagesc(adjMtx_PC{g12(subj)},[-1 1]); colormap('jet'); 
%     
%     subplot(ax(subj+10));
%     imagesc(adjMtx_CC{g12(subj)},[-1 1]); colormap('jet'); 
% end
% clear k subj
% set(gcf,'Position',[10 800 1900 180])

% % For g13
% for k = 1:34
%     ax(k) = subplot(2,17,k);
% end
% 
% for subj = 1:length(g13)
%     subplot(ax(subj));
%     imagesc(adjMtx_PC{g13(subj)},[-1 1]); colormap('jet'); 
%     
%     subplot(ax(subj+17));
%     imagesc(adjMtx_CC{g13(subj)},[-1 1]); colormap('jet'); 
% end
% clear k subj
% set(gcf,'Position',[10 800 1900 180])

% % For g23
% for k = 1:16
%     ax(k) = subplot(2,8,k);
% end
% 
% for subj = 1:length(g23)
%     subplot(ax(subj));
%     imagesc(adjMtx_PC{g23(subj)},[-1 1]); colormap('jet'); 
%     
%     subplot(ax(subj+8));
%     imagesc(adjMtx_CC{g23(subj)},[-1 1]); colormap('jet'); 
% end
% clear k subj
% set(gcf,'Position',[10 800 1900 180])

% % For g33
% for k = 1:10
%     ax(k) = subplot(2,5,k);
% end
% 
% for subj = 1:length(g33)
%     subplot(ax(subj));
%     imagesc(adjMtx_PC{g33(subj)},[-1 1]); colormap('jet'); 
%     
%     subplot(ax(subj+5));
%     imagesc(adjMtx_CC{g33(subj)},[-1 1]); colormap('jet'); 
% end
% clear k subj
% set(gcf,'Position',[10 800 1900 180])

% For g34
for k = 1:8
    ax(k) = subplot(2,4,k);
end

for subj = 1:length(g34)
    subplot(ax(subj));
    imagesc(adjMtx_PC{g34(subj)},[-1 1]); colormap('jet'); 
    
    subplot(ax(subj+4));
    imagesc(adjMtx_CC{g34(subj)},[-1 1]); colormap('jet'); 
end
clear k subj
set(gcf,'Position',[10 800 1900 180])

%% plot adjacency matrix

thr = .45; clear adjMtx
for group = 1:6
    for line = 1:size(Median_Group{group},1)
        for column = 1:size(Median_Group{group},2)
            if Median_Group{group}(line,column) >= thr 
                adjMtx{group}(line,column) = 1;
            else
                adjMtx{group}(line,column) = 0;
            end
        end
    end
end
clear thr group line column
        
subplot(1,6,1); imagesc(squeeze(adjMtx{1}),[-1 1]); colormap('jet'); colorbar
subplot(1,6,2); imagesc(squeeze(adjMtx{2}),[-1 1]); colormap('jet'); colorbar
subplot(1,6,3); imagesc(squeeze(adjMtx{3}),[-1 1]); colormap('jet'); colorbar
subplot(1,6,4); imagesc(squeeze(adjMtx{4}),[-1 1]); colormap('jet'); colorbar
subplot(1,6,5); imagesc(squeeze(adjMtx{5}),[-1 1]); colormap('jet'); colorbar
subplot(1,6,6); imagesc(squeeze(adjMtx{6}),[-1 1]); colormap('jet'); colorbar

set(gcf,'Position',[300 400 1500 150])

[sum(adjMtx{1}(:)) sum(adjMtx{2}(:)) sum(adjMtx{3}(:)) sum(adjMtx{4}(:)) sum(adjMtx{5}(:)) sum(adjMtx{6}(:))]/2


%% PLOT: for transit time network
clc

% % For lag
cte = 0.001;
ROIVector = [3 16 20 8 23 12 27 40 44 32 47 36];
g = graph(cte+Median_Group{1},'upper','omitselfloops');
ax_g12 = graph(cte+Median_Group{2},'upper','omitselfloops');
ax_g13 = graph(cte+Median_Group{3},'upper','omitselfloops');
ax_g23 = graph(cte+Median_Group{4},'upper','omitselfloops');
ax_g33 = graph(cte+Median_Group{5},'upper','omitselfloops');
ax_g34 = graph(cte+Median_Group{6},'upper','omitselfloops');

% % Para Cross Correlation
% axMtxOnes = ones(12,12)-.5;
% g = graph(cte+axMtxOnes,'upper','omitselfloops');
% ax_g12 = graph(cte+axMtxOnes,'upper','omitselfloops');
% ax_g13 = graph(cte+axMtxOnes,'upper','omitselfloops');
% ax_g23 = graph(cte+axMtxOnes,'upper','omitselfloops');
% ax_g33 = graph(cte+axMtxOnes,'upper','omitselfloops');
% ax_g34 = graph(cte+axMtxOnes,'upper','omitselfloops');

xmin = 0; xmax = 1;
% xmin = -1; xmax = 1;
axSizeNode = 0.01; axSizeLink = 3.5;%3.5; 
% factor = 5;
vecolor = [.5 .5 .5];

subplot(1,6,1)
p = plot(g,'-b','XData',ChaPos(ROIVector,1),'YData',ChaPos(ROIVector,2),'MarkerSize',axSizeNode)
p.EdgeCData = g.Edges.Weight;
p.LineWidth = axSizeLink;%axSizeLink*g.Edges.Weight; %  
p.NodeCData = diag(Median_Group{1}); 
caxis([xmin xmax]); labelnode(p,1:12,'')
axis('off'); hold on
plot(ChaPos(ROIVector,1),ChaPos(ROIVector,2),'ko','MarkerSize',15,'MarkerEdgeColor',vecolor,'Linewidth',1.0); 
clear p h g 

subplot(1,6,2)
p12 = plot(ax_g12,'-b','XData',ChaPos(ROIVector,1),'YData',ChaPos(ROIVector,2),'MarkerSize',axSizeNode)
p12.EdgeCData = ax_g12.Edges.Weight; 
p12.LineWidth = axSizeLink;% axSizeLink*ax_g12.Edges.Weight; % 
p12.NodeCData = diag(Median_Group{2}); 
caxis([xmin xmax]); labelnode(p12,1:12,'')
axis('off'); hold on
plot(ChaPos(ROIVector,1),ChaPos(ROIVector,2),'ko','MarkerSize',15,'MarkerEdgeColor',vecolor,'Linewidth',1.0); 
clear p12 ax_g12

subplot(1,6,3)
p13 = plot(ax_g13,'-b','XData',ChaPos(ROIVector,1),'YData',ChaPos(ROIVector,2),'MarkerSize',axSizeNode)
p13.EdgeCData = ax_g13.Edges.Weight; 
p13.LineWidth = axSizeLink;%axSizeLink*ax_g13.Edges.Weight; %  
p13.NodeCData = diag(Median_Group{3}); 
caxis([xmin xmax]); labelnode(p13,1:12,'')
axis('off'); hold on
plot(ChaPos(ROIVector,1),ChaPos(ROIVector,2),'ko','MarkerSize',15,'MarkerEdgeColor',vecolor,'Linewidth',1.0); 
clear p13 ax_g13

subplot(1,6,4)
p23 = plot(ax_g23,'-b','XData',ChaPos(ROIVector,1),'YData',ChaPos(ROIVector,2),'MarkerSize',axSizeNode)
p23.EdgeCData = ax_g23.Edges.Weight; 
p23.LineWidth = axSizeLink;% axSizeLink*ax_g23.Edges.Weight; %  
p23.NodeCData = diag(Median_Group{4}); 
caxis([xmin xmax]); labelnode(p23,1:12,'')
axis('off'); hold on
plot(ChaPos(ROIVector,1),ChaPos(ROIVector,2),'ko','MarkerSize',15,'MarkerEdgeColor',vecolor,'Linewidth',1.0); 
clear p23 ax_g23

subplot(1,6,5)
p33 = plot(ax_g33,'-b','XData',ChaPos(ROIVector,1),'YData',ChaPos(ROIVector,2),'MarkerSize',axSizeNode)
p33.EdgeCData = ax_g33.Edges.Weight; 
p33.LineWidth = axSizeLink;%  axSizeLink*ax_g33.Edges.Weight; %
p33.NodeCData = diag(Median_Group{5}); 
caxis([xmin xmax]); labelnode(p33,1:12,'')
axis('off'); hold on
plot(ChaPos(ROIVector,1),ChaPos(ROIVector,2),'o','MarkerSize',15,'MarkerEdgeColor',vecolor,'Linewidth',1.0); 
clear p33 ax_g33

subplot(1,6,6)
p34 = plot(ax_g34,'-b','XData',ChaPos(ROIVector,1),'YData',ChaPos(ROIVector,2),'MarkerSize',axSizeNode)
p34.EdgeCData = ax_g34.Edges.Weight; %colormap(flipud(autumn)) %colormap(copper) % 
p34.LineWidth = axSizeLink;% axSizeLink*ax_g34.Edges.Weight; % 
p34.NodeCData = diag(Median_Group{6}); colormap(flipud(copper)) %colormap(flipud(autumn)); %colormap(jet); %  
% colorbar
caxis([xmin xmax]); labelnode(p34,1:12,'')
axis('off'); hold on
plot(ChaPos(ROIVector,1),ChaPos(ROIVector,2),'o','MarkerSize',15,'MarkerEdgeColor',vecolor,'Linewidth',1.0); 
clear p34 ax_g34

set(gcf,'Position',[300 300 1500 200])

clear xmax xmin ROIVector axNodeSixe axSizeLink axSizeNode axMtxOnes vecolor

%% For graph considering median_Sub
% Compute for group mtx
% 3) Degree density (Medida de centralidade)  
% 4) Clustering coefficient (Medida de segregação)
% 5) Efficiency (Medida de integração)

% % Use Plot_Avg_Network.m to generate plots

clc
% cd 'G:\Meu Drive\Doctorado\PhD_Project\Tese'

threshold = [0.05:.05:0.95];
clear GlobalStrength KDensity AvgCC lambda eff
for subj = 1:length(median_Subj)
    for thr = 1:length(threshold)
        % Pearson
        clear adj_mtx
        conn_mtx = tanh(median_Subj{subj});
        for line = 1:size(conn_mtx,1)
            for column = 1:size(conn_mtx,2)
                if conn_mtx(line,column) > threshold(thr)
                    adj_mtx(line,column) = 1;
                else
                    adj_mtx(line,column) = 0;
                end
            end
        end
        clear line column 
                
        % Diagonal for 0
        [n, m] = size(adj_mtx);
        for ii = 1:min(n, m)
            adj_mtx(ii,ii) = 0;
        end
        clear m n ii
        
       
        % Average global node strength
        GlobalStrength(:,subj) = nanmedian(nansum(conn_mtx));
                    
        % Degree density
        K_density(subj,thr) = density_und(adj_mtx);  

        % Global clustering coefficient
        clear axCC
        axCC = clustering_coef_bu(adj_mtx);
        AvgCC(subj,thr) = median(axCC);
               
        % Global efficiency
        [lambda(subj,thr),eff(subj,thr)] = charpath(distance_bin(adj_mtx),0,1);
    end    
end
clear subj thr 
clear axCC_health axCC_occluded adj_lag adj_lag_health





