% Global Analysis
clear all
clc

cd 'C:\Users\User\OneDrive - Universidade Estadual de Campinas\Doutorado_2024\Stenose-Chapter\Processing_Data'
load Pearson_CrossCorrelation_Lag_Matrices.mat
load 'StenoseClassification.mat'

%% Complete Cross correlation
% % Not necesary
% clear MtxCrossCorr_hb2
% for subj = 1:length(MtxPearson_hb)
%     for run = 1:length(MtxPearson_hb{subj})
%         for hb = 1:3
%             for line = 1:size(MtxPearson_hb{subj}{run},1)
%                 for column = 1:size(MtxPearson_hb{subj}{run},2)
% %                     if MtxPearson_hb{subj}{run}(line,column,hb) >= 0.5 
%                     if MtxCrossCorr_hb{subj}{run}(line,column,hb) == 0.5
%                         MtxCrossCorr_hb2{subj}{run}(line,column,hb) = MtxPearson_hb{subj}{run}(line,column,hb);
%                     else
%                         MtxCrossCorr_hb2{subj}{run}(line,column,hb) = MtxCrossCorr_hb{subj}{run}(line,column,hb);
%                     end
%                 end
%             end
%         end
%     end
% end
% clear subj run line column MtxCrossCorr_hb hb
% 
% MtxCrossCorr_hb = MtxCrossCorr_hb2; clear MtxCrossCorr_hb2

%% PLOT: compare Pearson and cross-correlation
% subj = 68; run = 1;
% hb = 3;
% subplot(1,4,1); imagesc(MtxPearson_hb{subj}{run}(:,:,hb),[-1 1]); colormap('jet');
% subplot(1,4,2); imagesc(MtxCrossCorr_hb{subj}{run}(:,:,hb),[-1 1]); colormap('jet');
% % subplot(1,4,3); imagesc(MtxCrossCorr_hb2{subj}{run}(:,:,hb),[-1 1]); colormap('jet');
% % subplot(1,4,4); imagesc(MtxCrossCorr_hb3{subj}{run}(:,:,hb),[-1 1]); colormap('jet');
% 
% % histogram(MtxCrossCorr_hb2{subj}{run}(:,:,hb))
% % hold on
% % histogram(MtxCrossCorr_hb3{subj}{run}(:,:,hb))


%% For graph analysis
% Diagonal by Zero on Pearson, Cross-correlation, and Lag
clc

% MtxCrossCorr_hbNaN = MtxCrossCorr_hb;
for subj = 1:length(MtxCrossCorr_hb)
    for run = 1:length(MtxCrossCorr_hb{subj})
        for hb = 1:3
            [n, m] = size(MtxCrossCorr_hb{subj}{run}(:,:,hb));
            for i = 1:min(n, m)
                MtxPearson_hb{subj}{run}(i,i,hb) = 1;%NaN;%0;
                MtxCrossCorr_hb{subj}{run}(i,i,hb) = 1;%0;
%                 MtxCrossCorr_hbNaN{subj}{run}(i,i,hb) = 1; %NaN -> para Median_CC_Subj
                MtxLag_hb{subj}{run}(i,i,hb) = 0;%NaN; %0; 
            end
            clear i n m
        end
    end
end
clear subj run hb

%% Apply threshold on Pearson, cross-correlation, and lag mtxs
clc

thr = 0.4;
clear MtxPearson_hb_thr MtxCrossCorr_hb_thr MtxLag_hb_thr
for subj = 1:length(MtxCrossCorr_hb)
    for run = 1:length(MtxCrossCorr_hb{subj})
        for hb = 1:3
            clear pearson_mtx crosscorr_mtx conn_mtx lag_mtx
            pearson_mtx = MtxPearson_hb{subj}{run}(:,:,hb);
            crosscorr_mtx = MtxCrossCorr_hb{subj}{run}(:,:,hb);
            conn_mtx = MtxCrossCorr_hb{subj}{run}(:,:,hb);
            % p >= 0.5 -> tt = 0 (já garantí en la 1ra Matriz de lag) 
            % p < 0.5, cc >= 0.5 -> tt > 0
            % p < 0.5, cc < 0.5 -> tt = NaN
            lag_mtx = MtxLag_hb{subj}{run}(:,:,hb);
            
            % Creation Pearson mtx with significant correlations
            for line = 1:size(pearson_mtx,1)
                for column = 1:size(pearson_mtx,2)
                    if pearson_mtx(line,column) >= thr 
                        axPearson(line,column) = 1; %pearson_mtx(line,column);
                    else
                        axPearson(line,column) = 0;%NaN;% 
                    end
                end
            end
            MtxPearson_hb_thr{subj}{run}(:,:,hb) = axPearson;
            clear line column axPearson
            
            % Creation CrossCorrelation mtx with significant correlations
            for line = 1:size(crosscorr_mtx,1)
                for column = 1:size(crosscorr_mtx,2)
                    if crosscorr_mtx(line,column) >= thr 
                        axCrossCorr(line,column) = 1; %crosscorr_mtx(line,column);
                    else
                        axCrossCorr(line,column) = 0;%NaN; 
                    end
                end
            end
            MtxCrossCorr_hb_thr{subj}{run}(:,:,hb) = axCrossCorr;
            clear line column axCrossCorr
            
            % Creation of Lag mtx with significant correlations
            for line = 1:size(conn_mtx,1)
                for column = 1:size(conn_mtx,2)
                    if conn_mtx(line,column) >= thr 
                        axLagNaN(line,column) = lag_mtx(line,column);
                    else
                        axLagNaN(line,column) = -1; %NaN;
                    end
                end
            end
            MtxLag_hb_thr{subj}{run}(:,:,hb) = axLagNaN;
            clear line column axLagNaN     
        end
    end
end
clear thr subj run hb conn_mtx MtxCrossCorr_hb
clear pearson_mtx crosscorr_mtx conn_mtx lag_mtx

% MtxPearson_hb_thr = MtxPearson_hb;
% MtxCrossCorr_hb_thr = MtxCrossCorr_hbNaN;

%% Invert Hemisphere: 
% Right hemisphere has most occlusion
cd 'G:\Meu Drive\Doctorado\PhD_Project\Fisiologia do Cérebro\Analise-Data-NIRS\Raw Data\Resting State Stenosis Project\New_Cohort'
InvertGroups = [group_21; group_31; group_32; group_41; group_42; group_43]'; 
for subj = 1:length(InvertGroups)
    for run = 1:length(InvertGroups(subj))
        for hb = 1:3
            % Pearson: lag zero 
            axPearson = MtxPearson_hb_thr{InvertGroups(subj)}{run}(:,:,hb);
            MtxPearson_hb_thr{InvertGroups(subj)}{run}(:,:,hb) = InvertHemisphere48Channels(axPearson);
            
            % Cross correlation: lag diferent to zero
            axCrossCorr = MtxCrossCorr_hb_thr{InvertGroups(subj)}{run}(:,:,hb);
            MtxCrossCorr_hb_thr{InvertGroups(subj)}{run}(:,:,hb) = InvertHemisphere48Channels(axCrossCorr);
            
            % Lag matriz
            axMtxNaN = MtxLag_hb_thr{InvertGroups(subj)}{run}(:,:,hb);
            MtxLag_hb_thr{InvertGroups(subj)}{run}(:,:,hb) = InvertHemisphere48Channels(axMtxNaN);
                        
            clear axMtxNaN axPearson axCrossCorr
        end
    end
end
clear subj run hb axMtxNaN InvertGroups

clear group_12 group_21 group_13 group_31 group_23 group_32 
clear group_33 group_43 group_34 group_24 group_41 group_42

%% Compute ROI Global effect
clc

clear medianPearsonCorrROI medianCrossCorrROI medianLagROI 
for subj = 1:length(MtxLag_hb_thr)
    for run = 1:length(MtxLag_hb_thr{subj})
        for hb = 1:3
            % % Pearson correlation
            %conn_mtx_PC = MtxPearson_hb_thr{subj}{run}(:,:,hb);
            %z_mtx_PC = 0.5*( log(1+conn_mtx_PC) - log(1-conn_mtx_PC) );
            %medianPearsonCorrROI{subj}{run}(:,:,hb) = CreateROIMtx(z_mtx_PC);
            %clear conn_mtx_PC z_mtx_PC
            medianPearsonCorrROI{subj}{run}(:,:,hb) = CreateROIMtx(MtxPearson_hb_thr{subj}{run}(:,:,hb));
            
            % % Cross correlation
            %conn_mtx_CC = MtxCrossCorr_hb_thr{subj}{run}(:,:,hb);
            %z_mtx_CC = 0.5*( log(1+conn_mtx_CC) - log(1-conn_mtx_CC) );
            %medianCrossCorrROI{subj}{run}(:,:,hb) = CreateROIMtx(z_mtx_CC);
            %clear conn_mtx_CC z_mtx_CC
            medianCrossCorrROI{subj}{run}(:,:,hb) = CreateROIMtx(MtxCrossCorr_hb_thr{subj}{run}(:,:,hb));
            
            % Lag matrix
            medianLagROI{subj}{run}(:,:,hb) = CreateROIMtx(MtxLag_hb_thr{subj}{run}(:,:,hb));
        end
    end
end
clear subj run hb 

%% Median subject and group
clc

% Compute median subject across all runs
clear median_Subj_Pearson median_Subj_CrossCorr median_Subj_Delay 
for subj = 1:length(medianPearsonCorrROI)
    ax_medianP = NaN(12,12,4); ax_medianCC = NaN(12,12,4); ax_medianD = NaN(12,12,4);
    for run = 1:length(medianPearsonCorrROI{subj})
        % % Pearson
        ax_medianP(:,:,run) = medianPearsonCorrROI{subj}{run}(:,:,3);
        
        % % Cross correlation
        ax_medianCC(:,:,run) = medianCrossCorrROI{subj}{run}(:,:,3);
        
        % % Lag mtx
        ax_medianD(:,:,run) = medianLagROI{subj}{run}(:,:,3);
    end
    median_Subj_Pearson{subj} = nanmedian(ax_medianP,3); % *Best results for PLOTS
    median_Subj_CrossCorr{subj} = nanmedian(ax_medianCC,3); % *Best results for PLOTS
    median_SubjD{subj} = nanmedian(ax_medianD,3); % *Best results for PLOTS
    %medianLag_run{subj} = nanmean(ax_medianLag,3);
end
clear subj run ax_medianP ax_medianCC ax_medianD

% Compute median group across all subject of specifics groups
control_P = NaN(12,12,length(gControl)); control_CC = NaN(12,12,length(gControl)); control_D = NaN(12,12,length(gControl));
g12_P = NaN(12,12,length(g12)); g12_CC = NaN(12,12,length(g12)); g12_D = NaN(12,12,length(g12));
g13_P = NaN(12,12,length(g13)); g13_CC = NaN(12,12,length(g13)); g13_D = NaN(12,12,length(g13));
g23_P = NaN(12,12,length(g23)); g23_CC = NaN(12,12,length(g23)); g23_D = NaN(12,12,length(g23));
g33_P = NaN(12,12,length(g33)); g33_CC = NaN(12,12,length(g33)); g33_D = NaN(12,12,length(g33));
g34_P = NaN(12,12,length(g34)); g34_CC = NaN(12,12,length(g34)); g34_D = NaN(12,12,length(g34));

for subj = 1:length(gControl)
    control_P(:,:,subj) = median_Subj_Pearson{gControl(subj)};
    control_CC(:,:,subj) = median_Subj_CrossCorr{gControl(subj)};
    control_D(:,:,subj) = median_SubjD{gControl(subj)};
end
for subj = 1:length(g12)
    g12_P(:,:,subj) = median_Subj_Pearson{g12(subj)};
    g12_CC(:,:,subj) = median_Subj_CrossCorr{g12(subj)};
    g12_D(:,:,subj) = median_SubjD{g12(subj)};
end
for subj = 1:length(g13)
    g13_P(:,:,subj) = median_Subj_Pearson{g13(subj)};
    g13_CC(:,:,subj) = median_Subj_CrossCorr{g13(subj)};
    g13_D(:,:,subj) = median_SubjD{g13(subj)};
end
for subj = 1:length(g23)
    g23_P(:,:,subj) = median_Subj_Pearson{g23(subj)};
    g23_CC(:,:,subj) = median_Subj_CrossCorr{g23(subj)};
    g23_D(:,:,subj) = median_SubjD{g23(subj)};
end
for subj = 1:length(g33)
    g33_P(:,:,subj) = median_Subj_Pearson{g33(subj)};
    g33_CC(:,:,subj) = median_Subj_CrossCorr{g33(subj)};
    g33_D(:,:,subj) = median_SubjD{g33(subj)};
end
for subj = 1:length(g34)
    g34_P(:,:,subj) = median_Subj_Pearson{g34(subj)};
    g34_CC(:,:,subj) = median_Subj_CrossCorr{g34(subj)};
    g34_D(:,:,subj) = median_SubjD{g34(subj)};
end

% % For Pearson and cross-correlation
% Median_Group{1} = tanh(nanmedian(ax_lagControl,3));
% Median_Group{2} = tanh(nanmedian(ax_lag_g12,3));
% Median_Group{3} = tanh(nanmedian(ax_lag_g13,3));
% Median_Group{4} = tanh(nanmedian(ax_lag_g23,3));
% Median_Group{5} = tanh(nanmedian(ax_lag_g33,3));
% Median_Group{6} = tanh(nanmedian(ax_lag_g34,3));

% % For average mtx
avg_P{1} = nanmedian(control_P,3); avg_CC{1} = nanmedian(control_CC,3); avg_D{1} = nanmedian(control_D,3);
avg_P{2} = nanmedian(g12_P,3); avg_CC{2} = nanmedian(g12_CC,3); avg_D{2} = nanmedian(g12_D,3);
avg_P{3} = nanmedian(g13_P,3); avg_CC{3} = nanmedian(g13_CC,3); avg_D{3} = nanmedian(g13_D,3);
avg_P{4} = nanmedian(g23_P,3); avg_CC{4} = nanmedian(g23_CC,3); avg_D{4} = nanmedian(g23_D,3);
avg_P{5} = nanmedian(g33_P,3); avg_CC{5} = nanmedian(g33_CC,3); avg_D{5} = nanmedian(g33_D,3);
avg_P{6} = nanmedian(g34_P,3); avg_CC{6} = nanmedian(g34_CC,3); avg_D{6} = nanmedian(g34_D,3);

clear control_P g12_P g13_P g23_P g33_P g34_P
clear control_CC g12_CC g13_CC g23_CC g33_CC g34_CC
clear control_D g12_D g13_D g23_D g33_D g34_D

%% PLOT
for k = 1:18
    ax(k) = subplot(3,6,k);
end

for subj = 1:6
    subplot(ax(subj));
    imagesc(avg_P{(subj)},[-1 1]); colormap('jet'); 
    
    subplot(ax(subj+6));
    imagesc(avg_CC{(subj)},[-1 1]); colormap('jet'); 
    
    subplot(ax(subj+12));
    imagesc(avg_D{(subj)},[-1 1]); colormap('jet'); 
end
clear k ax subj 

%% aply threshold
clc
thr = 0.4; clear adj_mtx_PC adj_mtx_CC
for subj = 1:70
    % for Pearson
    for line = 1:size(medianPearsonCorrROI{subj}{1}(:,:,3),1)
        for column = 1:size(medianPearsonCorrROI{subj}{1}(:,:,3),2)
            if medianPearsonCorrROI{subj}{1}(line,column,3) > thr
                ax_adj_mtx_PC(line,column) = 1;
            else
                ax_adj_mtx_PC(line,column) = 0;
            end
        end
    end
    clear line column 
    adj_mtx_PC{subj} = ax_adj_mtx_PC; clear ax_adj_mtx_PC
    
    % for CrossCorrelation
    for line = 1:size(medianCrossCorrROI{subj}{1}(:,:,3),1)
        for column = 1:size(medianCrossCorrROI{subj}{1}(:,:,3),2)
            if medianCrossCorrROI{subj}{1}(line,column,3) > thr
                ax_adj_mtx_CC(line,column) = 1;
            else
                ax_adj_mtx_CC(line,column) = 0;
            end
        end
    end
    clear line column 
    adj_mtx_CC{subj} = ax_adj_mtx_CC; clear ax_adj_mtx_CC    
end
clear subj

%% PLOT g34
for k = 1:8
    ax(k) = subplot(2,4,k);
end

for subj = 1:4
    subplot(ax(subj));
    imagesc(adj_mtx_PC{g34(subj)},[-1 1]); colormap('jet'); 
    
    subplot(ax(subj+4));
    imagesc(adj_mtx_CC{g34(subj)},[-1 1]); colormap('jet');    
end
clear k ax subj 

%% PLOT control
for k = 1:40
    ax(k) = subplot(2,20,k);
end

for subj = 1:20
    subplot(ax(subj));
    imagesc(adj_mtx_PC{gControl(subj)},[-1 1]); colormap('jet'); 
    
    subplot(ax(subj+20));
    imagesc(adj_mtx_CC{gControl(subj)},[-1 1]); colormap('jet');    
end
clear k ax subj 

%% Compute ROIs Local effect
% revisar zspace
clc

ROI.FS.L = [1 2 3]; ROI.FS.R = [25 26 27];
ROI.FM.L = [4 5 7 15 16]; ROI.FM.R = [28 29 31 39 40];
ROI.FI.L = [14 17 18 19 20 21]; ROI.FI.R = [38 41 42 43 44 45];
ROI.Pre.L = [6 8 9]; ROI.Pre.R = [30 32 33];
ROI.PST.L = [13 22 23 24]; ROI.PST.R = [37 46 47 48];
ROI.OA.L = [10 11 12]; ROI.OA.R = [34 35 36];

hb = 3;
clear ROIMtxNetwork
for subj = 1:length(MtxLag_hb_thr)
    for run = 1:length(MtxLag_hb_thr{subj})
%         ax_mtx = MtxPearson_hb_thr{subj}{run}(:,:,hb);
%         ax_mtx = MtxCrossCorr_hb_thr{subj}{run}(:,:,hb);
        ax_mtx = MtxLag_hb_thr{subj}{run}(:,:,hb);
        
        % % Left
        ROIMtxNetwork{subj}{run}{1} = ax_mtx(ROI.FS.L,ROI.FS.L);
        ROIMtxNetwork{subj}{run}{2} = ax_mtx(ROI.FM.L,ROI.FM.L);
        ROIMtxNetwork{subj}{run}{3} = ax_mtx(ROI.FI.L,ROI.FI.L);
        ROIMtxNetwork{subj}{run}{4} = ax_mtx(ROI.Pre.L,ROI.Pre.L);
        ROIMtxNetwork{subj}{run}{5} = ax_mtx(ROI.PST.L,ROI.PST.L);
        ROIMtxNetwork{subj}{run}{6} = ax_mtx(ROI.OA.L,ROI.OA.L);

        % % Right
        ROIMtxNetwork{subj}{run}{7} = ax_mtx(ROI.FS.R,ROI.FS.R);
        ROIMtxNetwork{subj}{run}{8} = ax_mtx(ROI.FM.R,ROI.FM.R);
        ROIMtxNetwork{subj}{run}{9} = ax_mtx(ROI.FI.R,ROI.FI.R);
        ROIMtxNetwork{subj}{run}{10} = ax_mtx(ROI.Pre.R,ROI.Pre.R);
        ROIMtxNetwork{subj}{run}{11} = ax_mtx(ROI.PST.R,ROI.PST.R);
        ROIMtxNetwork{subj}{run}{12} = ax_mtx(ROI.OA.R,ROI.OA.R);
        clear ax_mtx
    end
end
clear subj run hb

% cd 'C:\Users\User\OneDrive - Universidade Estadual de Campinas\Doutorado\Stenose-Chapter\Processing_Data'
% save("Matrix_with_Global_Effect.mat","medianPearsonCorrROI","medianCrossCorrROI","medianLagROI") 
% save("Matrix_with_Local_Effect.mat","ROIMtxNetwork_P","ROIMtxNetwork_CC","ROIMtxNetwork_Lag") 
