%% Plot Generator for FP control images
% Version July 2018

% Linear regression and b1 calculation are now only in SpectralLeakageCorrector_3.mat.
% Save data from this section as 09_SpectralLeakageCorrector.mat.

clear all;
% load relevant data
analysis_kinetics = uigetdir; % choose folder DMSPkinetics analysis
analysis_fpctrl = uigetdir; % choose folder FPctrl analysis

load(strcat(analysis_fpctrl,'\RegressionSlope_b1_final.mat')); % load b1 from FPctrl folder
load(strcat(analysis_kinetics,'\06_FluoCellComputation_analysisonly.mat')); % load cell recognition data from DMSPkinetics
saveplace = strcat(analysis_kinetics,'/leakage figures'); mkdir(saveplace);

% specify colors (RGB) and plot styles
% http://shirt-ediss.me/matlab-octave-more-colours/
load H:\PlotAppearance_Variables


%% Perform spectral leakage correction
% Perform spectral leakage correction BEFORE YFP threshold.
% Data to feed into this part of the code should be "fintcells" (straight
% from cell mask generation). 

% measured signal = (matrix of slopes) x (actual signal)
% Y = b1*X  --> solve for X
% X = Y*b1 or X = b1\Y

allc_fintcells = {};
for t = 1:total_time;
    for xy = 1:total_xy;
        c = length(tfp_fintcells{t,xy});
        for cell = 1:c;
            
            % 3x1 vector containing the fluo value of a cell in each color
            allc_fintcells{t,xy}{1,cell}(1,1) = tfp_fintcells{t,xy}(1,cell);
            allc_fintcells{t,xy}{1,cell}(2,1) = yfp_fintcells{t,xy}(1,cell);
            allc_fintcells{t,xy}{1,cell}(3,1) = rfp_fintcells{t,xy}(1,cell);
            
            % solve the matrix
            leakcorrect_fintcells{t,xy}{1,cell} = b1 \ allc_fintcells{t,xy}{1,cell}; 
            
            % separate
            tfp_leakcorr_fintcells{t,xy}(1,cell) = leakcorrect_fintcells{t,xy}{1,cell}(1,1);
            yfp_leakcorr_fintcells{t,xy}(1,cell) = leakcorrect_fintcells{t,xy}{1,cell}(2,1);
            rfp_leakcorr_fintcells{t,xy}(1,cell) = leakcorrect_fintcells{t,xy}{1,cell}(3,1);
             
            % difference between pre-correction and post-correction
            tfp_leakcorrdiff_fintcells{t,xy}{1,cell} = leakcorrect_fintcells{t,xy}{1,cell}(1,1) - allc_fintcells{t,xy}{1,cell}(1,1);
            yfp_leakcorrdiff_fintcells{t,xy}{1,cell} = leakcorrect_fintcells{t,xy}{1,cell}(2,1) - allc_fintcells{t,xy}{1,cell}(2,1);
            rfp_leakcorrdiff_fintcells{t,xy}{1,cell} = leakcorrect_fintcells{t,xy}{1,cell}(3,1) - allc_fintcells{t,xy}{1,cell}(3,1);
            
        end
    end
end

% collect amount of correction of each cell by expt condition
for q = 1:total_cond
    cond = expt_conditions_xy{q};
    tfp_temp = []; % initialize
    yfp_temp = [];
    rfp_temp = [];
    for xy = cond(1):cond(end);
        for t = 1:total_time;
    
        t_temp = cell2mat(tfp_leakcorrdiff_fintcells{t,xy});   
        tfp_temp = [tfp_temp,t_temp];

        y_temp = cell2mat(yfp_leakcorrdiff_fintcells{t,xy});   
        yfp_temp = [yfp_temp,y_temp];

        r_temp = cell2mat(rfp_leakcorrdiff_fintcells{t,xy});   
        rfp_temp = [rfp_temp,r_temp];
        
        end
    end

    tfp_leakcorrdiff_fintcells_q{q} = tfp_temp; 
    yfp_leakcorrdiff_fintcells_q{q} = yfp_temp;
    rfp_leakcorrdiff_fintcells_q{q} = rfp_temp;
    
    clearvars t_temp y_temp r_temp 
    clearvars tfp_temp yfp_temp rfp_temp 
end


%% visualize leakage correction magnitude for each color
% visualize the amount of correction with histogram, separated by color and
% for all experimental conditions

close all;
cmap = hsv(total_cond);
SaveFig = 1; % 1 = save figure

for FPcolors = 1:3;
    if FPcolors == 1;
        data = rfp_leakcorrdiff_fintcells_q;
        xaxislim = [-1 0.5];
        color = 'RFP';
    elseif FPcolors == 2;
        data = yfp_leakcorrdiff_fintcells_q;
        xaxislim = [-0.8 0.1];
        color = 'YFP';
    elseif FPcolors == 3;
        data = tfp_leakcorrdiff_fintcells_q;
        xaxislim = [-80 5]; % should be biggest correction
        color = 'TFP';
    end
    fig_hist(FPcolors) = figure('units','normalized','outerposition',[0 0 1 1]); 
    
    for q = 1:total_cond;

        histdata_q(q) = histogram(data{q},'FaceColor',cmap(q,:)); hold on;
            centers = histdata_q(q).BinEdges + histdata_q(q).BinWidth/2; % center of each bar in histogram
            heights = [histdata_q(q).Values,0];        
                hold on
            line_q(q) = plot(centers,heights,'LineWidth',3,'Color',cmap(q,:));    

    end
    title(strcat('magnitude of spectral leakage correction (',color,'); final value - original'),'FontSize',20)
    xlabel('magnitude of spectral leakage correction (final value minus original)','FontSize',20);
    ylabel('frequency','FontSize',20)
    set(histdata_q(:),'facealpha',0,'edgecolor','none'); % turn OFF bars
    legend(line_q(:),{expt_conditions_string{:}},'FontSize',16,'Location','NorthWest'); 
    xlim(xaxislim);
    % save figure
    if SaveFig == 1;
        saveas(fig_hist(FPcolors),strcat(saveplace,'\',color,'_LeakCorrectionDeltaHistogram.png'));
        close all;
    end
end

clearvars cmap FPcolors SaveFig data xaxislim color fig_hist q 
clearvars ans histdata_q centers heights line_q
 

%% calculate mean, std, and sem with spectral leakage corrected values

for xy = 1:total_xy;
    for t = 1:total_time
    
    % rfp
    r_temp = rfp_leakcorr_fintcells{t,xy};    
    rfp_leakcorr_fintview(t,xy) = mean(r_temp);
    rfp_leakcorr_sem(t,xy) = std(r_temp)/sqrt(length(r_temp)); 
    rfp_leakcorr_std(t,xy) = std(r_temp);
    % yfp
    y_temp = yfp_leakcorr_fintcells{t,xy};    
    yfp_leakcorr_fintview(t,xy) = mean(y_temp);
    yfp_leakcorr_sem(t,xy) = std(y_temp)/sqrt(length(y_temp)); 
    yfp_leakcorr_std(t,xy) = std(y_temp);    
    % tfp
    t_temp = tfp_leakcorr_fintcells{t,xy};
    tfp_leakcorr_fintview(t,xy) = mean(t_temp);
    tfp_leakcorr_sem(t,xy) = std(t_temp)/sqrt(length(t_temp)); 
    tfp_leakcorr_std(t,xy) = std(t_temp);
    
    clearvars r_temp y_temp t_temp
    
    end
end

% clear variables before save 
clearvars ans c cell f final_timepoint h heights LeakCorr q 
clearvars start_timepoint t xy 

save(strcat(analysis_kinetics,'\09_SpectralLeakageCorrector.mat'))

%% NEXT: runYFPpositiveFunc.m







