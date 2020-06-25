%% Identify YFP-positive cells (clean up data)
% version JULY 2018

% with increasingly stringent YFP threshold, does the average look better?

% 1) collect YFP intensities by t and q
% 2) histogram of YFP intensities of all cells
% 3) histogram of YFP intensities of cells separated by t and q
% 4) apply YFP threshold
% 5) number of cells that are positive or negative in YFP, plots
% 6) Plot: average fluorescence by q (post-YFP-threshold)

% analysis_kinetics = uigetdir; % choose the analysis/DMSPkinetics folder
% (should be the same if running after SpectralLeakageCorrector_3.m)
clearvars -except analysis_kinetics
cd(analysis_kinetics);
load 09_SpectralLeakageCorrector.mat
saveplace = strcat(analysis_kinetics,'/YFP positive'); mkdir(saveplace); % where plots will be saved


%% 1) Collect YFP intensities by q and t

% define q and xy positions containing 3cv3 or 3cv4 
q_with_yfp = analyzed_q;

% collect all yfp intensities of all cells
inputdata = yfp_leakcorr_fintcells; % pay attention to before or after leakage correction
yfp_allcells = [];

for q = analyzed_q;
    xy_range = expt_conditions_xy{q};
    yfp_allcells_q_temp = [];
    for t = 1:total_time;    
        for xy = xy_range;
            yfp_allcells = [yfp_allcells inputdata{t,xy}];
            yfp_allcells_q_temp = [yfp_allcells_q_temp inputdata{t,xy}]; % deleted after every q
        end
        yfp_fintcells_leakcorr_q{t,q} = yfp_allcells_q_temp; % store in permanent folder
    end  
end
clearvars yfp_allcells_q_temp q xy_range t

% save collected yfp intensities
save(strcat(analysis_kinetics,'\10_YFPallcells_afterLeakCorr.mat'),'yfp_allcells','yfp_fintcells_leakcorr_q'); % save variable


%% 2) Visualize distribution of YFP intensities of all cells
% make a decision about YFP threshold 

thresh_test = [25 50 100 150 200]; % plotted on histogram
SaveIm = 1; % 1 = save image and work space
% histogram and visualize threshold
yfp_cells = figure('units','normalized','outerposition',[0 0 1 1]);

% plot
histogram(yfp_allcells,1000);
    ylim_temp = ylim; % y axis limit
    hold on
    
% plot threshold
for i = 1:length(thresh_test);
    thresh_plot = thresh_test(i);
    plot(thresh_plot,0:ylim_temp(2)/50:ylim_temp(2),'r.','MarkerSize',10); % dot-plot of threshold
end

% set axis limits
xlim([-20 800]);

% figure labels
title({strcat('distribution of YFP fluorescence intensities AFTER SPECTRAL LEAKAGE CORRECTION, YFP threshold =',num2str(thresh_test)),...
       strcat('Experiment date: ',expt_date)},'FontSize',20);
ylabel('frequency of cells','FontSize',20);
xlabel('YFP fluorescence intensity after background subtraction','FontSize',20);
set(gca, 'FontSize', 20);

% save figure
if SaveIm == 1;
    savename = strcat('\\YFPcells_distribution_threshes','_afterLeakCorr.png');
    print(yfp_cells,strcat(saveplace,savename),'-dpng','-r300'); % resolution 300

    clearvars t xy inputdata ylim_temp ylim thresh_test SaveIm yfp_cells; close all
end



%% 3) plot YFP intensities by t and q
% OPTIONAL section, to visualize YFP distribution evolution over time for each q.
% Line histogram of each condition overlaid, at each time point.

saveplace_thist = strcat(analysis_kinetics,'/YFP positive/allq_t_slice'); mkdir(saveplace_thist);
inputdata = yfp_fintcells_leakcorr_q;
cmap = hsv(total_cond);
nbins = 30;
SaveIm = 1; % 1 = save images

if SaveIm == 0;
    t_range = [1 10 20 30]; % a sample of time range to see if axes are ok
elseif SaveIm == 1;
    t_range = 1:total_time;
end

for t = t_range;
    fig(t) = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    for q = analyzed_q;
        clearvars centers heights data; % initialize for each line
        data = inputdata{t,q}; 

        histdata_q(q) = histogram(data,nbins,'Normalizatio','probability','FaceColor',cmap(q,:));
            centers = histdata_q(q).BinEdges + histdata_q(q).BinWidth/2; % center of each bar in histogram
            heights = [histdata_q(q).Values,0];        
        hold on
        line_q(q) = plot(centers,heights,'LineWidth',3,'Color',cmap(q,:));
        set(histdata_q(q),'facealpha',0,'edgecolor','none'); % turn OFF bars
    end
    
    % for each figure (at each time point)
    legend(line_q(analyzed_q),expt_conditions_string,'FontSize',16); 
    xlabel('YFP intensity (after leakage correction)','FontSize',20);
    ylabel('probability','FontSize',20);
    title(strcat('t = ',num2str(t),'  ; Distribution of YFP intensities of all cells, separated by q; nbins = ',num2str(nbins),'; expt = ',expt_date),...
          'FontSize',20);
    xlim([-20 800])
    ylim([0 0.3])
    if SaveIm == 1; % save image at each time point
        savename = strcat('\\YFPintensity_distribution_allq_t',num2str(t),'_afterLeakCorr.png');
        print(fig(t),strcat(saveplace_thist,savename),'-dpng','-r150'); % resolution 150
        close all
    end
         
end
  
%% 3) clean workspace  
clearvars t q nbins fig inputdata centers heights saveplace_thist SaveIm t_range; hold off;
close all


%% 4) Apply YFP threshold
close all;
thresh_YFP_range = [25 50 100 150 200]; % all the YFP thresholds to run
% q, xy, and t to be analyzed
q_range_analyzed = q_with_yfp;
xy_range_analyzed = [expt_conditions_xy{q_range_analyzed}];%xy_with_yfp;
t_range_analyzed = 1:total_time;
skip_txy = [0,0]; % in [t,xy] array format

for i = 1:length(thresh_YFP_range);
    
    thresh_YFP = thresh_YFP_range(i);

    % fluorescence values to threshold on
    rfp_cells_fluo = rfp_leakcorr_fintcells;
    yfp_cells_fluo = yfp_leakcorr_fintcells;
    tfp_cells_fluo = tfp_leakcorr_fintcells;

    % run script to threshold YFP
    [savename] = YFPpositiveFunc(skip_txy, thresh_YFP, xy_range_analyzed, t_range_analyzed, q_range_analyzed, ...
                                 rfp_cells_fluo, yfp_cells_fluo, tfp_cells_fluo,...
                                 cell_num, cell_size, expt_conditions_xy, analysis_kinetics);

   
end 

% clean up workplace before loading YFP-thresholded data
clearvars saveplace xy_range_analyzed time_range_analyzed q_range_analyzed 
clearvars yfp_fintcells_ypos yfp_fintcells_yneg tfp_fintcells_ypos tfp_fintcells_yneg rfp_fintcells_ypos rfp_fintcells_yneg;
clearvars yfp_fintview_yneg yfp_fintview_ypos tfp_fintview_yneg tfp_fintview_ypos rfp_fintview_yneg rfp_fintview_ypos;
clearvars cell_num_ypos cell_num_yneg percent_ypos percent_yneg
clearvars cell_num_ypos_avg cell_num_ypos_std cell_num_yneg_avg cell_num_yneg_std
clearvars rfp_fintcells_q_ypos tfp_fintcells_q_ypos
clearvars rfp_mean_q_ypos tfp_mean_q_ypos rfp_std_q_ypos tfp_std_q_ypos rfp_sem_q_ypos tfp_sem_q_ypos ans

% load YFP-thresholded data
% savename = '11_YFPon_thresh_25.mat'
% load(savename);                          
                

%% 5) number of cells that are positive or negative in YFP
% OPTIONAL, should be done for each YFP threshold

for xy = 1:total_xy;
    for t = 1:total_time;
    cell_num_ypos(t,xy) = sum(~isnan(rfp_fintcells_ypos{t,xy}));
    cell_num_yneg(t,xy) = sum(~isnan(rfp_fintcells_yneg{t,xy}));    
    end
end
% calculate ratios at each xy
cell_num_ypos_yneg_ratio_xy = cell_num_ypos ./ cell_num_yneg; % calculate ratios

% number of cells by expt condition
for q = 1:total_cond;
    xy_range = expt_conditions_xy{q};
    for t = 1:total_time;        
        cell_num_ypos_avg(t,q) = mean(cell_num_ypos(t,xy_range));
        cell_num_ypos_std(t,q) = std(cell_num_ypos(t,xy_range));
        cell_num_yneg_avg(t,q) = mean(cell_num_yneg(t,xy_range));
        cell_num_yneg_std(t,q) = std(cell_num_yneg(t,xy_range));
        % take the mean ratio for each q
        cell_num_ypos_yneg_ratio_q(t,q) =  mean(cell_num_ypos_yneg_ratio_xy(t,xy_range));
        cell_num_ypos_yneg_ratio_q_std(t,q) =  std(cell_num_ypos_yneg_ratio_xy(t,xy_range));
    end
end
clearvars ans h q xy xy_range t

% quick plot to visualize Ypos and Yneg number of cells per xy
plot_cellnum_ypos = cell_num_ypos
plot_cellnum_yneg = cell_num_yneg
title_text = strcat('number of YFP-positive and YFP-negative cells for each xy; expt = ',expt_date)
x_text = 'time (hours)'
y_text = 'number of cells'
leg_text = {'YFP positive','YFP negative'};
[h xy_t_cellprop] = plot_xy_t_cellnum(time_expt_hours, q_range,...
                                      expt_conditions_string, expt_conditions_xy,...
                                      plot_cellnum_ypos, plot_cellnum_yneg,...
                                      title_text, x_text, y_text, leg_text);
                                                                 
% plot ratios Ypos / Yneg (cell numbers)
figure('units','normalized','outerposition',[0 0 1 1]);
for q = 1:9;
    errorbar(time_expt_hours(:,1),cell_num_ypos_yneg_ratio_q(:,q),cell_num_ypos_yneg_ratio_q_std(:,q),...
        'color',cmap(q,:),'linewidth',2); hold on;
end
legend({expt_conditions_string{:}},'FontSize',16);
title(strcat('ratio of number of YFP-positive cells and number of YFP-negative cells; errorbars = std across ratios at 7 xy positions; Ythresh = ',num2str(thresh_YFP),...
    '; expt = ',expt_date))
xlabel('time (hours)')
ylabel('ratio (Ypos / Yneg cell number)')
set(gca,'FontSize',20)
% save figure
savename = strcat('\\cellnum_ratio_ypos_yneg_afterLeakCorr.png');
print(strcat(saveplace,savename),'-dpng','-r300'); % resolution 300
        
        
% save variables
save('cellnum_ypos_yneg.mat',...
    'cell_num_ypos','cell_num_yneg','cell_num_ypos_avg','cell_num_ypos_std',...
    'cell_num_yneg_avg','cell_num_yneg_std','thresh_YFP','cell_num_ypos_yneg_ratio_xy',...
    'cell_num_ypos_yneg_ratio_q','cell_num_ypos_yneg_ratio_q_std');
                         

%% 6) Plot: average fluorescence by q (post-YFP-threshold)
% OPTIONAL, should be done for each YFP threshold.

avgFluo = figure('units','normalized','outerposition',[0 0 1 1]);

q_range = q_range_analyzed;
plotYFP = 1; % 1 = plot 2 colors; 0 = plot only 1 color

plot_rfp = rfp_mean_q_ypos;
error_rfp = rfp_sem_q_ypos;
plot_tfp = tfp_mean_q_ypos;
error_tfp = tfp_sem_q_ypos;

if plotYFP == 1;
    plot_yfp = yfp_mean_q_ypos;
    error_yfp = yfp_sem_q_ypos;
end

% for marker legend
leg = {'RFP (dddW)',...
       'TFP (dmdA)'}; % legend
if plotYFP == 1;
   leg{3} ='YFP (constitutive)';
end

% the real plot
for q = q_range;
    i = find(q == q_range);
    h(i) = subplot(1,length(q_range),i);
    %h(i) = subplot(1,6,i);
    title(expt_conditions_string{q});
    hold on;
    % RFP
        errorbar(time_expt_hours(:,1),plot_rfp(:,q),error_rfp(:,q),...
                     marker_rfp,'LineWidth', linewidth_rfp, ...
                     'Color', color_rfp_in, ...
                     'MarkerEdgeColor', color_rfp_out, ...
                     'MarkerFaceColor',color_rfp_in, ...
                     'MarkerSize', markersize_rfp);
    % TFP
        errorbar(time_expt_hours(:,1),plot_tfp(:,q),error_tfp(:,q),...
                     marker_tfp,'LineWidth', linewidth_tfp, ...
                     'Color', color_tfp_in, ...
                     'MarkerEdgeColor', color_tfp_out, ...
                     'MarkerFaceColor',color_tfp_in, ...
                     'MarkerSize', markersize_tfp);
    if plotYFP == 1;
    % YFP
        errorbar(time_expt_hours(:,1),plot_yfp(:,q),error_yfp(:,q),...
                     marker_yfp,'LineWidth', linewidth_yfp, ...
                     'Color', color_yfp_in, ...
                     'MarkerEdgeColor', color_yfp_out, ...
                     'MarkerFaceColor',color_yfp_in, ...
                     'MarkerSize', markersize_yfp); 
    end
    
   % axis([0 15.3 0 18]);
    set(gca, 'FontSize', 10)
end

suptitle({sprintf('DMSP kinetics (errorbars = SEM); YFP threshold after leakage correction = %0.f',thresh_YFP),...
       strcat('Experiment date: ',expt_date)});
      
xlabel(h(5),'time from treatment (hours)','FontSize',16)
ylabel(h(1),'average fluorescence intensity per cell, background subtracted (A.U.)','FontSize',16);

% Generate dummy info for plots with handles "m" outside of range of last subplot
m = zeros(2,1);
m(1) = plot(-1,-1, marker_rfp,'LineWidth', linewidth_rfp, ...
                 'Color', color_rfp_in, ...
                 'MarkerEdgeColor', color_rfp_out, ...
                 'MarkerFaceColor',color_rfp_in, ...
                 'MarkerSize', markersize_rfp); hold on;
m(2) = plot(-1,-1,marker_tfp,'LineWidth', linewidth_tfp, ...
                 'Color', color_tfp_in, ...
                 'MarkerEdgeColor', color_tfp_out, ...
                 'MarkerFaceColor',color_tfp_in, ...
                 'MarkerSize', markersize_tfp+2);
if plotYFP == 1;
m(3) = zeros;
m(3) = plot(-1,-1,marker_yfp,'LineWidth', linewidth_yfp, ...
                 'Color', color_yfp_in, ...
                 'MarkerEdgeColor', color_yfp_out, ...
                 'MarkerFaceColor',color_yfp_in, ...
                 'MarkerSize', markersize_yfp); hold off;    
end
  
legend(m(:),leg,'FontSize',16);

%% axis limits
% axis(h(1),[0 26.5 0 35]);
xlim(h(:),[0 24]);
ylim(h(:),[0 320]);


%%
% saveplace = strcat(analysis_kineticsfolder,'/data plots'); mkdir(saveplace);

% save figure
print(avgFluo,strcat(saveplace,'\avgFluo_YFPpos_YFPthresh_',num2str(thresh_YFP),'.png'),'-dpng','-r300'); % resolution 300
% print(avgFluo,strcat(saveplace,'\avgFluo_YFPpos_YFPthresh_TFPonly',num2str(thresh_YFP),'_final.png'),'-dpng','-r300'); % resolution 300
 
% clearvars ans i h leg m q t xy avgFluo

%% NEXT: PlotGenerator.m



