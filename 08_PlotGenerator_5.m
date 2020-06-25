%% Plot Generator
% version March 2018
% UPDATED: 12/12/2019 [code from MultipleExpAggregated_PlotGenerator_2.m]

% 1) calculate time points
% 2) specify colors and plot styles (eliminate in the future)
% 3) plot xy vs. time - generic
% 4) plot q vs. time - generic
% 5) cell properties: # or % cells YFP-positive OR cell size
% 6) growth curves
% 7) YFP vs. YFP-positive cell size vs. YFP-positive cell number 
% 8) FP vs. FP scatter

% dot plot vs. time (each dot = cell)
% 8) make a plot with shaded error bars


% load analysis results
analysis = uigetdir; % choose folder with analysis results data
load(strcat(analysis,'/13_YFP_normalized_t2Y_Ythresh_50.mat'));
load(strcat(analysis,'/time_expt.mat'));
saveplace = uigetdir; % choose folder to save figures in 

% load publication colormap (Parula, 12/2018)
load('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/Matlab/plot color map/colormap_publication.mat')

expt_date = '11/16/2017 (3cv3 no Ab DMSP kinetics; for reviewer response)'; % for titles


%% 01) 2-panel figures for each experiment (dmdA and dddW) 
% all data are background-subtracted, spectral leakage-corrected, YFP-thresholded, YFP-normalized

close all;
SaveIm = 0; 

% define data and labels
leg_string = expt_conditions_string;
fp_color{1} = 'TFP'; % color of demethylation
fp_color{2} = 'RFP'; % color of cleavage
data_plot_A = tfp_fintview_q_leakcorr_ypos_ynorm;
data_plot_W = rfp_fintview_q_leakcorr_ypos_ynorm;
error_plot_A = tfp_sem_q_leakcorr_ypos_ynorm;
error_plot_W = rfp_sem_q_leakcorr_ypos_ynorm;
t_range = time_expt_hours(:,1);

% generate a figure
if SaveIm == 0
    fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'on'); % figure pops up
elseif SaveIm == 1
    fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off'); % no figure popup   
end

% plot (each subplot = one pathway)
for q = 1:9 % iterate through DMSP concentrations
    s(1) = subplot(121); % demethylation
        errorbar(t_range, data_plot_A(:,q), error_plot_A(:,q), '.-', 'markersize', 20, 'color', cmap_final(q,:)); hold on;       

    s(2) = subplot(122); % cleavage
        errorbar(t_range, data_plot_W(:,q), error_plot_W(:,q), '.-', 'markersize', 20, 'color', cmap_final(q,:)); hold on; 
end

% labels
xlabel(s(1), 'time since treatment (hours)'); 
xlabel(s(2), 'time since treatment (hours)'); 
ylabel(s(1), strcat('fluorescence (',fp_color{1},')'))
ylabel(s(2), strcat('fluorescence (',fp_color{2},')'))
title(s(1),'dmdA (demethylation)')
title(s(2),'dddW (cleavage)')
sgtitle({expt_date,'fluorescence is background-subtracted, NOT spectral leakage-corrected, YFP-thresholded, and YFP-normalized (by t2 average)'})
legend(s(1),leg_string,'FontSize',15,'location','northwest')
ax = findall(gcf,'type','axes'); set(ax,'fontsize', 20); % set all axis labels the same size

if SaveIm == 1
    savename = strcat('/t_vs_fluo_2panel_noleakcorr');
    print(fig,strcat(saveplace,savename),'-depsc'); % resolution 300
end
    



%% BELOW: old
% above = from MultipleExpAggregated_PlotGenerator.m

%% calculate time points
% for // expt

% start_timepoint = 15; % minutes after treatment
% time_interval = 45; % minutes between time points
% final_timepoint = start_timepoint + (time_interval * (length(bfImgs(:,1)) - 1 ) ); % total minutes after treatment
% time_expt = [start_timepoint : time_interval : final_timepoint] % calculate time points (should be same size as avgFluo)
% time_expt_hours = time_expt ./ 60 % in hours
% time_expt_hours = repmat(time_expt_hours,total_xy_per_cond,1)' % for plotting # xy positions per condition
% 
% clearvars start_timepoint time_interval final_timepoint time_expt

% save analysis results workspace
% save(strcat(analysis,'\06_FluoCellComputation_analysisonly.mat')); % overwrite previous version


%% specify colors (RGB) and plot styles
% http://shirt-ediss.me/matlab-octave-more-colours/

% color_rfp_in = [246 175 66] ./ 255;
% color_rfp_out = [231 99 12] ./ 255;
% color_yfp_in = [94 250 81] ./ 255;
% color_yfp_out = [12 195 82] ./ 255;
% color_tfp_in = [8 180 238] ./ 255;
% color_tfp_out = [1 17 181] ./ 255;
% 
% 
% marker_rfp = 'gs';
% marker_yfp = 'go';
% marker_tfp = 'gd';
% 
% markersize_rfp = 9;
% markersize_yfp = 7;
% markersize_tfp = 5;
% 
% linewidth_rfp = 2;
% linewidth_yfp = 2;
% linewidth_tfp = 1.5;
% 
% plotstyles = {color_rfp_in, color_rfp_out, color_yfp_in, color_yfp_out, color_tfp_in, color_tfp_out,...
%     marker_rfp, marker_yfp, marker_tfp,...
%     markersize_rfp, markersize_yfp, markersize_tfp,...
%     linewidth_rfp, linewidth_yfp, linewidth_tfp}

%% 3) plot xy vs. time - generic
% For scenario where there are multiple curves (for multiple xy positions) per expt condition.
    % i.e. 7 curves per color per subplot

PlotYFP = 1; % 1 = plot YFP data; 0 = no YFP data plotted but still need to give dummy data
ErrorBarOn = 1; % 1 = with errorbar; 0 = no errorbar, but still need to give dummy error data
SaveIm = 0; 

q_range = 1:9;%q_range_analyzed; % will determine xy data from expt_conditions_xy

% define data to plot (multiple data array per expt condition)
plot_rfp =  rfp_fintview;
plot_tfp = tfp_fintview;
plot_yfp = yfp_fintview;
error_rfp = rfp_sem;
error_tfp =tfp_sem;
error_yfp = yfp_sem;

% plot labels
title_text = {sprintf('3cv3 DMSP kinetics (leakage corrected; YFP threshold = %.f); each point = avg fluo of cells in image',thresh_YFP)...
              strcat('errorbars = SEM;     expt date:   ',expt_date,')')};
x_text = 'time from treatment (hours)';
y_text = 'fluorescence (A.U.)';
leg_text = {'RFP (dmdA)',...
            'TFP (dddW)',...
            'YFP (constitutive)'}; % legend, in the order of rfp -> tfp ( -> yfp )

[h, xy_t] = plot_xy_t(time_expt_hours, PlotYFP, ErrorBarOn, q_range,...
                     expt_conditions_string, expt_conditions_xy,...
                     plot_rfp, plot_tfp, plot_yfp, error_rfp, error_tfp, error_yfp,...
                     title_text, x_text, y_text, leg_text);

                                         
% alter the axis range
xlim(h(:),[0 24]);  
%ylim(h(:),[300 370]);
                                         
% save figure
if SaveIm == 1;
    fig_savename = strcat('\avgFluo_xy_DMSPkinetics_YFPthresh_',num2str(thresh_YFP),'.png')
    print(xy_t,strcat(saveplace,fig_savename),'-dpng','-r300'); % resolution 300
end 

clearvars PlotYFP ErrorBarOn plot_rfp plot_tfp plot_yfp error_rfp error_tfp error_yfp
clearvars title_text x_text y_text leg_text fig_savename SaveIm

%% 3) continued - investigate outlier txy positions
value = 126.5; % from datatip
dec_place = 1; % number of decimal places
q = 7; % which subplot panel
t = 1:total_time;
temp = round(yfp_fintview_ypos, dec_place); % round to number of dec places in data tip so that equivalency can be determined
[t_out, xy_out] = find(temp(t,expt_conditions_xy{q})==value)
xy_out = expt_conditions_xy{q}(xy_out)

% investigate why this is an outlier
figure('units','normalized','outerposition',[0 0 1 1]);
phasename = bfImgs(t_out,xy_out).name;
phase = imread(phasename);
ph_adj = imadjust(phase);

yname = YFP(t_out,xy_out).name;
y_im = imread(yname);
y_adj = imadjust(y_im);
ax(1) = subplot(1,3,1); imshowpair(ph_adj,bwf_post_thresh{t_out,xy_out}); title('overlay');
ax(2) = subplot(1,3,2); imshowpair(y_adj,bwf_post_thresh{t_out,xy_out}); title('overlay');
ax(3) = subplot(1,3,3); imshow(y_adj); title(sprintf('phase original; (%s)',phasename));
linkaxes(ax,'xy');

%% aggregate

% suspicious_txy = [];

suspicious_txy = [suspicious_txy; t_out, xy_out];

% clearvars t_out xy_out

%% 4) plot q vs. time - generic
% For scenario where there is one curve per expt condition.
    % i.e. 1 curve per color per subplot

PlotYFP = 1; % 1 = plot YFP data; 0 = no YFP data plotted but still need to give dummy data
ErrorBarOn = 1; % 1 = with errorbar; 0 = no errorbar, but still need to give dummy error data

q_range = 1:total_cond; % will determine xy data from expt_conditions_xy

% subtract values of negative C control from each time point
rfp_fintview_q_ypos_negsub = []; tfp_fintview_q_ypos_negsub = [];
for t = 1:total_time;
    for q = 2:total_cond;
        rfp_fintview_q_leakcorr_ypos_negsub(t,q) = rfp_mean_q_ypos(t,q) - rfp_mean_q_ypos(t,1);
        tfp_fintview_q_leakcorr_ypos_negsub(t,q) = tfp_mean_q_ypos(t,q) - tfp_mean_q_ypos(t,1);
    end
end
rfp_fintview_q_leakcorr_ypos_negsub(:,1) = rfp_mean_q_ypos(:,1);
tfp_fintview_q_leakcorr_ypos_negsub(:,1) = tfp_mean_q_ypos(:,1);

% take the ratio of dddW/dmdA at each timepoint
% for t = 1:total_time;
%     for q = 2:7;
%         ratio(t,q) = rfp_fintview_q_ypos_ynorm_negsub(t,q)./tfp_fintview_q_ypos_ynorm_negsub(t,q);
%     end
% end


% define data to plot (1 data array per expt condition)
plot_rfp = rfp_fintview_q_leakcorr_ypos_negsub;% ratio;
plot_tfp = tfp_fintview_q_leakcorr_ypos_negsub;%tfp_fintview_q_ypos_ynorm_negsub*36;
plot_yfp = yfp_mean_q_ypos;
error_rfp = rfp_sem_q_ypos;
error_tfp = tfp_sem_q_ypos;
error_yfp = yfp_sem_q_ypos;


% plot labels
% title_text = {sprintf('RATIO dddW / dmdA (normalized by YFP of each cell; subtracted succinate data at each time point (YFP threshold = %.f)',thresh_YFP)...
%               'TFP values NOT multiplied by 22 for ratio calculation',...
%                strcat('expt date:   ',expt_date,')')};

title_text = {sprintf('average fluorescence after leakage correction, YFP thresholding, subtraction of negative C control at each time point (YFP threshold = %.f)',thresh_YFP)...
             % 'TFP values NOT multiplied by 36 for ratio calculation',...
               strcat('each point = average fluo of pooled cells;  errorbars = SEM;  expt date:   ',expt_date)};
           
% title_text = {sprintf('average normalized fluorescence (YFP threshold = %.f; each cell normalized by its own YFP value); each point = avg fluo of pooled cells',thresh_YFP)...
%              % 'each data point subtracted by corresponding data point from succinate condition; all TFP values multiplied by 22',...
%                strcat('errorbars = SEM (may be too small to be seen);     expt date:   ',expt_date)};
x_text = 'time from treatment (hours)';
y_text = 'average fluorescence (normalized by YFP) (A.U.)';
% y_text = 'RATIO RFP/TFP (dddW/dmdA)';
leg_text = {'RFP (dmdA)',...
            'TFP (dddW)',...
            'YFP (constitutive)'}; % legend, in the order of rfp -> tfp ( -> yfp )
        

[h q_t] = plot_q_t(time_expt_hours, PlotYFP, ErrorBarOn, q_range,...
                   expt_conditions_string, expt_conditions_xy,...
                   plot_rfp, plot_tfp, plot_yfp, error_rfp, error_tfp, error_yfp,...
                   title_text, x_text, y_text, leg_text);

axis(h(:),[0 24 -10 400]);   

% save figure
fig_savename = '\avgFluo_RFPTFP_leakcorr_Ypos_negsub.png';
print(q_t,strcat(saveplace,fig_savename),'-dpng','-r300'); % resolution 300
                            
clearvars PlotYFP ErrorBarOn plot_rfp plot_tfp plot_yfp error_rfp error_tfp error_yfp
clearvars title_text x_text y_text leg_text fig_savename


%% 5) CELL PROPERTY: # or % cells YFP positive OR cell size
% [need to fix]
close all

q_range = q_range_analyzed; % will determine xy data from expt_conditions_xy
plot_number = 0; % plot cell number YFP positive or negative
plot_percent = 0; % plot percentage YFP positive or negative
plot_size = 1;
SaveIm = 0; 


% define data to plot (1 data array per expt condition)
if plot_number == 1;
    plot_percent = 0; plot_size = 0; % set other flags to 0
        plot_cellprop_ypos = cell_num_ypos; % cell # ypos
        plot_cellprop_yneg = cell_num_yneg; % cell # yneg
    title_text = {strcat('YFP-positive and YFP-negative NUMBER of cells (leakage corrected; YFP-positive threshold = ',num2str(thresh_YFP),'); each point = one image'),...
              expt_date};
    x_text = 'time from treatment (hours)';
    y_text = '# of cells';
    fig_savename = strcat('\cellNum_xy_YFPpos_YFPneg_YFPthresh_',num2str(thresh_YFP),'.png')
    axislim = [0 24 0 800];
elseif plot_percent == 1;   
    plot_number = 0; plot_size = 0; % set other flags to 0
        plot_cellprop_ypos = percent_ypos; % percent_ypos
        plot_cellprop_yneg = percent_yneg; % percent_yneg
    title_text = {strcat('YFP-positive and YFP-negative PERCENTAGE of total cells (leakage corrected; YFP-positive threshold = ',num2str(thresh_YFP),'); each point = one image'),...
              expt_date};
    x_text = 'time from treatment (hours)';
    y_text = '% of cells';        
    fig_savename = strcat('\cellPercent_xy_YFPpos_YFPneg_YFPthresh_',num2str(thresh_YFP),'.png')
    axislim = [0 24 0 110];
elseif plot_size == 1;
    plot_number = 0; plot_percent = 0;
        plot_cellprop_ypos = cell_size_ypos_avg; % cell size ypos
        plot_cellprop_yneg = cell_size_yneg_avg; % cell size yneg
    title_text = {strcat('YFP-positive and YFP-negative SIZE of cells (leakage corrected; YFP-positive threshold = ',num2str(thresh_YFP),'); each point = one image'),...
              expt_date};
    x_text = 'time from treatment (hours)';
    y_text = 'cell size';
    fig_savename = strcat('\cellSize_xy_YFPpos_YFPneg_YFPthresh_',num2str(thresh_YFP),'.png')
    axislim = [0 24 0 110];
end


leg_text = {'YFP-positive','YFP-negative'};
        

[h xy_t_cellprop] = plot_xy_t_cellprop(time_expt_hours, q_range,...
                                     expt_conditions_string, expt_conditions_xy,...
                                     plot_cellprop_ypos, plot_cellprop_yneg,...
                                     title_text, x_text, y_text, leg_text);
                              

                                 
axis(h(:),axislim);

% save figure
if SaveIm == 1;
    print(xy_t_cellprop,strcat(saveplace,fig_savename),'-dpng','-r300'); % resolution 300
end

% clearvars ans h leg m q t xy avgFluo cellNumber


%% 6) growth curves
figure('units','normalized','outerposition',[0 0 1 1]);
cmap = hsv(total_cond);

for q = 1:total_cond;
    errorbar(time_expt_hours(:,1),cell_num_avg(:,q),cell_num_std(:,q),'.-','MarkerSize',15,'LineWidth',3,'Color',cmap(q,:))
    hold on
end
legend({expt_conditions_string{1:total_cond}},'FontSize',16,'Location','NorthWest');
title('Growth Curves (average of total number of cells recognized in an image)')
ylabel('cell number'); xlabel('time, hours since treatment');
set(gca,'FontSize',16)
clearvars cmap


%% 7) YFP vs. YFP-positive cell size vs. YFP-positive cell number
close all;
YvsCellProp = figure('units','normalized','outerposition',[0 0 1 1]);

q_range = q_range_analyzed;

YvsSize = 0;
YvsNum = 0;
SizevsNum = 0;
YvsYnegNum = 1; 

if YvsSize == 1;
    YvsNum = 0; SizevsNum = 0; YvsYnegNum = 0;
    title_text = {sprintf('YFP vs. cell SIZE (YFP threshold after leakage correction = %0.f',thresh_YFP),...
                 strcat('Experiment date: ',expt_date)};
        plot_axis1 = yfp_mean_q_ypos;
            ytext_axis1 = 'average YFP fluorescence intensity per cell, background subtracted (A.U.)';    
        plot_axis2 = cell_size_ypos_avg;
            ytext_axis2 = 'average cell SIZE per image (YFP-positive only)';
elseif YvsNum == 1;
    YvsSize = 0; SizevsNum = 0; YvsYnegNum = 0;
    title_text = {sprintf('YFP vs. cell NUMBER (YFP threshold after leakage correction = %0.f',thresh_YFP),...
                 strcat('Experiment date: ',expt_date)};
        plot_axis1 = yfp_mean_q_ypos;
            ytext_axis1 = 'average YFP fluorescence intensity per cell, background subtracted (A.U.)';    
        plot_axis2 = cell_num_ypos;
            ytext_axis2 = 'average cell NUMBER per image (YFP-positive only)';
elseif SizevsNum == 1;
    YvsSize = 0; YvsNum = 0; YvsYnegNum = 0;
    title_text = {sprintf('cell SIZE vs. cell NUMBER (YFP threshold after leakage correction = %0.f',thresh_YFP),...
                 strcat('Experiment date: ',expt_date)};
        plot_axis1 = cell_size_ypos_avg;
            ytext_axis1 = 'average cell SIZE per image (YFP-positive only)';  
        plot_axis2 = cell_num_ypos;
            ytext_axis2 = 'average cell NUMBER per image (YFP-positive only)';
elseif YvsYnegNum == 1;
    YvsSize = 0; YvsNum = 0; SizevsNum = 0;
    title_text = {sprintf('YFP vs. YFP-negative cell NUMBER (due to YFP threshold after leakage correction = %0.f',thresh_YFP),...
                 strcat('Experiment date: ',expt_date)};
        plot_axis1 = yfp_mean_q_ypos;
            ytext_axis1 = 'YFP intensity (YFP-positive only)';  
        plot_axis2 = cellnum_yneg_q;
            ytext_axis2 = 'average cell NUMBER per image (eliminated due to YFP threshold)';
end

    
    


% yyaxis left; ylim(h(:),[200 380]);
% yyaxis right; ylim(h(:),[30 60]);
   
%   fig_savename = strcat('\CellNumvsCellSize_YFPpos_YFPthresh_',num2str(thresh_YFP),'.png')
%   print(YvsCellProp,strcat(saveplace,fig_savename),'-dpng','-r300'); % resolution 300



% the real plot
for q = q_range;
    i = find(q == q_range);
    h(i) = subplot(1,length(q_range),i);
    title(expt_conditions_string{q});
    hold on;
      
    % YFP
        yyaxis left        
       % ylim([200 380]); 
        plot(time_expt_hours(2:end,1),plot_axis1(2:end,q),...%error_yfp(:,q),...
                     '.-','LineWidth', linewidth_yfp)%, ...
             %        'Color', color_yfp_in, ...
             %        'MarkerEdgeColor', color_yfp_out, ...
             %        'MarkerFaceColor',color_yfp_in, ...
             %        'MarkerSize', markersize_yfp); 
             
             if i == 1;
                 ylabel(ytext_axis1,'FontSize',16);
             end

    % cell size
        yyaxis right
       %  ylim([30 52]);        
        plot(time_expt_hours(2:end,1),plot_axis2(2:end,q),...%error_cellsize(:,q),...
                     '.-','LineWidth', linewidth_rfp)%, ...
             %        'Color', color_rfp_in, ...
              %       'MarkerEdgeColor', color_rfp_out, ...
              %       'MarkerFaceColor',color_rfp_in, ...
              %       'MarkerSize', markersize_rfp); 
                
            if i == 9;
                ylabel(ytext_axis2,'FontSize',16);
            end
%    set(gca, 'FontSize', 12)
end
xlim(h(:),[0 24])
xlabel(h(5),'time from treatment (hours)','FontSize',16)
suptitle(title_text);


%% added 3/20/2018
% histogram of fluorescence intensity distributions, separated by time
    % time as slices

% each q will be a separate figure
t_range = 1:total_time;
q_range = [9];
color = hsv(length(t_range));
SaveIm = 1; 

% histogram properties
nbin = 100; % number of bins in histogram
histplot_data = yfp_fintcells_leakcorr_q;% yfp_fintcells_q_ypos; % data to plot
normalization_method = 'probability';
x_text = 'fluorescence intensity';
y_text = 'probability';
leg_text = arrayfun(@num2str, t_range, 'UniformOutput', false);

for q = q_range;
%    HistInt = figure('units','normalized','outerposition',[0 0 1 1]);
    for t = t_range;
        HistInt = figure('units','normalized','outerposition',[0 0 0.5 0.9]); % square shape
       set(HistInt,'Visible','off')
        % i = find(t_range == t); % for defining bar plot color
       i = 1;
        histogram(histplot_data{t,q},nbin,...
                  'Normalization',normalization_method,...
                  'FaceColor',color(i,:),...
                  'FaceAlpha',0.4);
        hold on;
   % end % end of time
    
    % plot labeles
%    title_text = {strcat('Probability distribution of YFP intensities at different time points; BEFORE YFP threshold  expt condition: ',...
%                    expt_conditions_string{q}),expt_date,strcat('nbins = ',num2str(nbin))};%...
                    %strcat('YFP threshold = ',num2str(thresh_YFP))};
     title_text = {strcat('t',num2str(t)),expt_conditions_string{q}}
    xlabel(x_text,'FontSize',20); ylabel(y_text,'FontSize',20); suptitle(title_text);
%    legend(leg_text,'FontSize',20); % time points as legend
    hold off;
 
    % save image
    if SaveIm == 1;
        axis([-50 1000 0 0.06]);
        fig_savename = strcat('/q',num2str(q),'t',num2str(t),'_histogram_YFPintensityProbabilityDistribution_beforeYthresh.png');
        print(HistInt,strcat(saveplace,fig_savename),'-dpng','-r300'); % resolution 300
        close all; clearvars HistInt
    end
    end % alternative end of time
end
clearvars HistInt q t i nbin t_range q_range

% xlim([0 1000]);


%% look at the number of cells that are eliminated due to YFP threshold over time for all q
% 3/20/2018
% count the number in yfp_fintcells_yneg

for q = 1:total_cond;
    xy_range = expt_conditions_xy{q};
    for t = 1:total_time;
        for xy = xy_range;
            cellnum_yneg_q(t,q) = length(yfp_fintcells_yneg{t,xy});            
        end      
    end
end

% plot number of YFP-negative cells over time
YnegCellNum = figure('units','normalized','outerposition',[0 0 1 1]);
cmap = hsv(total_cond);
for q = 1:total_cond;
    plot(time_expt_hours(:,1),cellnum_yneg_q(:,q),'o-','MarkerSize',15,'Color',cmap(q,:),'LineWidth',2);
    hold on;
end
legend({expt_conditions_string{:,:}}, 'FontSize',20)
suptitle({'Number of cells eliminated due to YFP threshold = 100', expt_date});
xlabel('time since treatment (hours)','FontSize',20);
ylabel('number of cells eliminated due to YFP threshold','FontSize',20);

clearvars YnegCellNum q xy_range xy t cmap

%% 8) FP vs. FP scatter plot (dots = cells)
% Separate time by colors.

fitPerform = 0; % 1= to perform fit for each time point ('poly1') ==> NEED TO FIX LEGEND
% define x and y data points (e.g. RFP and TFP of each cell)
    % input data should be all cell data, pooled for each q
% x_data = tfp_fintcells_q_leakcorr_ypos_ynorm;
% y_data = rfp_fintcells_q_leakcorr_ypos_ynorm;
x_data = tfp_fintcells_q_ypos;
y_data = rfp_fintcells_q_ypos;

q_range = 2; % define the expt condition to plot (all cells in the xy condition will be 

time_range = [22:32]; 
color = jet(length(time_range));

markersize = 15;
title_text = {strcat(expt_conditions_string{q},'  (YFP positive cells only; YFP threshold = ',num2str(thresh_YFP),...
    '; different colors = time points;'),expt_date};


% title_text = {expt_conditions_string{q}, strcat('TFP vs. RFP of each cell (YFP-positive threshold = ',num2str(thresh_YFP),'); each point = one cell'),...
%               expt_date};
          
x_text = 'TFP';
y_text = 'RFP'

% SPECIFY LEGEND (NEED FIX)
if fitPerform == 1;
    leg_text = num2str(time_range(:)); % need to fix
else
    leg_text = num2str(time_range(:));
end

[scatterPlot] = scatter_2D(q_range, time_range, x_data, y_data, color, markersize, fitPerform,...
                                           title_text, x_text, y_text, leg_text);
         

% alternatively
figure;
for t = 22:32;
    for q = 2;
    x_data = tfp_fintcells_q_ypos{t,q};
    y_data = rfp_fintcells_q_ypos{t,q};

    scatter(x_data,y_data,'MarkerFaceColor','k','MarkerFaceAlpha',.05,'MarkerEdgeColor','none'); hold on;

    end
end
xlim([-10 100])
xlim([-10 600])
ylim([-10 600])
ylim([100 1000])
ylim([-10 100])
%ylim([-10 50])

suptitle({'all cells that are recognized between timepoints 22-32 in 1 mM DMSP condition (each dot = one cell at one timepoint)',...
    'spectral leakage corrected; YFP threshold = 100; marker alpha = 0.05',...
    '3/3/2018 (3cv4 [DMSP])'}) %expt_date})

% suptitle({'all cells that are tracked between timepoints 22-32 (each dot = one cell at one timepoint)',...
%     'spectral leakage corrected; YFP threshold = 100 for DMSP & gluc; YFP threshold = 50 for PA-YFP; marker alpha = 0.05',...
%     'red o = negative control cells (1 mM glucose) at t = 32;    blue o = PA-YFP leak corrected signal in RFP and TFP channels;         green o = 1 mM DMSP cell',...
%    '3/3/2018 (3cv4 [DMSP])'}) %expt_date})
xlabel('dmdA, RFP')
xlabel('dddW, TFP')
ylabel('dddW, TFP')
ylabel('YFP')
set(gca,'FontSize',20)


%% 7) dot plot vs. time (each dot = cell)
% plot all cells at each time point
% like a dot plot at each time point

allCells = figure('units','normalized','outerposition',[0 0 1 1]);
time_range = time_expt_hours(:,1);



plot_data = tfp_fintcells_q_ypos;
plot_data_mean_= tfp_mean_q_ypos;
plot_error = tfp_std_q_ypos;

plot_data_r = rfp_fintcells_q_ypos;
plot_data_mean_r = rfp_mean_q_ypos;
plot_error_r = rfp_std_q_ypos;

         plot_data = plot_data_r;
            plot_data_mean = plot_data_mean_r;
            plot_error = plot_error_r;
            marker_cell = 'r.';
            marker_mean = 'go';

            
%%            
% data to feed into function
q_range = 1:9;%1:total_cond;
t_range = 1:32;%1:total_time;
PlotMean = 0; % plot mean or not
ContourPlot = 1; % to plot contour (colormap) or not; 0 = scatter plot
plot_data = tfp_fintcells_q_ypos; %tfp_fintcells_q_leakcorr_ypos_ynorm;
plot_data_mean= 0;%tfp_mean_q_ypos;
plot_error = 0;% tfp_std_q_ypos;
            marker_cell = 'b.';
            marker_mean = 'go';
            markersize_cell = 1;
            markersize_mean = 3;
% title_text = {strcat('TFP fluorescence (dddW pathway) of each cell; each dot = 1 cell; errorbar = std; (YFP threshold = ',num2str(thresh_YFP),')'),...
%     expt_date}; % suptitle
title_text = {'TFP fluorescence (dddW pathway) of each cell; colorbar indicates number of cells in each fluorescence intensity bin (YFP threshold = 100)',expt_date}; % suptitle

x_text = 'time since treatment (hours)';
y_text = 'fluorescence intensity (A.U.), background subtracted, normalized by YFP';
leg_text = {'test'}; 

[scatter_t_fig, h] = scatter_t(q_range, t_range, expt_conditions_string, ContourPlot, PlotMean, ...
                                        plot_data, plot_data_mean, plot_error, ...
                                        marker_cell, marker_mean, markersize_cell, markersize_mean, ...              
                                        title_text, x_text, y_text, leg_text); 
                           
% for some reason this script leads to different sized subplots
for i = 1:9;
    temppos = get(h(i),'Position');
    temppos(3) = 0.0672; % width, for 9-panel subplot
    temppos(4) = 0.8150; % height, for 9-panel subplot
     set(h(i),'Position',temppos);
end

suptitle(title_text); % readjusts the height of the figure
  %%              contour plots           
axis(h(:),[0 32 -50 900]);
% axis(h(:),[0 32 -10 120]);  
axis(h(:),[0 32 -0.5 3]);
axis(h(:),[0 32 -0.5 0.5]);
% print(scatter_t_fig,strcat(saveplace,'\TFP_CellFluoInt_Ynorm_Ythresh100.png'),'-dpng','-r300'); % resolution 300

% change color bar limit for all
for i = 1:9
    h(i).CLim = [0 300];
end
xlim(h(:),[0 32]);
ylim(h(:),[-20 inf]);
     ylim(h(:),[-20 500]);
     
     
     
axis(h(:),[0 15.3 -20 3000])
axis(h(2),[0 15.3 -10 2500])
axis(h(3),[0 15.3 -10 3000])

suptitle({'dddW pathway: fluorescence value of all cells YFP-positive (thresh = 100)',...
        'each dot = fluo value of one cell (red = RFP; blue = TFP); overlaid points = mean with std  (11/16/2017 expt)'});

xlabel(h(2),'time from treatment (hours)','FontSize',16)
ylabel(h(1),'fluorescence value of each cell, background subtracted (A.U.)','FontSize',16);


% print(allCells,strcat(saveplace,'\RFPcell_Ypos_nomean_thresh100.png'),'-dpng','-r300'); % resolution 300

%% tease out the 2 populations that emerge in dddW

t_range = 13:32;
q = 2;
thresh_TFP = 50;
index_tneg = {}; index_tpos = {};
r_tneg_cum = []; t_tneg_cum = []; y_tneg_cum = [];
for t = t_range;
    
    index_tneg{t,q} = find(tfp_fintcells_q_ypos{t,2} < thresh_TFP);
    index_tpos{t,q} = find(tfp_fintcells_q_ypos{t,2} >= thresh_TFP);

    
    r_tneg{t,q} = rfp_fintcells_q_ypos{t,q}(index_tneg{t,q});
    t_tneg{t,q} = tfp_fintcells_q_ypos{t,q}(index_tneg{t,q});
    y_tneg{t,q} = yfp_fintcells_q_ypos{t,q}(index_tneg{t,q});
    
    r_tpos{t,q} = rfp_fintcells_q_ypos{t,q}(index_tpos{t,q});
    t_tpos{t,q} = tfp_fintcells_q_ypos{t,q}(index_tpos{t,q});
    y_tpos{t,q} = yfp_fintcells_q_ypos{t,q}(index_tpos{t,q});
    
    t_tneg_ynorm{t,q} = t_tneg{t,q} ./ y_tneg{t,q};
    r_tneg_ynorm{t,q} = r_tneg{t,q} ./ y_tneg{t,q};
    
    r_tneg_cum = [r_tneg_cum r_tneg{t,q}]; 
    t_tneg_cum = [t_tneg_cum t_tneg{t,q}];
    y_tneg_cum = [y_tneg_cum y_tneg{t,q}];
end


figure;
for t = t_range;
    plot(rfp_fintcells_q_ypos{t,q},yfp_fintcells_q_ypos{t,q},'gs'); hold on;
    plot(r_tneg{t,q},y_tneg{t,q},'r.');
    hold on
 %   plot(t_alt,r_alt,'b.')

  yyaxis right
  
    plot(tfp_fintcells_q_ypos{t,q},yfp_fintcells_q_ypos{t,q},'ms')
end

% box plot
x = [r_tneg{t_range,q}];
y = [r_tpos{t_range,q}];
y = y(~isnan(y));

x2 = [y_tneg{t_range,q}];
y2 = [y_tpos{t_range,q}];
y2 = y2(~isnan(y));

% z = rand(15,1);
% group = [repmat({'RFP (cells TFP < 50)'}, length(x), 1); repmat({'RFP (cells TFP >= 50)'}, length(y), 1),...
%     repmat({'YFP (cells TFP < 50)'}, length(x2), 1); repmat({'YFP (cells TFP >= 50)'}, length(y2), 1)];%; repmat({'Third'}, 15, 1)];
% boxplot([x;y;x2;y2], group)

G = [ones(size(x))  2*ones(size(y)) 3*ones(size(x2)) 4*ones(size(y2))];

X = [x, y, x2, y2];
boxplot(X,G,'notch','on','colors',[0 0 0],'symbol','','labels',{'RFP (cells TFP < 50)','RFP (cells TFP >= 50)',...
    'YFP (cells TFP < 50)','YFP (cells TFP >= 50)'});

suptitle({'dddW pathway: fluorescence value of all cells YFP-positive (thresh = 100)',...
        'each dot = fluo value of one cell (red = RFP; blue = TFP); overlaid points = mean with std  (11/16/2017 expt)'});

suptitle({'TFP-positive cells, average RFP and YFP (TFP thresh = 50), t13-32, to tease apart the 2 populations that emerge in dddW response',expt_date})
ylabel('average fluorescence intensity (A.U.)');
get(gca, 'XTick');
set(gca, 'FontSize', 16)


box_data = {};
box_data{1} = [y_tneg{t_range,q}];
box_data{2} = [yfp_fintcells_q_ypos{t_range,q}];
figure;
boxplot([y_tneg{t_range,q}],'test')
hold on;
boxplot([yfp_fintcells_q_ypos{t_range,q}],'2')

[h,p] = ttest2(x,y)


% histograms
figure
histogram(x,100,'Normalization','probability');
hold on;
histogram(y,100,'Normalization','probability');
legend({'RFP (cells TFP < 50)','RFP (cells TFP > 50)'},'FontSize',16)
xlabel('fluorescence intensity'); ylabel('normalized counts');
suptitle({'fluorescence intensity distribution of RFP in low or high TFP cell populations',expt_date});
get(gca, 'XTick');
set(gca, 'FontSize', 16)

figure
histogram(x2,100,'Normalization','probability');
hold on;
histogram(y2,100,'Normalization','probability');
legend({'YFP (cells TFP < 50)','YFP (cells TFP > 50)'},'FontSize',16)
xlabel('fluorescence intensity'); ylabel('normalized counts');
suptitle({'fluorescence intensity distribution of YFP in low or high TFP cell populations',expt_date});
get(gca, 'XTick');
set(gca, 'FontSize', 16)

%% yfp values 100-500 -> divergent populations?
q = 2;
for t = t_range;

    y_tneg_index{t,q} = find(y_tneg{t,q} > 100 & y_tneg{t,q} < 500);
    y_tpos_index{t,q} = find(y_tpos{t,q} > 100 & y_tpos{t,q} < 500);

    y_tneg_values{t,q} = y_tneg{t,q}(y_tneg_index{t,q});
    y_tpos_values{t,q} = y_tpos{t,q}(y_tpos_index{t,q});

    r_tneg_values{t,q} = r_tneg{t,q}(y_tneg_index{t,q});
    r_tpos_values{t,q} = r_tpos{t,q}(y_tpos_index{t,q});
    
    t_tneg_values{t,q} = t_tneg{t,q}(y_tneg_index{t,q});
    t_tpos_values{t,q} = t_tpos{t,q}(y_tpos_index{t,q});
end

figure;
for t = 25:30;%t_range;
scatter(t_tneg_values{t,q},y_tneg_values{t,q},'bo');hold on;
scatter(t_tpos_values{t,q},y_tpos_values{t,q},'ro')
end

xlabel('rfp')
ylabel('yfp')
suptitle({'red = TFP more than 50; blue = TFP less than 50; 100 > YFP > 500 cells only',expt_date})


%% make a plot with shaded errorbars
avgFluo_shaded = figure('units','normalized','outerposition',[0 0 1 1]);

q_range = q_with_yfp;%[1:total_cond];

plot_rfp = rfp_mean_q_ypos_ynorm_t;
plot_tfp = tfp_mean_q_ypos_ynorm_t;

error_rfp = rfp_sem_q_ypos_ynorm_t;
error_tfp = tfp_sem_q_ypos_ynorm_t;



% the real plot
for q = q_range;
    h(q) = subplot(1,total_cond,q); hold on;
    
    % specify data
    x = time_expt_hours(:,1);
    yr = plot_rfp(:,q);
    dyr = error_rfp(:,q);    
    yt = plot_tfp(:,q);
    dyt = error_tfp(:,q);
    

    % shaded errorbars
        % flipud -> can make continuous line that encapsulates a shape 
        % in this case, shape = -ve and +ve errorbars to be filled
    %fill([x;flipud(x)],[yr-dyr;flipud(yr+dyr)],[.9 .9 .9],'linestyle','none');
    % plot
    errorbar(x,yr,dyr,...
         marker_rfp,'LineWidth', linewidth_rfp, ...
         'Color', color_rfp_in, ...
         'MarkerEdgeColor', color_rfp_out, ...
         'MarkerFaceColor',color_rfp_in, ...
         'MarkerSize', markersize_rfp)
  
    
        
    errorbar(x,yt,dyt,...
         marker_tfp,'LineWidth', linewidth_tfp, ...
         'Color', color_tfp_in, ...
         'MarkerEdgeColor', color_tfp_out, ...
         'MarkerFaceColor',color_tfp_in, ...
         'MarkerSize', markersize_tfp)

    % fill([x;flipud(x)],[yt-dyt;flipud(yt+dyt)],[.9 .9 .9],'linestyle','none');

  
   % axis([0 10.02 0 18]);
    title(expt_conditions_string{q});
end

suptitle({'average fluorescence (for 3cv3 and 3cv4, YFP-positive cells only (thresh = 3), each cell YFP-normalized); for other strains all recognized cells included',...
    'each point = avg fluo taken of all cells in a given expt condition and time',...
          'errorbars = SEM  (12/14/2017 expt)'})
      
xlabel(h(5),'time from treatment (hours)','FontSize',16)
ylabel(h(1),'average fluorescence intensity, background subtracted (A.U.)','FontSize',16);


% Generate dummy info for plots with handles "m" outside of range of last subplot
m = zeros(3,1);
m(1) = plot(-1,-1, marker_rfp,'LineWidth', linewidth_rfp, ...
                 'Color', color_rfp_in, ...
                 'MarkerEdgeColor', color_rfp_out, ...
                 'MarkerFaceColor',color_rfp_in, ...
                 'MarkerSize', markersize_rfp); hold on;
m(2) = plot(-1,-1,marker_yfp,'LineWidth', linewidth_yfp, ...
                 'Color', color_yfp_in, ...
                 'MarkerEdgeColor', color_yfp_out, ...
                 'MarkerFaceColor',color_yfp_in, ...
                 'MarkerSize', markersize_yfp);
m(3) = plot(-1,-1,marker_tfp,'LineWidth', linewidth_tfp, ...
                 'Color', color_tfp_in, ...
                 'MarkerEdgeColor', color_tfp_out, ...
                 'MarkerFaceColor',color_tfp_in, ...
                 'MarkerSize', markersize_tfp+2); hold off;
% for marker legend
leg = {'RFP (dddW in 3cv3; dmdA in 3cv4)',...
       'YFP (constitutive)',...
       'TFP (dmdA in 3cv3; dddW in 3cv4)'}; % legend
legend(m(:),leg,'FontSize',16);

%axis(h(:),[0 15.4 0 10]);
hold off;

axis(h(4),[0 15.4 0 50]);

% save figure
% print(avgFluo_shaded,strcat(saveplace,'\avgFluo_YFPpos_YFPnorm_cellavg_axis2.png'),'-dpng','-r300'); % resolution 300

%clearvars ans h leg m q t xy avgFluo

axis(h([1, 2, 5, 6]),[0 15.4 0 0.8]);

%% plot t vs. fluorescence for each pathway (not multi-panel)
% [needs to be fixed]

plot_data = tfp_mean_q_ypos;
plot_error = tfp_sem_q_ypos;



figure('units','normalized','outerposition',[0 0 1 1]);
cmap = hsv(9);

normalization_factor = plot_data(total_time,8) - plot_data(1,8);
for q = 1:total_cond;
    
    plot_q = plot_data(:,q);
    plot_q = plot_q - plot_q(1); % bring down to 0
    plot_q = plot_q / normalization_factor;
 
        plot(time_expt_hours(:,1),plot_q,...
        'o-','markersize',15,'color',cmap(q,:),'linewidth',2); hold on;
    
  %  errorbar(time_expt_hours(:,1),plot_data(:,q),plot_error(:,q),...
  %      'o-','markersize',15,'color',cmap(q,:),'linewidth',2); hold on;
  clearvars plot_q
end
suptitle({'TFP (dmdA) fluorescence, YFPthresh = 50','subtract initial value and divide by final timepoint value of 1 mM DMSP (q = 8)'...
    ,expt_date})

legend(expt_conditions_string,'FontSize',16,'location','northwest')
xlim([0 15.5])
ylabel('normalized average fluorescence')
xlabel('time from treatment (hours)')
set(gca,'FontSize',20)

