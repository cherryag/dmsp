%% Identify FP-positive cells in FPctrl cells
% version July 2018
% Run this code right after SpectralLeakageCorrector_FPctrl.m

% In this version, order of spectral leakage and FP thresholding was
% changed. With this change, runFPpositiveFunc_FPctrl.m comes AFTER  
% SpectralLeakageCorrector_FPctrl.m to be consistent with the DMSPkinetics analysis
% (where I also do spectral leakage correction before YFP thresholding).

% Summary of changes: 
% 1. swapped the order of analysis s.t. SpectralLeakageCorrector_FPctrl comes before runFPpositiveFunc_FPcontrol.
% 2. Calculate R/T conversion rate here, and feed to MultipleExpAggregated_DataLoader.m

% 00) load relevant data
% 01) collect FP intensities of all cells and save (10_FPallcells_leakcorr_collected.mat)
% 02) HISTOGRAMS (of all cells; colors separate and together)
% 03) apply FP threshold
% 04) calculate the RFP-TFP conversion factor

% % in the future consider refining code: 05) look at images of YFP [needs to be refined]


%% 00) load relevant data
clearvars -except analysis_fpctrl 
cd(analysis_fpctrl); % FPctrl analysis folder
load 09_SpectralLeakageCorrector_FPctrl.mat; % from FPctrl
saveplace = strcat(analysis_fpctrl,'/FP-positive figures'); mkdir(saveplace);


%% 1) Collect FP intensities
% feed spectral leakage-corrected data.
% separate colors, only taking into account relevant cells.

% collect all FP intensities of all cells
fp_allcells_leakcorr_colors = {}; % initialize final container
for q = 1:total_cond;
    fp_allcells_temp = []; % initialize
    
        % define colors
        if q == 1;
           fp_color = rfp_leakcorr_fintcells;
        elseif q == 2;
           fp_color = yfp_leakcorr_fintcells;
        elseif q == 3;
           fp_color = tfp_leakcorr_fintcells;
        end
    
    for xy = expt_conditions_xy{q};
         fp_allcells_temp = [fp_allcells_temp fp_color{xy}];
    end
    
    fp_allcells_leakcorr_colors{q} = fp_allcells_temp;
    
end
clearvars ans q xy fp_allcells_temp fp_color

save(strcat(analysis_fpctrl,'\10_FPallcells_leakcorr_collected.mat'),'fp_allcells_leakcorr_colors');


%% 2) HISTOGRAMS: Visualize distribution of FP intensities - colors separate
% Plot histograms and assess the appropriateness of FP threshold. 
close all;
SaveIm = 1; % save in the leaka
thresh_test = [200 30 200]; % usually [200 30 200]; plotted on histogram (rfp, yfp, tfp thresholds respectively)
nbin = 100;
color_bars = {'r','k','b'}; % in the order of q
color_string = {'RFP','YFP','TFP'};

for q = 1:3;    

    % first, plot single colors and visualize 
    fp_cells_separate = figure('units','normalized','outerposition',[0 0 1 1]);
    histogram(fp_allcells_leakcorr_colors{q},nbin,'EdgeColor',color_bars{q},'FaceColor','none');
        ylim_temp = ylim; % y axis limit
        hold on
    plot(thresh_test(q),0:ylim_temp(2)/50:ylim_temp(2),'.','Color',color_bars{q},'MarkerSize',10); % dot-plot of threshold

    % figure labels
    title({sprintf('distribution of fluorescence intensities (after spectral leakage correction); threshold = %0.f',thresh_test(q)),...
           strcat('Experiment date: ',expt_date, '   ;', expt_conditions_string{q})},'FontSize',20);
    ylabel('number of cells','FontSize',20);
    xlabel('fluorescence intensity after background subtraction','FontSize',20);

    % record the threshold
    thresh_allcolors(q) = thresh_test(q);

    if SaveIm == 1;
       savename = strcat(sprintf('\\FPallcells_distribution_thresh_%0.f',thresh_test(q)),'_',color_string{q},'.png');
       print(fp_cells_separate,strcat(saveplace,savename),'-dpng','-r300'); % resolution 300
    end
    
end


%% 2) HISTOGRAMS: Visualize distribution of FP intensities - all colors on the same plot
% Plot all histograms on the same figure. 

close all;
SaveIm = 1; % 1 = save figure
fp_cells = figure('units','normalized','outerposition',[0 0 1 1]);
for q = 1:total_cond;
    histogram(fp_allcells_leakcorr_colors{q},nbin,'EdgeColor',color_bars{q},'FaceColor','none');
        ylim_temp = ylim; % y axis limit
        hold on
    plot(thresh_allcolors(q),0:ylim_temp(2)/50:ylim_temp(2),'.','Color',color_bars{q},'MarkerSize',10); % dot-plot of threshold
end

% set x-axis limits
xlim([-200 10000]);

% figure labels
title({sprintf('distribution of fluorescence intensities (after spectral leakage correction); thresholds (RFP (red), YFP (black), TFP (blue)) = %0.f, %0.f, %0.f',thresh_allcolors(1),thresh_allcolors(2),thresh_allcolors(3)),...
       strcat('Experiment date: ',expt_date,'   ; # bins = ',num2str(nbin))},'FontSize',16);
ylabel('number of cells','FontSize',16);
xlabel('fluorescence intensity after background subtraction','FontSize',16);

% save figure
if SaveIm == 1;
    savename = strcat('\FPallcells_IntensityHistogram_nbin_',num2str(nbin),'.png');
    print(fp_cells,strcat(saveplace,savename),'-dpng','-r300'); % resolution 300
end

%% 2) clean
close all;
clearvars ans q ylim_temp ylim fp_cells nbin SaveIm savename thresh_test



%% 3) Apply FP thresholds

close all;
saveplace = analysis_fpctrl;

% q, xy, and t to be analyzed
q_range_analyzed = analyzed_q;
xy_range_analyzed = analyzed_xy;
time_range_analyzed = 1; % put 1 if no time component
skip_txy = [0,0]; % [0,0] if no image to skip; [0,xy] coordinates if there are images to skip in FPctrl

% data to feed into function
rfp_cell_data = rfp_leakcorr_fintcells;
yfp_cell_data = yfp_leakcorr_fintcells;
tfp_cell_data = tfp_leakcorr_fintcells;

% thresholds (in the order of RFP, YFP, TFP) 
thresh_allcolors = thresh_allcolors;

% run script to threshold FP 
    % in this script
        % q = 1 ==> RFP
        % q = 2 ==> YFP
        % q = 3 ==> TFP
[savename] = FPpositiveFunc_FPctrl(skip_txy, thresh_allcolors, q_range_analyzed, xy_range_analyzed, time_range_analyzed,...
                                   expt_conditions_xy, rfp_cell_data, yfp_cell_data, tfp_cell_data,...
                                   cell_num, cell_size, saveplace);
                        
                        
% clean up workplace before loading YFP-thresholded data
clearvars saveplace xy_range_analyzed t_range_analyzed q_range_analyzed skip_txy
clearvars yfp_fintcells_FPpos yfp_fintcells_FPneg tfp_fintcells_FPpos tfp_fintcells_FPneg rfp_fintcells_FPpos rfp_fintcells_FPneg;
clearvars yfp_fintview_FPneg yfp_fintview_FPpos tfp_fintview_FPneg tfp_fintview_FPpos rfp_fintview_FPeg rfp_fintview_FPpos;
clearvars cell_num_FPpos cell_num_FPneg percent_FPpos percent_FPneg
clearvars cell_num_FPpos_avg cell_num_FPpos_std cell_num_FPneg_avg cell_num_FPneg_std
clearvars cell_size_FPpos cell_size_FPneg cell_size_FPpos_avg cell_size_FPpos_std cell_size_FPneg_avg cell_size_FPneg_std
clearvars rfp_fintcells_q_FPpos tfp_fintcells_q_FPpos
clearvars rfp_mean_q_FPpos tfp_mean_q_FPpos rfp_std_q_FPpos tfp_std_q_FPpos rfp_sem_q_FPpos tfp_sem_q_FPpos

% load FP-thresholded data
    % in the form of 11_FPon_thresh_200   30  200.mat
load(savename);     


%% 4) calculate the RFP-TFP conversion factor
% Calculate the mean RFP or TFP signals of all RFP or TFP cells (after spectral leakage
% correction and after FP thresholding). Then calculate the RFP/TFP ratio using these mean values.

analysis_kinetics = uigetdir; % folder where DMSPkinetics is located, for saving the final variable
for colors = [1,3]; % RFP and TFP only
    xy_range = expt_conditions_xy{colors};
    
    if colors == 1; % RFP
        r_mean = nanmean([rfp_fintcells_FPpos{xy_range}])
    elseif colors == 3; % TFP
        t_mean = nanmean([tfp_fintcells_FPpos{xy_range}])
    end
          
end
RT_conversion = r_mean / t_mean
save(strcat(analysis_kinetics,'\RT_conversion.mat'),'RT_conversion')

%% NEXT STEP: MultipleExpAggregated_DataLoader.m
    % OR
% SpectralLeakageCorrector.m (for DMSPkinetics)
                  


