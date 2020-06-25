%% Identify FP-positive cells in FPctrl cells (AGGREGATED FPctrl data)
% version 1/1/2019

% Run this code right after SpectralLeakageCorrector_FPctrl_aggregatedb1.m
% In this script, calculate R/T conversion rate here, and feed to MultipleExpAggregated_DataLoader.m

% Note: runFPpositiveFunc_FPctrl.m comes AFTER SpectralLeakageCorrector_FPctrl.m 
%       to be consistent with the DMSPkinetics analysis (where I also do spectral 
%       leakage correction before YFP thresholding).

% SUMMARY: 
% 01) load relevant data
% 02) Collect FP intensities (separately for each expt)
% 03) HISTOGRAMS: Visualize distribution of FP intensities - colors separate
% 04) HISTOGRAMS: Visualize distribution of FP intensities - all colors on the same plot [[ NOT UPDATED for aggregated data format ]]
% 05) count number of cells in each expt at each t,xy
% 06) Apply FP thresholds (function FPpositiveFunc_FPctrl.m)
% 07) load each expt file and store relevant variables into fpagg
% 08) visualize FP thresholding results with histogram
% 09) calculate the RFP-TFP conversion factor
% 10) FINAL SAVE

% NEXT STEP: MultipleExpAggregated_DataLoader.m
    % OR
% SpectralLeakageCorrector.m (for DMSPkinetics)



%% 01) load relevant data

% load aggregated spectral leakage corrected FPctrl data 
    % data generated from 03_SpectralLeakageCorrector_FPctrl_aggregatedb1.m
load('/Volumes/Latte 4TB/_DMSPAvailabilityHypothesis_Datasets/aggregated_b1/all expts_b1 final/03_AggregatedData_SpectralLeakageCorrected.mat')

% change directory
cd('/Volumes/Latte 4TB/_DMSPAvailabilityHypothesis_Datasets/aggregated_b1/all expts_b1 final')


%% 02) Collect FP intensities (separately for each expt)
% collect all FP intensities of all cells (for making histogram in next step)
    % feed spectral leakage-corrected data
    % separated by colors

for expt = 1:length(expts_desired)
    
    % define data
    rfp_leakcorr_fintcells = fpagg(expt).rfp_leakcorr_fintcells;
    yfp_leakcorr_fintcells = fpagg(expt).yfp_leakcorr_fintcells;
    tfp_leakcorr_fintcells = fpagg(expt).tfp_leakcorr_fintcells;
    expt_conditions_xy = fpagg(expt).expt_conditions_xy;
    total_cond = 3; % number of FP colors
    
    % initialize
    fp_allcells_leakcorr_colors = {};
    
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
    
    % store in container
    fpagg_leakcorr_details(expt).fp_allcells_leakcorr_colors = fp_allcells_leakcorr_colors;
    
end % end of expts
clearvars ans q xy fp_allcells_temp fp_color


%% 03) HISTOGRAMS: Visualize distribution of FP intensities - colors separate
% Plot histograms and assess the appropriateness of FP threshold. 
% Separate histogram for each experiment.

% saveplace = strcat(pwd,'/FPpositive threshold histograms'); mkdir(saveplace)

close all;
SaveIm = 1; % save in the leaka
thresh_test = [200 30 200]; % usually [200 30 200]; plotted on histogram (rfp, yfp, tfp thresholds respectively)
nbin = 100;
color_bars = {'r','k','b'}; % in the order of q
color_string = {'RFP','YFP','TFP'};

for expt = 1:length(expts_desired)
    
    % define data (from previous section)
    fp_allcells_leakcorr_colors = fpagg_leakcorr_details(expt).fp_allcells_leakcorr_colors;
    expt_date = fpagg_leakcorr_details(expt).expt_date;
    expt_conditions_string = fpagg_leakcorr_details(expt).expt_conditions_string;
    
    for q = 1:3; % for each FP color    

        % make figure
        figure('units','normalized','outerposition',[0 0 1 1]);

        % plot single colors 
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
           savename = strcat(sprintf('/FPallcells_distribution_thresh_%0.f',thresh_test(q)),'_',color_string{q},'_expt_',num2str(expt),'.png');
           print(strcat(saveplace,savename),'-dpng','-r300'); % resolution 300
           close all
        end

    end % end of q
    
end % end of expt


%% 04) HISTOGRAMS: Visualize distribution of FP intensities - all colors on the same plot
% [[ NOT UPDATED ]]
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


%% 05) count number of cells in each expt at each t,xy
% needed for running FPpositive function

for expt = 1:length(expts_desired)
    
    rfp_fintcells = fpagg(expt).rfp_fintcells; % just to count total number of cells
    
    for t = 1
        for xy = 1:length(rfp_fintcells)
            
            fpagg(expt).cell_num(t,xy) = length(rfp_fintcells{t,xy});
            fpagg_leakcorr_details(expt).cell_num(t,xy) = length(rfp_fintcells{t,xy});           
            
        end % end of xy
    end % end of t

end


%% 06) Apply FP thresholds (function FPpositiveFunc_FPctrl.m)

% make a folder with all
saveplace_master = strcat(analysis,'/FPpositive threshold results'); mkdir(saveplace_master)
close all;

% define thresholds (in the order of RFP, YFP, TFP) 
thresh_allcolors = [200 30 200];

for expt = 1:length(expts_desired)
    
    % make a new folder for each expt to save data in (easier)
    saveplace = strcat(saveplace_master,'/expt_',num2str(expt),'_11_FPon result file'); mkdir(saveplace)

    % define data (from previous section)    
    rfp_leakcorr_fintcells = fpagg_leakcorr_details(expt).rfp_leakcorr_fintcells;
    yfp_leakcorr_fintcells = fpagg_leakcorr_details(expt).yfp_leakcorr_fintcells;
    tfp_leakcorr_fintcells = fpagg_leakcorr_details(expt).tfp_leakcorr_fintcells;
    expt_conditions_xy = fpagg_leakcorr_details(expt).expt_conditions_xy;
    
    % define q, xy, and t to be analyzed
    q_range_analyzed = 1:3; % for each FP color
    xy_range_analyzed = 1:length(rfp_leakcorr_fintcells); 
    time_range_analyzed = 1; % put 1 if no time component
    skip_txy = [0,0]; % [0,0] if no image to skip; [0,xy] coordinates if there are images to skip in FPctrl

    % data to feed into function
    rfp_cell_data = rfp_leakcorr_fintcells;
    yfp_cell_data = yfp_leakcorr_fintcells;
    tfp_cell_data = tfp_leakcorr_fintcells;

    cell_num =  fpagg_leakcorr_details(expt).cell_num;
    cell_size = fpagg_leakcorr_details(expt).rfp_fintcells; % THIS is a dummy variable

    % run script to threshold FP 
        % in this script
            % q = 1 ==> RFP
            % q = 2 ==> YFP
            % q = 3 ==> TFP
    [savename] = FPpositiveFunc_FPctrl(skip_txy, thresh_allcolors, q_range_analyzed, xy_range_analyzed, time_range_analyzed,...
                                       expt_conditions_xy, rfp_cell_data, yfp_cell_data, tfp_cell_data,...
                                       cell_num, cell_size, saveplace);

    savename_all{expt} = savename; % store savename for easy loading in the next step                             
                        
end % end of expt

% save savename_all
save(strcat(saveplace_master,'/FPon_results_directories.mat'),'savename_all')


%% 07) load each expt file and store relevant variables into fpagg

for expt = 1:length(expts_desired)

    clearvars -except expt analysis b1_agg expts_desired fpagg fpagg_leakcorr_details motherdir savename_all saveplace saveplace_master
    
    % load each experiment FPpos results file one by one
    load(savename_all{expt}, 'rfp_fintcells_FPpos','rfp_fintview_FPpos','rfp_sem_FPpos','rfp_std_FPpos',...
                             'yfp_fintcells_FPpos','yfp_fintview_FPpos','yfp_sem_FPpos','yfp_std_FPpos',...
                             'tfp_fintcells_FPpos','tfp_fintview_FPpos','tfp_sem_FPpos','tfp_std_FPpos',...
                             'cell_num_FPpos','thresh_allcolors')
                         
    % store in permanent container                   
    fpagg(expt).rfp_fintcells_FPpos = rfp_fintcells_FPpos;
    fpagg(expt).rfp_fintview_FPpos = rfp_fintview_FPpos;
    fpagg(expt).rfp_sem_FPpos = rfp_sem_FPpos;
    fpagg(expt).rfp_std_FPpos = rfp_std_FPpos;

    fpagg(expt).yfp_fintcells_FPpos = yfp_fintcells_FPpos;
    fpagg(expt).yfp_fintview_FPpos = yfp_fintview_FPpos;
    fpagg(expt).yfp_sem_FPpos = yfp_sem_FPpos;
    fpagg(expt).yfp_std_FPpos = yfp_std_FPpos;
    
    fpagg(expt).tfp_fintcells_FPpos = tfp_fintcells_FPpos;
    fpagg(expt).tfp_fintview_FPpos = tfp_fintview_FPpos;
    fpagg(expt).tfp_sem_FPpos = tfp_sem_FPpos;
    fpagg(expt).tfp_std_FPpos = tfp_std_FPpos;
    
    fpagg(expt).cell_num_FPpos = cell_num_FPpos;
    fpagg(expt).thresh_allcolors = thresh_allcolors;
    
end
clearvars -except analysis b1_agg expts_desired fpagg fpagg_leakcorr_details motherdir savename_all saveplace saveplace_master


%% 08) visualize FP thresholding results with histogram

expt = 10; % MANUAL just for test
nbin = 100;
rfp_fintcells_FPpos = [fpagg(expt).rfp_fintcells_FPpos{fpagg(expt).expt_conditions_xy{1}}];
yfp_fintcells_FPpos = [fpagg(expt).yfp_fintcells_FPpos{fpagg(expt).expt_conditions_xy{2}}];
tfp_fintcells_FPpos = [fpagg(expt).tfp_fintcells_FPpos{fpagg(expt).expt_conditions_xy{3}}];

figure('units','normalized','outerposition',[0 0 1 1]); 
histogram(tfp_fintcells_FPpos,nbin,'EdgeColor','b','FaceColor','none'); hold on
histogram(rfp_fintcells_FPpos,nbin,'EdgeColor','r','FaceColor','none');
histogram(yfp_fintcells_FPpos,nbin,'EdgeColor','k','FaceColor','none');


%% 09) calculate the RFP-TFP conversion factor (all cells pooled)

% Calculate the mean RFP or TFP signals of all RFP or TFP cells (after spectral leakage
% correction and after FP thresholding). Then calculate the RFP/TFP ratio using these mean values.

% pool all cell fluorescence (spectral leakage corrected and FP positive) of all experiments
r_cells_pooled = []; t_cells_pooled = [];
for expt = 1:length(expts_desired)
    all_xy = 1:length(fpagg(expt).rfp_fintcells); 
    
    for colors = [1,3]; % RFP and TFP only
        
        % define the xy positions corresponding to the color
        xy_range = all_xy(fpagg(expt).expt_conditions_xy{colors}); 
        
        % initialize
        data = [];
        
        if colors == 1; % RFP
            data = [fpagg(expt).rfp_fintcells_FPpos{xy_range}];
            r_cells_pooled = [r_cells_pooled data];
        elseif colors == 3; % TFP
            data = [fpagg(expt).tfp_fintcells_FPpos{xy_range}];
            t_cells_pooled = [t_cells_pooled data];
        end

    end % end of colors

end % end of expt

% calculate mean
t_mean = nanmean(t_cells_pooled)
r_mean = nanmean(r_cells_pooled)
t_mean_std = nanstd(t_cells_pooled)
r_mean_std = nanstd(r_cells_pooled)

% calculate RT conversion
RT_conversion_agg = r_mean / t_mean

% save
save(strcat(analysis,'/RTconversion_factor_aggregated_pooledcells.mat'),'RT_conversion_agg','r_cells_pooled','t_cells_pooled','r_mean','t_mean')


%% 09) ALTERNATIVE: calculate the RFP-TFP conversion factor (mean of individual experiments)

% Better method than previous section: First calculate mean fluo within each experiment, 
% then calculate mean across experiments, then take the ratio. 

% pool all cell fluorescence (spectral leakage corrected and FP positive) of all experiments
r_mean_expt = []; t_mean_expt = []; r_std_expt = []; t_std_expt = [];
for expt = 1:length(expts_desired)
    all_xy = 1:length(fpagg(expt).rfp_fintcells); 
    
    for colors = [1,3]; % RFP and TFP only
        
        % define the xy positions corresponding to the color
        xy_range = all_xy(fpagg(expt).expt_conditions_xy{colors}); 
        
        % initialize
        data = [];
        
        if colors == 1; % RFP
            data = [fpagg(expt).rfp_fintcells_FPpos{xy_range}]; % all cells in a single experiment
            r_mean_expt(expt) = nanmean(data); 
            r_std_expt(expt) = nanstd(data);
        elseif colors == 3; % TFP
            data = [fpagg(expt).tfp_fintcells_FPpos{xy_range}]; % all cells in a single experiment
            t_mean_expt(expt) = nanmean(data);
            t_std_expt(expt) = nanstd(data);
        end

    end % end of colors

end % end of expt

% calculate mean
t_mean = mean(t_mean_expt)
r_mean = mean(r_mean_expt)
t_std = std(t_mean_expt)
r_std = std(r_mean_expt)

% calculate RT conversion
RT_conversion_agg_exptmean = r_mean / t_mean

% save
save(strcat(analysis,'/RTconversion_factor_aggregated_exptmean.mat'),'RT_conversion_agg_exptmean','r_mean_expt','r_std_expt','t_mean_expt','t_std_expt','t_mean','r_mean',...
    't_std','r_std','expts_desired')


%% 10) FINAL SAVE
clearvars -except analysis b1_agg expts_desired fpagg fpagg_leakcorr_details motherdir savename_all saveplace_master RT_conversion_agg

% save the final analysis results
save(strcat(analysis,'/04_AggregatedData_FPon_RTconv_pooledcells.mat'))


%% NEXT STEP: MultipleExpAggregated_DataLoader.m
    % OR
% SpectralLeakageCorrector.m (for DMSPkinetics)
                 
