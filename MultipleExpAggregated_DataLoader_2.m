%% aggregate multiple experiments: load data into structure array
% Version November 2018 
% UPDATED: 1/1/2019 to accommodate after mask_shift_fix
% FINALIZED as of 1/1/2019

% Load data from multiple experiments and store in data structures.

% In this version, load data BEFORE spectral leakage correction. 
% Spectral leakage correction is performed after loading data using
% aggregated b1. Also in this script, YFP thresholding and YFP
% normalization are performed.


% SUMMARY:
% 01) define directories

% 02) load and aggregate data [[ 01_aggdata.mat ]]
    % aggregate datasets into structures (data BEFORE spectral leakage correction)

  % 01_ParsedFileNames.mat
    % .date = date of experiment
    % .time = time points (hours since treatment)
    % .Ythresh = YFP threshold used
    
  % 06_FluoCellComputation_analysisonly (after background subtraction, before spectral leakage correction)
    % .Rc_raw = rfp_fintcells; raw RFP fluorescence of all cells by t,xy (background subtracted)
    % .Tc_raw = tfp_fintcells; raw TFP fluorescence of all cells by t,xy (background subtracted)
    % .Yc_raw = yfp_fintcells; raw YFP fluorescence of all cells by t,xy (background subtracted)
    % .Rm_raw = rfp_fintview; mean raw RFP fluo by t,xy (background subtracted)
    % .Tm_raw = tfp_fintview; mean raw TFP fluo by t,xy (background subtracted)
    % .Ym_raw = yfp_fintview; mean raw YFP fluo by t,xy (background subtracted)
    % .Rstd_raw = rfp_std; std of mean raw RFP fluo of cells by t,xy (background subtracted)
    % .Tstd_raw = tfp_std; std of mean raw TFP fluo of cells by t,xy (background subtracted)
    % .Ystd_raw = yfp_std; std of mean raw YFP fluo of cells by t,xy (background subtracted)
    % .Rsem_raw = rfp_sem; sem of mean raw RFP fluo of cells by t,xy (background subtracted)
    % .Tsem_raw = tfp_sem; sem of mean raw TFP fluo of cells by t,xy (background subtracted)
    % .Ysem_raw = yfp_sem; sem of mean raw YFP fluo of cells by t,xy (background subtracted)  
    % .R_bckg = rfp_bckg; mean of background pixels by t,xy (used for bckg subtraction)
    % .T_bkcg = tfp_bckg; mean of background pixels by t,xy (used for bckg subtraction)
    % .Y_bckg = yfp_bckg; mean of background pixels by t,xy (used for bckg subtraction)
    % .R_bckg_std = rfp_bckg_std; std of mean of background pixels by t,xy
    % .T_bkcg_std = tfp_bckg_std; std of mean of background pixels by t,xy
    % .Y_bckg_std = yfp_bckg_std; std of mean of background pixels by t,xy

  % 07_CellNumber_preYFPthresh
    % .cell_num{expt} = cell_num; % number of cells 
   
   % 08_CellSize_preYFPthresh
    % .cell_size{expt} = cell_size;
  
% 03) [PART 1] perform spectral leakage correction & visualize magnitude of spectral leakage correction using aggregated b1
% 03) [PART 2] perform spectral leakage correction & visualize magnitude of spectral leakage correction
% 03) [PART 3] SAVE spectral leakage corrected variables [[ 09_aggdata_SpectralLeakageCorrected.mat ]]

    % load aggregated b1
    % loop through each experiment's c_raw or fintcells (check) and apply aggregated b1
    % save leak corrected values 
        % c_leakcorr = fp_leakcorr_fintcells
        % m_leakcorr = fp_leakcorr_fintview
        % std_leakcorr = fp_leakcorr_std
        % sem_leakcorr = fp_leakcorr_sem
        
% 04) [PART 1] YFP positive threshold (collect YFP intensities by q and t)
% 04) [PART 2] save collected YFP intensities [[ 10_aggdata_YFPallcells_afterLeakCorr.mat ]]
% 04) [PART 3] Visualize distribution of YFP intensities of all cells
% 04) [PART 4] Apply YFP threshold (YFPpositiveFunc.m)
% 04) [PART 5] store YFP thresholding in permanent container in t,q format
% 04) [PART 6] save YFP thresholded values [[ 11_aggdata_YFPon_thresh_50.mat ]]

% 05) [PART 1] YFP normalize
% 05) [PART 2] save YFP normalized values [[12_aggdata_YFPnorm.mat]]



    
    
%% 01) define directories

motherdir = '/Volumes/Latte 4TB/_DMSPAvailabilityHypothesis_Datasets'; % folder containing all experiment data
cd(motherdir)


%% 02) load and aggregate data
% NOTE: ADD cell_num and cell_size (for YFP thresholding)
% UPDATED: 1/1/2019, for after mask_shift_fix

% list & describe desired experiments [MANUAL ENTRY]
expt_list = {'2018-02-24 (3cv3)',... % 3cv3 (unevenly focused)
             '2018-02-25 (3cv3)',... % 3cv3 (less cells in glucose initially)
             '2018-05-09 (3cv3)',... % 3cv3 (Kan25 instead of Kan10)
             '2018-05-10 (3cv3)',... % 3cv3 (Kan25 instead of Kan10)
             '2018-08-10 (3cv3)',... % 3cv3 (true rep 1)
             '2018-08-14 (3cv3)',... % 3cv3 (true rep 2)
             '2018-08-16 (3cv3)',... % 3cv3 (true rep 3; 6 xy positions per condition)
             '2018-02-28 (3cv4)',... % 3cv4 
             '2018-03-01 (3cv4)',... % 3cv4
             '2018-03-03 (3cv4)'};   % 3cv4

% extract the dates from the descriptions (for reading folders)
expt_date_folder = {}; % initialize
for i = 1:length(expt_list)
    expt_date_folder{i} = expt_list{i}(1:10); % the first 10 characters are dates
end

% list variables of interest          
vars_load = {'expt_date','time_expt_hours','analyzed_xy','expt_conditions_xy','expt_conditions_string',...
             'rfp_fintcells','tfp_fintcells','yfp_fintcells','rfp_fintview','tfp_fintview','yfp_fintview',...
             'rfp_std','tfp_std','yfp_std','rfp_sem','tfp_sem','yfp_sem',...
             'rfp_bckg','tfp_bckg','yfp_bckg','rfp_bckg_std','tfp_bckg_std','yfp_bckg_std','cell_num','cell_size'};
         
expt_range_processed = 5:10; % define the experiments desired for processing; normally 1:length(expt_date_folder)

aggdata = {}; % initialize

for expt = expt_range_processed
    
    cd(motherdir)
    expt_folder = dir(strcat(expt_date_folder{expt},'*')); % a string
    data_folder = strcat(motherdir,'/',expt_folder.name,'/analysis/DMSPkinetics'); % a string of directory containing desired data files
    cd(data_folder)
    
    load('01_ParsedFileNames.mat',vars_load{:});
    load('07_CellNumber_preYFPthresh',vars_load{:});
    load('08_CellSize_preYFPthresh',vars_load{:});
    
   % comment out if NOT after shift_mask_fixed
       if expt == 5 % for 8/10/2018 experiment, folder structure is a little different
        cd('shifted_mask_fix_2018-12-21/final')
       elseif expt == 10 % for 3/3/2018 experiment, folder structure is a little different
        cd('shifted_mask_fix_2018-12-21')
       else
        cd('shifted_mask_fix')       
       end
   
   
    load('06_FluoCellComputation_all',vars_load{:});
    
   % 01_ParsedFileNames.mat   
    aggdata.date{expt} = expt_list{expt};
    aggdata.time{expt} = time_expt_hours(:,1);
    aggdata.xy_all{expt} = analyzed_xy;
    aggdata.xy_q{expt} = expt_conditions_xy;
    aggdata.string_q{expt} = expt_conditions_string;
    
   % Cell recognition step: 06_FluoCellComputation_analysisonly (after background subtraction)
    aggdata.Rc_raw{expt} = rfp_fintcells;
    aggdata.Tc_raw{expt} = tfp_fintcells; 
    aggdata.Yc_raw{expt} = yfp_fintcells; 
    aggdata.Rm_raw{expt} = rfp_fintview; 
    aggdata.Tm_raw{expt} = tfp_fintview; 
    aggdata.Ym_raw{expt} = yfp_fintview; 
    aggdata.Rstd_raw{expt} = rfp_std; 
    aggdata.Tstd_raw{expt} = tfp_std;
    aggdata.Ystd_raw{expt} = yfp_std; 
    aggdata.Rsem_raw{expt} = rfp_sem; 
    aggdata.Tsem_raw{expt} = tfp_sem; 
    aggdata.Ysem_raw{expt} = yfp_sem; 
    aggdata.R_bckg{expt} = rfp_bckg; 
    aggdata.T_bkcg{expt} = tfp_bckg;
    aggdata.Y_bckg{expt} = yfp_bckg; 
    aggdata.R_bckg_std{expt} = rfp_bckg_std; 
    aggdata.T_bkcg_std{expt} = tfp_bckg_std; 
    aggdata.Y_bckg_std{expt} = yfp_bckg_std; 
    
   % cell number: 07_CellNumber_preYFPthresh
    aggdata.cell_num{expt} = cell_num;
   
   % cell size: 08_CellSize_preYFPthresh
    aggdata.cell_size{expt} = cell_size;
  
    clearvars -except i motherdir vars_load expt_list expt_date_folder aggdata expt_range_processed
   
end

clearvars ans i 

% MANUAL: need to fix 2/28/2018 expt_conditions_string (mislabeled as 3cv3, when it's actually 3cv4)
aggdata.string_q{8} = aggdata.string_q{9}


%% 03) [PART 1] perform spectral leakage correction & visualize magnitude of spectral leakage correction
% This step is time-consuming because need to correct spectral leakage for
% each cell.

% specify where figures will be saved (of magnitude of spectral leakage correction)
% cd(motherdir); saveplace = 'aggregated_analysis/spectral_leakage_correction'; mkdir(saveplace); 
saveplace_master = uigetdir % manually choose
saveplace = strcat(saveplace_master,'/spectral_leakage_correction'); mkdir(saveplace)

% save workspace
save(strcat(saveplace_master,'/01_aggdata.mat'))

% load aggregated b1
load('/Volumes/Latte 4TB/_DMSPAvailabilityHypothesis_Datasets/aggregated_b1/all expts_b1 final/RegressionSlope_b1_aggregated.mat')


%% 03) [PART 2] perform spectral leakage correction & visualize magnitude of spectral leakage correction

% SPECTRAL LEAKAGE CORRECTION
    % measured signal = (matrix of slopes) x (actual signal)
    % Y = b1*X  --> solve for X
    % X = Y*b1 or X = b1\Y

% specify handles
SaveFig = 1; % visualize the magnitude of spectral leakage correction for each color

h = waitbar(0,'out of total experiments'); % waitbar    

for expt = expt_range_processed % run for each experiment
    
 % WAITBAR
    waitbar(expt / length(expt_list));
 
 % INITIALIZE variables & temporary containers
    tfp_fintcells = []; yfp_fintcells = []; rfp_fintcells = [];  
    allc_fintcells = {}; leakcorrect_fintcells = {};
    tfp_leakcorr_fintcells = {}; yfp_leakcorr_fintcells = {}; rfp_leakcorr_fintcells = {};
    tfp_leakcorrdiff_fintcells = {}; yfp_leakcorrdiff_fintcells = {}; rfp_leakcorrdiff_fintcells = {};
    rfp_leakcorr_fintview = []; yfp_leakcorr_fintview = []; tfp_leakcorr_fintview = [];
    rfp_leakcorr_std = []; yfp_leakcorr_std = []; tfp_leakcorr_std = [];
    rfp_leakcorr_sem = []; yfp_leakcorr_sem = []; tfp_leakcorr_sem = [];
    
 % DEFINE variables for each experiment
    tfp_fintcells = aggdata.Tc_raw{expt};
    yfp_fintcells = aggdata.Yc_raw{expt};
    rfp_fintcells = aggdata.Rc_raw{expt};
    t_range = 1:length(aggdata.time{expt});
    xy_range = aggdata.xy_all{expt};
    q_range = 1:length(aggdata.xy_q{expt});
    
 % SPECTRAL LEAKAGE CORRECTION for each cell
    for t = t_range
        for xy = xy_range
            c = length(tfp_fintcells{t,xy});
            
            % if there is no image corresponding to the t,xy (e.g. due to evaporation), skip to next loop and assign NaN
                if c == 0                  
                   rfp_leakcorr_fintcells{t,xy} = [NaN]; tfp_leakcorr_fintcells{t,xy} = [NaN]; yfp_leakcorr_fintcells{t,xy} = [NaN];
                   rfp_leakcorr_fintview(t,xy) = [NaN]; tfp_leakcorr_fintview(t,xy) = [NaN]; yfp_leakcorr_fintview(t,xy) = [NaN];
                   rfp_leakcorr_std(t,xy) = [NaN]; tfp_leakcorr_std(t,xy) = [NaN]; yfp_leakcorr_std(t,xy) = [NaN];
                   rfp_leakcorr_sem(t,xy) = [NaN]; tfp_leakcorr_sem(t,xy) = [NaN]; yfp_leakcorr_sem(t,xy) = [NaN]; 
                    continue            
                end
            
            for cell = 1:c

                % 3x1 vector containing the fluo value of a cell in each color
                allc_fintcells{t,xy}{1,cell}(1,1) = tfp_fintcells{t,xy}(1,cell);
                allc_fintcells{t,xy}{1,cell}(2,1) = yfp_fintcells{t,xy}(1,cell);
                allc_fintcells{t,xy}{1,cell}(3,1) = rfp_fintcells{t,xy}(1,cell);

                % solve the matrix (SPECTRAL LEAKAGE CORRECTION)
                leakcorrect_fintcells{t,xy}{1,cell} = b1_agg \ allc_fintcells{t,xy}{1,cell}; 

                % separate
                tfp_leakcorr_fintcells{t,xy}(1,cell) = leakcorrect_fintcells{t,xy}{1,cell}(1,1);
                yfp_leakcorr_fintcells{t,xy}(1,cell) = leakcorrect_fintcells{t,xy}{1,cell}(2,1);
                rfp_leakcorr_fintcells{t,xy}(1,cell) = leakcorrect_fintcells{t,xy}{1,cell}(3,1);

                % difference between pre-correction and post-correction
                tfp_leakcorrdiff_fintcells{t,xy}{1,cell} = leakcorrect_fintcells{t,xy}{1,cell}(1,1) - allc_fintcells{t,xy}{1,cell}(1,1);
                yfp_leakcorrdiff_fintcells{t,xy}{1,cell} = leakcorrect_fintcells{t,xy}{1,cell}(2,1) - allc_fintcells{t,xy}{1,cell}(2,1);
                rfp_leakcorrdiff_fintcells{t,xy}{1,cell} = leakcorrect_fintcells{t,xy}{1,cell}(3,1) - allc_fintcells{t,xy}{1,cell}(3,1);

            end
            
            % calculate mean, std, sem values across cells for each t,xy and store temporarily
            tfp_leakcorr_fintview(t,xy) = mean(tfp_leakcorr_fintcells{t,xy});
            yfp_leakcorr_fintview(t,xy) = mean(yfp_leakcorr_fintcells{t,xy});
            rfp_leakcorr_fintview(t,xy) = mean(rfp_leakcorr_fintcells{t,xy});
            
            tfp_leakcorr_std(t,xy) = std(tfp_leakcorr_fintcells{t,xy});
            yfp_leakcorr_std(t,xy) = std(yfp_leakcorr_fintcells{t,xy});
            rfp_leakcorr_std(t,xy) = std(rfp_leakcorr_fintcells{t,xy});
            
            tfp_leakcorr_sem(t,xy) = std(tfp_leakcorr_fintcells{t,xy}) / sqrt(length(tfp_leakcorr_fintcells{t,xy}));
            yfp_leakcorr_sem(t,xy) = std(yfp_leakcorr_fintcells{t,xy}) / sqrt(length(yfp_leakcorr_fintcells{t,xy}));
            rfp_leakcorr_sem(t,xy) = std(rfp_leakcorr_fintcells{t,xy}) / sqrt(length(rfp_leakcorr_fintcells{t,xy}));            
                
        end % of t
    end % of xy
    
    
 % VISUALIZE leakage correction magnitude for each color for each experiment
    if SaveFig == 1
        
       % i) COLLECT leakage difference by q
            tfp_leakcorrdiff_fintcells_q = {}; yfp_leakcorrdiff_fintcells_q = {}; rfp_leakcorrdiff_fintcells_q = {}; % initialize
            for q_plot = q_range
                tfp_temp = []; yfp_temp = []; rfp_temp = [];

                for xy_plot = aggdata.xy_q{expt}{q_plot}
                    for t_plot = t_range
                        
                    t_temp = cell2mat(tfp_leakcorrdiff_fintcells{t_plot,xy_plot});   
                    tfp_temp = [tfp_temp,t_temp];

                    y_temp = cell2mat(yfp_leakcorrdiff_fintcells{t_plot,xy_plot});   
                    yfp_temp = [yfp_temp,y_temp];

                    r_temp = cell2mat(rfp_leakcorrdiff_fintcells{t_plot,xy_plot});   
                    rfp_temp = [rfp_temp,r_temp];
                    
                    end
                end

                % store in containers
                tfp_leakcorrdiff_fintcells_q{q_plot} = tfp_temp; 
                yfp_leakcorrdiff_fintcells_q{q_plot} = yfp_temp;
                rfp_leakcorrdiff_fintcells_q{q_plot} = rfp_temp;

                clearvars t_temp y_temp r_temp tfp_temp yfp_temp rfp_temp 

            end % of q_plot

       % ii) PLOT leakage correction magnitude histogram
            cmap = hsv(length(q_range));
            for FPcolors = 1:3
                % specify color
                if FPcolors == 1
                    data = rfp_leakcorrdiff_fintcells_q;
                    color = 'RFP';
                elseif FPcolors == 2
                    data = yfp_leakcorrdiff_fintcells_q;
                    color = 'YFP';
                elseif FPcolors == 3
                    data = tfp_leakcorrdiff_fintcells_q;
                    color = 'TFP';
                end

                % set figure to invisible
                fig_hist = figure('units','normalized','outerposition',[0 0 1 1],'visible','off'); 
                
                % cleanup previous histogram iteration variables
                clearvars histdata_q line_q q_hist 

                % plot for each conditions
                for q_hist = q_range
                    histdata_q(q_hist) = histogram(data{q_hist},'FaceColor',cmap(q_hist,:)); hold on;
                        centers = histdata_q(q_hist).BinEdges + histdata_q(q_hist).BinWidth/2; % center of each bar in histogram
                        heights = [histdata_q(q_hist).Values,0];        
                            hold on
                        line_q(q_hist) = plot(centers,heights,'LineWidth',3,'Color',cmap(q_hist,:));    
                end

                % figure settings
                title(strcat('magnitude of spectral leakage correction (',color,'); final value - original'),'FontSize',20)
                    xlabel('magnitude of spectral leakage correction (final value minus original)','FontSize',20);
                    ylabel('frequency','FontSize',20)
                set(histdata_q(:),'facealpha',0,'edgecolor','none'); % turn OFF bars
                legend(line_q(:),aggdata.string_q{expt},'FontSize',16,'Location','NorthWest'); 

                % save figure
                print(fig_hist,strcat(saveplace,'/LeakCorrMagnitude','_expt_',aggdata.date{expt},'_fp_',color),'-dpng'); close all;

            end % of FPcolors

    end % of SaveFig == 1;    
    
    
 % STORE calculations in permanent data container
    aggdata.Rc_leakcorr{expt} = rfp_leakcorr_fintcells; 
    aggdata.Tc_leakcorr{expt} = tfp_leakcorr_fintcells; 
    aggdata.Yc_leakcorr{expt} = yfp_leakcorr_fintcells;
    aggdata.Rm_leakcorr{expt} = rfp_leakcorr_fintview;
    aggdata.Tm_leakcorr{expt} = tfp_leakcorr_fintview;
    aggdata.Ym_leakcorr{expt} = yfp_leakcorr_fintview;
    aggdata.Rstd_leakcorr{expt} = rfp_leakcorr_std;
    aggdata.Tstd_leakcorr{expt} = tfp_leakcorr_std;
    aggdata.Ystd_leakcorr{expt} = yfp_leakcorr_std;
    aggdata.Rsem_leakcorr{expt} = rfp_leakcorr_sem;
    aggdata.Tsem_leakcorr{expt} = tfp_leakcorr_sem;     
    aggdata.Ysem_leakcorr{expt} = yfp_leakcorr_sem; 
    
end % of experiment

close(h)


%% 03) [PART 3] SAVE spectral leakage corrected variables 

% clean up after spectral leakage loops
clearvars -except aggdata b1_agg expt_date_folder expt_list vars_load motherdir expt_range_processed saveplace_master

% cd(motherdir); saveplace = 'aggregated_analysis'; % where figures will be saved (of magnitude of spectral leakage correction)
save(strcat(saveplace_master,'/09_aggdata_SpectralLeakageCorrected.mat'))


%% 04) YFP positive threshold
% [PART 1]  Collect YFP intensities by q and t

yfp_fintcells_leakcorr_q = {};
for expt = expt_range_processed
    
 % DEFINE variables for each experiment
    t_range = 1:length(aggdata.time{expt});
    q_range = 1:length(aggdata.xy_q{expt});
    xy_q_range = aggdata.xy_q{expt};
    inputdata = aggdata.Yc_leakcorr{expt};

    for q = q_range
        xy_range = xy_q_range{q};
        yfp_allcells_q_temp = []; yfp_allcells = [];

        for t = t_range
            for xy = xy_range
                yfp_allcells = [yfp_allcells inputdata{t,xy}];
                yfp_allcells_q_temp = [yfp_allcells_q_temp inputdata{t,xy}]; % deleted after every q
            end        
            yfp_fintcells_leakcorr_q{expt}{t,q} = yfp_allcells_q_temp; % store in permanent folder
        end  

    end % end of q

    yfp_fintcells_leakcorr{expt} = yfp_allcells; % for histogram plotting
            
end % end of expt
clearvars yfp_allcells_q_temp xy_range t q xy


%% 04) [PART 2] save collected YFP intensities

% cd(motherdir); saveplace = 'aggregated_analysis'; % where figures will be saved (of magnitude of spectral leakage correction)
save(strcat(saveplace_master,'/10_aggdata_YFPallcells_afterLeakCorr.mat'),'yfp_fintcells_leakcorr')


%% 04) [PART 3] Visualize distribution of YFP intensities of all cells
% make a decision about YFP threshold 
% NOTE: need to finalize the part for saving the figures

% cd(motherdir); 
% saveplace = strcat(saveplace_master,'/YFP_distribution_histogram'); mkdir(saveplace);
close all
thresh_test = [10 30 50 100 150 200]; % plotted on histogram
SaveIm = 1; % 1 = save image and work space
SameAxis = 1; % 1 = consistent x-axis limits for all experiments

% define figure popup
if SaveIm == 1
    vis_on_off = 'off';
else vis_on_off = 'on';
end


for expt = expt_range_processed
    
    % histogram and visualize threshold
    yfp_cells = figure('units','normalized','outerposition',[0 0 1 1],'visible',vis_on_off);

    % plot YFP distribution of the experiment
    histogram(yfp_fintcells_leakcorr{expt},1000);
        ylim_temp = ylim; % y axis limit
        hold on

    % plot threshold
    for i = 1:length(thresh_test)
        thresh_plot = thresh_test(i);
        plot(thresh_plot,0:ylim_temp(2)/50:ylim_temp(2),'r.','MarkerSize',10); % dot-plot of threshold
    end

    % figure labels
    title({'distribution of YFP fluorescence intensities of all xy of all time AFTER SPECTRAL LEAKAGE CORRECTION',...
           strcat('YFP threshold =',num2str(thresh_test)),...
           strcat('Experiment date: ',aggdata.date{expt})},'FontSize',20);
    ylabel('frequency of cells','FontSize',20);
    xlabel('YFP fluorescence intensity after background subtraction','FontSize',20);
    set(gca, 'FontSize', 20);
    
    % set axis limits
    if SameAxis == 1
        xlim([-50 1500]);
    end
    
    % save figure
    if SaveIm == 1 && SameAxis == 0
        savename = strcat('\\YFPcells_distribution_threshes','_afterLeakCorr_',expt_date_folder{expt},'.png');
        print(yfp_cells,strcat(saveplace,savename),'-dpng','-r300'); % resolution 300
    elseif SaveIm == 1 && SameAxis == 1
        savename = strcat('\\YFPcells_distribution_threshes','_afterLeakCorr_',expt_date_folder{expt},'_sameaxis.png');
        print(yfp_cells,strcat(saveplace,savename),'-dpng','-r300'); % resolution 300
    end
end
clearvars i expt SameAxis SaveIm savename thresh_plot thresh_test vis_on_off xy_q_range yfp_cells ylim_temp


%% 04) [PART 4] Apply YFP threshold
close all;
thresh_YFP_range = [50]; % all the YFP thresholds to run

% cd(motherdir); 
% saveplace_master_yfpthresh = strcat(saveplace_master,'/YFPthreshold_50'); mkdir(saveplace_master_yfpthresh);

for expt = expt_range_processed
    
    % q, xy, and t to be analyzed
    q_range = 1:length(aggdata.xy_q{expt});
    xy_range = aggdata.xy_all{expt};
    xy_q_range = aggdata.xy_q{expt}; 
    t_range = 1:length(aggdata.time{expt});
        
    % cell data (to gain stats about the number and size of cells that get rejected)
    cell_num = aggdata.cell_num{expt};
    cell_size = aggdata.cell_size{expt};
    
    % for each expt, define images that are skipped (in [t,xy] array format)
    [t_zero,xy_zero] = find(aggdata.Rm_raw{expt} == 0); % find t,xy indices where the values are 0 (i.e. eliminated t,xy images)
    skip_txy = [t_zero,xy_zero];
    
    % make a new save folder for each experiment (because savename is the same for all, 11_YFPon_thresh_...mat
    saveplace = strcat(saveplace_master_yfpthresh,'/',expt_date_folder{expt},'_Ythresh_results'); mkdir(saveplace)
    
    for thresh_applied = 1:length(thresh_YFP_range) % loop is relevant if there are more than 1 thresholds to apply

        thresh_YFP = thresh_YFP_range(thresh_applied);

        % fluorescence values to threshold on
        rfp_cells_fluo = aggdata.Rc_leakcorr{expt};
        yfp_cells_fluo = aggdata.Yc_leakcorr{expt};
        tfp_cells_fluo = aggdata.Tc_leakcorr{expt};
        
        % run script to threshold YFP
        [savename] = YFPpositiveFunc(skip_txy, thresh_YFP, xy_range, t_range, q_range, ...
                                     rfp_cells_fluo, yfp_cells_fluo, tfp_cells_fluo,...
                                     cell_num, cell_size, xy_q_range, saveplace);

    end 
    
    % cleanup after each experiment
    clearvars saveplace q_range xy_range xy_q_range t_range cell_num cell_size t_zero xy_zero 
    clearvars skip_txy thresh_applied thresh_YFP rfp_cells_fluo yfp_cells_fluo tfp_cells_fluo
end              
                

%% 04) [PART 5] store YFP thresholding in permanent container in t,q format
% this part is fine until I test it on the real set of experiments

% cd(motherdir); folder_dir = 'aggregated_analysis/YFP threshold'; % where the folder containing Ythresh results are located
cd(saveplace_master_yfpthresh)
Ythresh_of_interest = 50;

% record YFP threshold in permanent data container
aggdata.Ythresh = Ythresh_of_interest;

% define variables of interest
Ythresh_vars_of_interest = {'rfp_fintcells_q_ypos','tfp_fintcells_q_ypos','yfp_fintcells_q_ypos',...
                            'rfp_mean_q_ypos','tfp_mean_q_ypos','yfp_mean_q_ypos',...
                            'rfp_std_q_ypos','tfp_std_q_ypos','yfp_std_q_ypos',...
                            'rfp_sem_q_ypos','tfp_sem_q_ypos','yfp_sem_q_ypos',...
                            'cell_num_ypos','cell_num_ypos_avg_q','cell_num_ypos_std_q',...
                            'cell_size_ypos','cell_size_ypos_avg','cell_size_ypos_std'};

for expt = expt_range_processed

    % navigate to expt date directory
    % cd(motherdir); cd(folder_dir);
    cd(saveplace_master_yfpthresh);
    date_name = dir(strcat(expt_date_folder{expt},'*'));
    cd(date_name.name)
    
    % load variables on interest from experiment
    load(strcat('11_YFPon_thresh_',num2str(Ythresh_of_interest)),Ythresh_vars_of_interest{:});
    
    % load data into permanent containers
    aggdata.Rc_Ypos{expt} = rfp_fintcells_q_ypos; % raw RFP fluorescence of all cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata.Tc_Ypos{expt} = tfp_fintcells_q_ypos; % raw TFP fluorescence of all cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata.Yc_Ypos{expt} = yfp_fintcells_q_ypos; % raw YFP fluorescence of all cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied) -> this way, can do a histogram across all days
    aggdata.Rm_Ypos{expt} = rfp_mean_q_ypos; % mean % raw RFP fluo by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata.Tm_Ypos{expt} = tfp_mean_q_ypos; % mean raw TFP fluo by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata.Ym_Ypos{expt} = yfp_mean_q_ypos; % mean raw YFP fluo by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata.Rstd_Ypos{expt} = rfp_std_q_ypos; % std of mean raw RFP fluo of cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata.Tstd_Ypos{expt} = tfp_std_q_ypos; % std of mean raw TFP fluo of cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata.Ystd_Ypos{expt} = yfp_std_q_ypos; % std of mean raw YFP fluo of cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata.Rsem_Ypos{expt} = rfp_sem_q_ypos; % sem of mean raw RFP fluo of cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata.Tsem_Ypos{expt} = tfp_sem_q_ypos; % sem of mean raw TFP fluo of cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata.Ysem_Ypos{expt} = yfp_sem_q_ypos; %  sem of mean raw YFP fluo of cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)  
    aggdata.cell_num_ypos{expt} = cell_num_ypos; % number of cells that passed YFP threshold
    aggdata.cell_num_ypos_avg_q{expt} = cell_num_ypos_avg_q; % average number of cells that passed YFP threshold, by q
    aggdata.cell_num_ypos_std_q{expt} = cell_num_ypos_std_q; % std of number of cells that passed YFP threshold, by q
    aggdata.cell_size_ypos{expt} = cell_size_ypos; % size of cells that passed YFP threshold (cells that didn't pass have NaN)
    aggdata.cell_size_ypos_avg_xy{expt} = cell_size_ypos_avg; % average size of cells that passed YFP threshold, by xy
    aggdata.cell_size_ypos_std_xy{expt} = cell_size_ypos_std; % std of size of cells that passed YFP threshold, by xy
      
    % clean up workspace
    clearvars(Ythresh_vars_of_interest{:})
    
end


%% 04) [PART 6] save YFP thresholded values

% cd(motherdir); saveplace = 'aggregated_analysis'; % where figures will be saved (of magnitude of spectral leakage correction)
cd(saveplace_master)
save(strcat(saveplace_master,'/11_aggdata_YFPon_thresh_50.mat'),'aggdata','Ythresh_vars_of_interest','Ythresh_of_interest')


%% 05) [PART 1] YFP normalize (each cell by its own YFP)
% goal: each cell is normalized by its own YFP, and averaged by q for each expt

% [PART 1] divide RFP and TFP value of each cell by its YFP value
close all

for expt = expt_range_processed
    
    % define t,q to be analyzed
    t_range = 1:length(aggdata.time{expt});
    q_range = 1:length(aggdata.xy_q{expt});
  
    % initialize semi-permanent containers for each expt
    rfp_fintcells_q_leakcorr_ypos_ynorm = {}; tfp_fintcells_q_leakcorr_ypos_ynorm = {};
    rfp_fintview_q_leakcorr_ypos_ynorm = []; tfp_fintview_q_leakcorr_ypos_ynorm = []; 
    rfp_std_q_leakcorr_ypos_ynorm = []; tfp_std_q_leakcorr_ypos_ynorm = [];
    rfp_sem_q_leakcorr_ypos_ynorm = []; tfp_sem_q_leakcorr_ypos_ynorm = [];
    
    for q = q_range
        for t = t_range
            
            % initialize before each image
            rfp_cells_temp = []; yfp_cells_temp = []; tfp_cells_temp = [];
            rfp_ynorm_temp = []; tfp_ynorm_temp = [];
            
            % define fluorescence values of cells (bckg-subtracted, spectral leakage-corrected, YFP thresholded)
            rfp_cells_temp = aggdata.Rc_Ypos{expt}{t,q};
            yfp_cells_temp = aggdata.Yc_Ypos{expt}{t,q};
            tfp_cells_temp = aggdata.Tc_Ypos{expt}{t,q};

                for cell = 1:length(yfp_cells_temp) % refer to each cell  

                    % normalize each cell's fluorescence by its own YFP value
                    if isnan(yfp_cells_temp(cell)) == 0 % if the cell passed YFP threshold (i.e. has non-NaN value), normalize by YFP
                       rfp_ynorm_temp(cell) = rfp_cells_temp(cell)./yfp_cells_temp(cell);
                       tfp_ynorm_temp(cell) = tfp_cells_temp(cell)./yfp_cells_temp(cell);
                    elseif isnan(yfp_cells_temp(cell)) == 1 % if the cell did not pass YFP threshold, set normalized value to NaN
                       rfp_ynorm_temp(cell) = NaN;
                       tfp_ynorm_temp(cell) = NaN;
                    end      

                end
        
            % semi-permanent containers (calculate using non-NaN values)
            rfp_fintcells_q_leakcorr_ypos_ynorm{t,q} = rfp_ynorm_temp;
            tfp_fintcells_q_leakcorr_ypos_ynorm{t,q} = tfp_ynorm_temp;
            rfp_fintview_q_leakcorr_ypos_ynorm(t,q) = nanmean(rfp_ynorm_temp);
            tfp_fintview_q_leakcorr_ypos_ynorm(t,q) = nanmean(tfp_ynorm_temp);        
            rfp_std_q_leakcorr_ypos_ynorm(t,q) = nanstd(rfp_ynorm_temp);
            tfp_std_q_leakcorr_ypos_ynorm(t,q) = nanstd(tfp_ynorm_temp);
            rfp_sem_q_leakcorr_ypos_ynorm(t,q) = nanstd(rfp_ynorm_temp) / sqrt(sum(~isnan(rfp_ynorm_temp)));
            tfp_sem_q_leakcorr_ypos_ynorm(t,q) = nanstd(tfp_ynorm_temp) / sqrt(sum(~isnan(tfp_ynorm_temp)));
        
        end % of t
    end % of q
    
    % store in permanent container (by q,t)      
    aggdata.Rc_Ynorm{expt} = rfp_fintcells_q_leakcorr_ypos_ynorm; % YFP-normalized RFP fluorescence of each cell by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata.Tc_Ynorm{expt} = tfp_fintcells_q_leakcorr_ypos_ynorm; % YFP-normalized TFP fluorescence of each cell by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata.Rm_Ynorm{expt} = rfp_fintview_q_leakcorr_ypos_ynorm; % mean YFP-normalized RFP fluorescence of each cell by t,q
    aggdata.Tm_Ynorm{expt} = tfp_fintview_q_leakcorr_ypos_ynorm; % mean YFP-normalized TFP fluorescence of each cell by t,q
    aggdata.Rstd_Ynorm{expt} = rfp_std_q_leakcorr_ypos_ynorm; % std of mean YFP-normalized RFP fluo of cells, by t,q
    aggdata.Tstd_Ynorm{expt} = tfp_std_q_leakcorr_ypos_ynorm; % std of mean YFP-normalized TFP fluo of cells, by t,q
    aggdata.Rsem_Ynorm{expt} = rfp_sem_q_leakcorr_ypos_ynorm; % sem of mean YFP-normalized RFP fluo of cells, by t,q
    aggdata.Tsem_Ynorm{expt} = tfp_sem_q_leakcorr_ypos_ynorm; % sem of mean YFP-normalized TFP fluo of cells, by t,q
        
end

clearvars yfp_cells_temp rfp_cells_temp tfp_cells_temp
clearvars t q expt cell rfp_ynorm_temp tfp_ynorm_temp t_range
clearvars rfp_fintcells_q_leakcorr_ypos_ynorm tfp_fintcells_q_leakcorr_ypos_ynorm 
clearvars rfp_fintview_q_leakcorr_ypos_ynorm tfp_fintview_q_leakcorr_ypos_ynorm 
clearvars rfp_std_q_leakcorr_ypos_ynorm tfp_std_q_leakcorr_ypos_ynorm
clearvars rfp_sem_q_leakcorr_ypos_ynorm tfp_sem_q_leakcorr_ypos_ynorm
    
    

%% 05) [PART 2] save YFP normalized values
% and final aggdata

% cd(motherdir); saveplace = 'aggregated_analysis'; % where figures will be saved (of magnitude of spectral leakage correction)
save(strcat(saveplace_master,'/12_aggdata_YFPnorm.mat'),'aggdata')


%% 06) [PART 1] YFP normalize (each cell by average YFP at that time point)
% GOAL:   each cell is normalized by average YFP at that time point, and averaged by q for each expt
% ADDED:  12/2/2019 (in response to Reviewer #1)

load('14_aggdata_YFPnorm_t2avg.mat')

% [PART 1] divide RFP and TFP value of each cell by its YFP value
close all

for expt = 5:10 % expt_range_processed
    
    % define t,q to be analyzed
    t_range = 1:length(aggdata.time{expt});
    q_range = 1:length(aggdata.xy_q{expt});
  
    % initialize semi-permanent containers for each expt
    rfp_fintcells_q_leakcorr_ypos_ynorm = {}; tfp_fintcells_q_leakcorr_ypos_ynorm = {};
    rfp_fintview_q_leakcorr_ypos_ynorm = []; tfp_fintview_q_leakcorr_ypos_ynorm = []; 
    rfp_std_q_leakcorr_ypos_ynorm = []; tfp_std_q_leakcorr_ypos_ynorm = [];
    rfp_sem_q_leakcorr_ypos_ynorm = []; tfp_sem_q_leakcorr_ypos_ynorm = [];
    
    for q = q_range
        for t = t_range
            
            % initialize before each image
            rfp_cells_temp = []; yfp_cells_temp = []; tfp_cells_temp = [];
            rfp_ynorm_temp = []; tfp_ynorm_temp = [];
            
            % define fluorescence values of cells (bckg-subtracted, spectral leakage-corrected, YFP thresholded)
            rfp_cells_temp = aggdata.Rc_Ypos{expt}{t,q};
            yfp_cells_temp = aggdata.Yc_Ypos{expt}{t,q};
                yfp_norm_factor{expt}(t,q) = nanmean(yfp_cells_temp); % average YFP fluorescence at that q and t
            tfp_cells_temp = aggdata.Tc_Ypos{expt}{t,q};

                for cell = 1:length(yfp_cells_temp) % refer to each cell  

                    % normalize each cell's fluorescence by its own YFP value
                    if isnan(yfp_cells_temp(cell)) == 0 % if the cell passed YFP threshold (i.e. has non-NaN value), normalize by YFP
                       rfp_ynorm_temp(cell) = rfp_cells_temp(cell)./yfp_norm_factor{expt}(t,q); % yfp_cells_temp(cell);
                       tfp_ynorm_temp(cell) = tfp_cells_temp(cell)./yfp_norm_factor{expt}(t,q); % yfp_cells_temp(cell);
                    elseif isnan(yfp_cells_temp(cell)) == 1 % if the cell did not pass YFP threshold, set normalized value to NaN
                       rfp_ynorm_temp(cell) = NaN;
                       tfp_ynorm_temp(cell) = NaN;
                    end      

                end
        
            % semi-permanent containers (calculate using non-NaN values)
            rfp_fintcells_q_leakcorr_ypos_ynorm{t,q} = rfp_ynorm_temp;
            tfp_fintcells_q_leakcorr_ypos_ynorm{t,q} = tfp_ynorm_temp;
            rfp_fintview_q_leakcorr_ypos_ynorm(t,q) = nanmean(rfp_ynorm_temp);
            tfp_fintview_q_leakcorr_ypos_ynorm(t,q) = nanmean(tfp_ynorm_temp);        
            rfp_std_q_leakcorr_ypos_ynorm(t,q) = nanstd(rfp_ynorm_temp);
            tfp_std_q_leakcorr_ypos_ynorm(t,q) = nanstd(tfp_ynorm_temp);
            rfp_sem_q_leakcorr_ypos_ynorm(t,q) = nanstd(rfp_ynorm_temp) / sqrt(sum(~isnan(rfp_ynorm_temp)));
            tfp_sem_q_leakcorr_ypos_ynorm(t,q) = nanstd(tfp_ynorm_temp) / sqrt(sum(~isnan(tfp_ynorm_temp)));
        
        end % of t
    end % of q
    
    % store in permanent container (by q,t)      
    aggdata.Rc_Ynorm_eacht{expt} = rfp_fintcells_q_leakcorr_ypos_ynorm; % YFP-normalized RFP fluorescence of each cell by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata.Tc_Ynorm_eacht{expt} = tfp_fintcells_q_leakcorr_ypos_ynorm; % YFP-normalized TFP fluorescence of each cell by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata.Rm_Ynorm_eacht{expt} = rfp_fintview_q_leakcorr_ypos_ynorm; % mean YFP-normalized RFP fluorescence of each cell by t,q
    aggdata.Tm_Ynorm_eacht{expt} = tfp_fintview_q_leakcorr_ypos_ynorm; % mean YFP-normalized TFP fluorescence of each cell by t,q
    aggdata.Rstd_Ynorm_eacht{expt} = rfp_std_q_leakcorr_ypos_ynorm; % std of mean YFP-normalized RFP fluo of cells, by t,q
    aggdata.Tstd_Ynorm_eacht{expt} = tfp_std_q_leakcorr_ypos_ynorm; % std of mean YFP-normalized TFP fluo of cells, by t,q
    aggdata.Rsem_Ynorm_eacht{expt} = rfp_sem_q_leakcorr_ypos_ynorm; % sem of mean YFP-normalized RFP fluo of cells, by t,q
    aggdata.Tsem_Ynorm_eacht{expt} = tfp_sem_q_leakcorr_ypos_ynorm; % sem of mean YFP-normalized TFP fluo of cells, by t,q
        
end

clearvars yfp_cells_temp rfp_cells_temp tfp_cells_temp
clearvars t q expt cell rfp_ynorm_temp tfp_ynorm_temp t_range
clearvars rfp_fintcells_q_leakcorr_ypos_ynorm tfp_fintcells_q_leakcorr_ypos_ynorm 
clearvars rfp_fintview_q_leakcorr_ypos_ynorm tfp_fintview_q_leakcorr_ypos_ynorm 
clearvars rfp_std_q_leakcorr_ypos_ynorm tfp_std_q_leakcorr_ypos_ynorm
clearvars rfp_sem_q_leakcorr_ypos_ynorm tfp_sem_q_leakcorr_ypos_ynorm

%% 06) [PART 2] SAVE YFP normalize (each cell by average YFP at that time point)
save('15_aggdata_YFPnorm_eacht_avg.mat','aggdata','yfp_norm_factor')


%% below: ARCHIVE


%% calculate RT converted RFP after Ypos Ynorm
RT_conv_mean = mean([aggdata.RT_convfactor(:)])
RT_conv_std = std([aggdata.RT_convfactor(:)])

RT_converted = [aggdata.RT_convfactor(:) ./ RT_conv_mean]
RT_converted_mean = mean(RT_converted)
RT_converted_std = std(RT_converted)

aggdata.Rc_RTconv = {};

for i = 1:6;
    rcells_conv = {};
    for t = 1:32;
        for q = 1:9;
            rcells = [];
            rcells = aggdata.Rc_Ynorm{i}{t,q};
            rcells_conv{t,q} = [rcells ./ RT_conv_mean];
        end
    end
    aggdata.Rc_RTconv{i} = rcells_conv;
end


%% panels of scatter plots for each step of analysis (Vicente's idea)
% subplot with 4 rows and 3 columns
    % row 1 = raw data
    % row 2 = spectral leakage corrected
    % row 3 = YFP thresholding
    % row 4 = RT conversion
% only do this for glucose condition

% saveplace = uigetdir; % directory to save figures
i_range = 1:6; % experiment replicate to plot
q = 1;
xy_range = 1:7;
t = 5;
num_panel = 9; % number of panels
color = 'r';
markersize = 5;
SaveFig = 1; % save figure if 1
sameax = 1; % 1 = use same axis limits across experiments

x_data_source = {aggdata.Tc_raw,aggdata.Yc_raw,aggdata.Yc_raw,...
                 aggdata.Tc_leakcorr,aggdata.Yc_leakcorr,aggdata.Yc_leakcorr,...
                 aggdata.Tc_Ypos,aggdata.Tc_Ynorm,aggdata.Tc_Ynorm};
             
y_data_source = {aggdata.Rc_raw,aggdata.Tc_raw,aggdata.Rc_raw,...
                 aggdata.Rc_leakcorr,aggdata.Tc_leakcorr,aggdata.Rc_leakcorr,...
                 aggdata.Rc_Ypos,aggdata.Rc_Ynorm,aggdata.Rc_RTconv};
   
x_label = {'TFP raw','YFP raw','YFP raw',...
           'TFP leak corrected','YFP leak corrected','YFP leak corrected',...
           'TFP Ythresholded','TFP Ynorm','TFP Ynorm'};
       
y_label = {'RFP raw','TFP raw','RFP raw',...
           'RFP leak corrected','TFP leak corrected','RFP leak corrected',...
           'RFP Ythresholded','RFP Ynorm','RFP Ynorm & mean RT converted'};
       
axis_lim = { [-10 120 -10 20],...
    [-10 900 -10 120],...
    [-10 900 -10 20],...
    [-10 120 -10 20],...
    [-10 900 -10 120],...
    [-10 900 -10 20],...
    [-10 20 -10 120],...
    [-0.15 0.15 -0.15 0.15],...
    [-0.15 0.15 -0.15 0.15]};
       
for i = i_range;
    fig{i} = figure('units','normalized','outerposition',[0 0 1 1]);
    
    if ismember(i,[1]) == 1;
       alpha = 0.3; % louder dots for 2/25 expt only, where there are less cells
    else 
       alpha = 0.1;
    end
            
    for panel = 1:num_panel;
        subplot(3,3,panel)

        if ismember(panel,[1:6])==1;
            range = xy_range;
        else
            range = q;
        end

        % specify data source
        x_data = x_data_source{panel};
            x_data_plot = [x_data{i}{t,range}];

        y_data = y_data_source{panel};
            y_data_plot = [y_data{i}{t,range}];

        % plot
        scatter(x_data_plot,y_data_plot,markersize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',alpha);
        s(panel) = gca;
        
        xlabel(x_label{panel});
        ylabel(y_label{panel});
        
        if sameax == 1;
            axis(axis_lim{panel})     
        end
        
    end
    
    suptitle(strcat('Scatter panels of each data processing step (raw; spectral leakage correction; YFP thresholding; YFP normalization; RT conversion); t',...
        num2str(t),';',expt_date_list{i}))
        
    set(s(:),'FontSize',16)
    
    if SaveFig == 1;
       savename = strcat('\ScatterPanels_DataProcessingSteps_SameAxes_t',num2str(t),'_',expt_date_list{i})
       print(fig{i},strcat(saveplace,savename),'-dpng','-r300'); % resolution 300
       close all;
    end
    
    
end




%% general data loading
clearvars -except wrkDir
% wrkdir = uigetdir; cd(wrkdir); % directory where the .mat of different dates are stored
vars_load = {'expt_date',...
             'rfp_fintcells_q_leakcorr_ypos_ynorm',...
             'rfp_fintview_q_leakcorr_ypos_ynorm',...
             'rfp_std_q_leakcorr_ypos_ynorm',...
             'tfp_fintcells_q_leakcorr_ypos_ynorm',...
             'tfp_fintview_q_leakcorr_ypos_ynorm',...
             'tfp_std_q_leakcorr_ypos_ynorm',...   
             'rfp_sem_q_leakcorr_ypos_ynorm',...
             'tfp_sem_q_leakcorr_ypos_ynorm',...
             'time_expt_hours'};
expt_date_list_v3 = {'2018-02-25',... % 3cv3
                     '2018-03-02'};   % 3cv3
expt_date_list_v4 = {'2018-02-28',... % 3cv4
                     '2018-03-01',... % 3cv4
                     '2018-03-03'};   % 3cv4
Ythresh = 50;        
expt_date_list_temp = expt_date_list_v3; % MANUAL CHANGE TO INDICATE WHICH STRAIN to aggregate data for

for i = 1:length(expt_date_list_temp)
    
    % load
    mat2load_1 = strcat('12_YFP_normalized_Ythresh_',num2str(Ythresh),'_',expt_date_list_temp{i},'.mat')
    mat2load_2 = strcat('01_ParsedFileNames_',expt_date_list_temp{i},'.mat')
    load(mat2load_1,vars_load{:}); load(mat2load_2,vars_load{:});
    
  % 01_ParsedFileNames.mat   
    aggdata_temp.date{i} = expt_date;
    aggdata_temp.time{i} = time_expt_hours;
    aggdata_temp.Ythresh{i} = Ythresh; % manual entry above
    
  % 11_YFPon_thresh_50.mat    
    aggdata_temp.Rc_raw{i} = rfp_fintcells_q_ypos; % raw RFP fluorescence of all cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata_temp.Tc_raw{i} = tfp_fintcells_q_ypos; % raw TFP fluorescence of all cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata_temp.Yc_raw{i} = yfp_fintcells_q_ypos; % raw YFP fluorescence of all cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied) -> this way, can do a histogram across all days
    aggdata_temp.Rm_raw{i} = rfp_mean_q_ypos; % mean raw RFP fluo by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata_temp.Tm_raw{i} = tfp_mean_q_ypos; % mean raw TFP fluo by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata_temp.Ym_raw{i} = yfp_mean_q_ypos; % mean raw YFP fluo by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata_temp.Rstd_raw{i} = rfp_std_q_ypos; % std of mean raw RFP fluo of cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata_temp.Tstd_raw{i} = tfp_std_q_ypos; % std of mean raw TFP fluo of cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata_temp.Ystd_raw{i} = yfp_std_q_ypos; % std of mean raw YFP fluo of cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata_temp.Rsem_raw{i} = rfp_sem_q_ypos; % sem of mean raw RFP fluo of cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata_temp.Tsem_raw{i} = rfp_sem_q_ypos; % sem of mean raw TFP fluo of cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata_temp.Ysem_raw{i} = rfp_sem_q_ypos; % sem of mean raw YFP fluo of cells by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
            
  % from 13_YFP_normalized_Ythresh_50   
    aggdata_temp.Rc_Ynorm{i} = rfp_fintcells_q_leakcorr_ypos_ynorm; % YFP-normalized RFP fluorescence of each cell by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata_temp.Tc_Ynorm{i} = tfp_fintcells_q_leakcorr_ypos_ynorm; % YFP-normalized TFP fluorescence of each cell by t,q (background subtracted, spectral leakage corrected, YFP threshold applied)
    aggdata_temp.Rm_Ynorm{i} = rfp_fintview_q_leakcorr_ypos_ynorm; % mean YFP-normalized RFP fluorescence of each cell by t,q
    aggdata_temp.Tm_Ynorm{i} = tfp_fintview_q_leakcorr_ypos_ynorm; % mean YFP-normalized TFP fluorescence of each cell by t,q
    aggdata_temp.Rstd_Ynorm{i} = rfp_std_q_leakcorr_ypos_ynorm; % std of mean YFP-normalized RFP fluo of cells, by t,q
    aggdata_temp.Tstd_Ynorm{i} = tfp_std_q_leakcorr_ypos_ynorm; % std of mean YFP-normalized TFP fluo of cells, by t,q
    aggdata_temp.Rsem_Ynorm{i} = rfp_sem_q_leakcorr_ypos_ynorm; % sem of mean YFP-normalized RFP fluo of cells, by t,q
    aggdata_temp.Tsem_Ynorm{i} = tfp_sem_q_leakcorr_ypos_ynorm; % sem of mean YFP-normalized TFP fluo of cells, by t,q

  % RT_conversion.mat   
    aggdata_temp.RT_convfactor(i) = RT_conversion % 1 value per expt, R/T conversion factor from PA-FP signals

    clearvars -except i wrkDir vars_load Ythresh expt_date_list_temp aggdata_temp expt_date_list_v3 expt_date_list_v4
    
end

% assign to permanent variable
if isequal([expt_date_list_temp{:}],[expt_date_list_v3{:}]);
    aggdata_v3 = {};
    aggdata_v3 = aggdata_temp;
    save(strcat('aggdata_Ythresh_',num2str(Ythresh),'_v3.mat'),'aggdata_v3')
elseif isequal([expt_date_list_temp{:}],[expt_date_list_v4{:}]);
    aggdata_v4 = {}
    aggdata_v4 = aggdata_temp;
    save(strcat('aggdata_Ythresh_',num2str(Ythresh),'_v4.mat'),'aggdata_v4')
end

clearvars i ans aggdata_temp mat2load expt_date_list_temp

%% For FPctrl datasets
vars_load = {'yfp_fintcells_q_FPpos',...
             'rfp_fintview_q_leakcorr_ypos_ynorm',...
             'rfp_std_q_leakcorr_ypos_ynorm',...
             'tfp_fintcells_q_leakcorr_ypos_ynorm',...
             'tfp_fintview_q_leakcorr_ypos_ynorm',...
             'tfp_std_q_leakcorr_ypos_ynorm',...   
             'rfp_sem_q_leakcorr_ypos_ynorm',...
             'tfp_sem_q_leakcorr_ypos_ynorm',...
             'time_expt_hours'};
fp_fintcells_q_leakcorr_ypos_ynorm_s = {}; % structure for containing data from different expt dates

RegressionSlope_b1

clearvars -except wrkdir vars_load fp_fintcells_q_leakcorr_ypos_ynorm_s

load('10_FPon_thresh_200   30  200.mat',vars_load{:})
load('RegressionSlope_b1')
    fp_fintcells_q_leakcorr_ypos_ynorm_s.v4_Rc = rfp_fintcells_q_leakcorr_ypos_ynorm;
 
%% For box plots 

vars_load = {'rfp_leakcorr_fintcells',...
             'yfp_leakcorr_fintcells',...
             'tfp_leakcorr_fintcells',...
             'expt_conditions_xy',...
             'expt_date'};
         
expt_date_list_v3 = {'2018-02-25',... % 3cv3
                     '2018-03-02'};   % 3cv3
expt_date_list_v4 = {'2018-02-28',... % 3cv4
                     '2018-03-01',... % 3cv4
                     '2018-03-03'};   % 3cv4       
expt_date_list_temp = expt_date_list_v3; % MANUAL CHANGE TO INDICATE WHICH STRAIN to aggregate data for

for i = 1:length(expt_date_list_temp)
    
    % load
    mat2load_1 = strcat('09_SpectralLeakageCorrector_FPctrl_',expt_date_list_temp{i},'.mat')
    %mat2load_2 = strcat('01_ParsedFileNames_',expt_date_list_temp{i},'.mat')
    load(mat2load_1,vars_load{:});% load(mat2load_2,vars_load{:});
    xyrange_YFP = expt_conditions_xy{2};
    
    aggdata_PA_YFP_temp.xyrange_YFP{i} = xyrange_YFP;
    aggdata_PA_YFP_temp.date{i} = expt_date;
    aggdata_PA_YFP_temp.Rc{i} = [rfp_leakcorr_fintcells{xyrange_YFP}];
    aggdata_PA_YFP_temp.Yc{i} = [yfp_leakcorr_fintcells{xyrange_YFP}];
    aggdata_PA_YFP_temp.Tc{i} = [tfp_leakcorr_fintcells{xyrange_YFP}];    
    % normalize by YFP
    aggdata_PA_YFP_temp.Rc_Ynorm{i} = aggdata_PA_YFP_temp.Rc{i} ./ aggdata_PA_YFP_temp.Yc{i};
    aggdata_PA_YFP_temp.Tc_Ynorm{i} = aggdata_PA_YFP_temp.Tc{i} ./ aggdata_PA_YFP_temp.Yc{i};
       
    
    clearvars -except i wrkDir vars_load Ythresh expt_date_list_temp aggdata_PA_YFP_temp expt_date_list_v3 expt_date_list_v4
    
end

% assign to permanent variable
if isequal([expt_date_list_temp{:}],[expt_date_list_v3{:}]);
    aggdata_PAyfp_v3 = {};
    aggdata_PAyfp_v3 = aggdata_PA_YFP_temp;
    save(strcat('aggdata_PAyfp_v3.mat'),'aggdata_PAyfp_v3')
elseif isequal([expt_date_list_temp{:}],[expt_date_list_v4{:}]);
    aggdata_PAyfp_v4 = {}
    aggdata_PAyfp_v4 = aggdata_PA_YFP_temp;
    save(strcat('aggdata_PAyfp_v4'),'aggdata_PAyfp_v4')
end

clearvars i ans aggdata_temp mat2load expt_date_list_temp aggdata_PA_YFP_temp
