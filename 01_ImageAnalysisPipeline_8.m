%% Image Analysis Pipeline
% version May 2019 (for Mac)
% convert images to double before fluorescence calculation (doesn't change from before, but for formality)

% 01) define directories 
% 02) parse file names & expt conditions
% 03) recognize cells
% 04) create phase-mask images for visual inspection   [PARTS 1-2]
% 05) segment cells and calculate area and UQ
% 06) Cell Property Collector & Histogram Generator   [PARTS 1-3]
% 07) thresholds on area and UQ   
% 08) make thresholded mask figures for visual inspection
% 09) visualize the number of eliminated cells   [PARTS 1-2]
% 10) create background mask   [PARTS 1-4]
% 11) calculate fluorescence   [PARTS 1-2]
% 12) extract other information - number of cells
% 13) extract other information - cell size
% 14) final save
% NEXT STEP: SPECTRAL LEAKAGE CORRECTION using aggregated b1



%% 01) define directories 
% reverse slash direction for Mac

wrkDir = uigetdir(); cd(wrkDir); % folder with TIFF images
analysis = uigetdir(); % where analysis files will be saved
expt_date = '5/16/2017 (3cv3 sub-1uM DMSP kinetics)'; % for titles


%% 02) parse file names & expt conditions

imgs = dir('*.tif');
lastImg = imgs(end).name
filename_base = 'DMSPconc_3cv3old_kinetics_'; % base of tiff file names
parsedFileName = strsplit(lastImg,{filename_base,'t','xy','c','.tif'})
total_time = str2num(parsedFileName{2}) % total number of time points
total_xy = str2num(parsedFileName{3}) % total number of positions
total_channel = str2num(parsedFileName{4}) % total number of channels

% matrices of TIFF image file names with t rows x p columns 
for t = 1:total_time;
    for xy = 1:total_xy;
     bfImgs(t,xy) = dir(strcat(filename_base,'t',num2str(t,'%02i'),'xy',num2str(xy,'%02i'),'c1.tif'));
     RFP(t,xy) = dir(strcat(filename_base,'t',num2str(t,'%02i'),'xy',num2str(xy,'%02i'),'c2.tif'));
     YFP(t,xy) = dir(strcat(filename_base,'t',num2str(t,'%02i'),'xy',num2str(xy,'%02i'),'c3.tif'));
     TFP(t,xy) = dir(strcat(filename_base,'t',num2str(t,'%02i'),'xy',num2str(xy,'%02i'),'c4.tif'));
    end
end
clearvars imgs lastImg parsedFileName t xy ans

% parse through xy positions (expt conditions)
total_xy_per_cond = 7;
total_cond = total_xy / total_xy_per_cond
xy = [1:total_xy_per_cond:total_xy];
difference = total_xy_per_cond - 1;

for i = 1:total_cond
    expt_conditions_xy{i} = [xy(i):xy(i)+difference];  
end
clearvars ans xy difference i

% ALTERNATIVE: manual parse through xy positions (expt conditions)
    % used for experiments with different # of xy positions per condition
expt_conditions_xy = {[1:4] [5:13] [14:23] [24:31] [32:41] [42:51] [52:60] [61:69] [70:77]}; % MANUAL entry
total_cond = length(expt_conditions_xy) 
total_xy_per_cond = 10; % MANUAL entry, the largest xy positions

    % collapse expt_conditions_xy (if some q are the same)
    % expt_conditions_xy{6} = [expt_conditions_xy{6},expt_conditions_xy{9}]; % append xy positions
    % expt_conditions_xy{9} = []; % delete xy contents of q9

% analyzed xy - if only analyzing certain expt conditions (q) 
analyzed_q = 1:total_cond;
analyzed_xy = [];
for q = analyzed_q       
    analyzed_xy = [analyzed_xy expt_conditions_xy{q}];
end
clearvars q
    
% experimental conditions in words (MANUAL CHANGE)
% 5/16/2017 experiment
expt_conditions_string = {'3cv3 - succinate 1 mM',...
                          '3cv3 - DMSP 1 mM',...
                          '3cv3 - DMSP 500 \muM',...
                          '3cv3 - DMSP 100 \muM',...
                          '3cv3 - DMSP 50 \muM',...
                          '3cv3 - DMSP 25 \muM',...
                          '3cv3 - DMSP 1 \muM',...
                          '3cv3 - DMSP 500 nM',...
                          '3cv3 - DMSP 100 nM'}; % for legend
                         
                      
% define time range
start_timepoint = 24; % minutes after treatment (12:36 treat; 13:00 first acquisition)
time_interval = 30; % minutes between time points
final_timepoint = start_timepoint + (time_interval * (length(bfImgs(:,1)) - 1 ) ); % total minutes after treatment
time_expt = [start_timepoint : time_interval : final_timepoint] % calculate time points (should be same size as avgFluo)
time_expt_hours = time_expt ./ 60 % in hours
time_expt_hours = repmat(time_expt_hours,total_xy_per_cond,1)' % for plotting # xy positions per condition

clearvars start_timepoint time_interval final_timepoint time_expt

% save workplace
save(strcat(analysis,'/01_ParsedFileNames.mat')); clearvars saveplace;
save(strcat(analysis,'/time_expt.mat'),'time_expt_hours'); 


%% 03) recognize cells
% run script CellRecognizer.m
% this step takes a long time if saving mask overlay images

% cd(wrkDir);
% saveplace = strcat(analysis,'/masks_2'); mkdir(saveplace);

close all;
threshold = 820; % for testing phase
SaveFlag = 1; % 1 = save images without figures popping up; 0 = figures pop up
redo = 1; % 1 = doesn't clear bwf and bwf_5, for re-analyzing masks

% determine which images to run through for-loop
if SaveFlag == 0; % for tests, only look at representative images
    q_test = 3;
    time_range = [1 10];
    xy_range = expt_conditions_xy{q_test}(5);
elseif SaveFlag == 1 & redo == 0; % all images
    
    time_range = 1:total_time; % all time points    
        time_range1 = time_range;
        time_range2 = time_range;
        time_range3 = time_range;
        time_range4 = time_range;
        
    xy_range = analyzed_xy;
        xy_range1 = expt_conditions_xy{1};
        xy_range2 = [expt_conditions_xy{2:9}];
        xy_range3 = NaN;
        xy_range4 = NaN;
        
    threshold1 = 1000; 
    threshold2 = 800;
    threshold3 = NaN;
    threshold4 = NaN; 
   
    bwf_5 = {}; % initialize mask container, BEFORE border clearing
    bwf = {}; % initialize mask container, FINAL mask
    
    h = waitbar(0,'out of total time points'); steps = length(time_range); % waitbar    
elseif SaveFlag == 1 & redo == 1; % don't clear bwf_5 and bwf if redo
    
    time_range = 1:total_time; % all time points    
        time_range1 = time_range;
        time_range2 = time_range;
        time_range3 = time_range;
        time_range4 = time_range;
        
    xy_range = expt_conditions_xy{3}; % MANUAL CHANGE for xy positions that need to be re-analyzed
        xy_range1 = expt_conditions_xy{3};
        xy_range2 = NaN;
        xy_range3 = NaN;
        
    threshold1 = 820; 
    threshold2 = NaN;
    threshold3 = NaN;
    threshold4 = NaN; 
    h = waitbar(0,'out of total time points'); steps = length(time_range); % waitbar    
end

for t = time_range;
        if SaveFlag == 1;
            waitbar(t / length(time_range));
        end
    for xy = xy_range;
        
        % if there is no image corresponding to the t,xy
        if isempty(bfImgs(t,xy).name) == 1;
            bwf_5{t,xy} = [];
            bwf{t,xy} = [];
            continue            
        end
        
        % specify threshold according to time range and xy range
        if SaveFlag == 1;
            
           if ismember(xy,xy_range1) == 1;            
                if ismember(t,time_range1) == 1;
                    threshold = threshold1;
                elseif ismember(t,time_range2) == 1;
                    threshold = threshold1;
                end               
           elseif ismember(xy,xy_range2) == 1; 
                if ismember(t,time_range3) == 1;
                    threshold = threshold2;
                elseif ismember(t,time_range4) == 1;
                    threshold = threshold2;
                end
           end % end of ismember(xy)
           
        end % end of saveflag
                           
        phasename = bfImgs(t,xy).name;
        [bwf_5{t,xy} bwf{t,xy}] = CellRecognizer(SaveFlag,saveplace,phasename,threshold);

    end
end
if SaveFlag == 1;
    close(h)
end
clearvars xy_range time_range q ans
clearvars h t xy phasename steps 


% save workspace
if SaveFlag == 1;
    save(strcat(analysis,'/02_CellRecognizer_1.mat')); % includes bwf
    if redo == 1;
        save(strcat(analysis,'/02_CellRecognizer_2.mat')); % includes bwf
    end
    save(strcat(analysis,'/bwf_final.mat'),'bwf'); % save the final masks variable
    save(strcat(analysis,'/bwf_5_final.mat'),'bwf_5'); % save the final masks variable
    clearvars SaveFlag redo threshold
    clearvars time_range1 time_range2 threshold1 threshold2
    clearvars time_range3 time_range4 threshold3 threshold4
    clearvars xy_range1 xy_range2 xy_range3
end


%% 04) [PART 1] create phase-mask images (ADDED MARCH 2018) 
% visually check ensure mask formation success

close all;
SaveIm = 1; % 1 = save image; 0 = don't save image, set to 0 for redo if only a few xy positions
redo = 1; % make sure that SaveIm = 0 if redo = 1 and only a few xy positions

% each figure is a time-lapse of one xy position (max. 35 time points per image)
subplot_row = 4; 
subplot_col = 8; % normally 8 
total_im_per_subplot = subplot_row * subplot_col; 

% define xy range for the loop
if redo == 0;
    xy_range_loop = [1:total_xy total_xy+1]; % 1 larger than total_xy
elseif redo == 1;
    xy_range = expt_conditions_xy{3}; % MANUAL ENTRY for redo
    xy_range_loop = [xy_range];
end


fig_count = 0; % initialize fig_count (always start with 0)
if SaveIm == 1 & redo == 0; % no waitbar for redo
     h = waitbar(0,'out of total xy positions'); steps = total_xy; % waitbar    
end
for xy = xy_range_loop;
    for t = 1:total_time;
        if SaveIm == 1;
            waitbar(xy / total_xy);
        end    
        
    % SAVE FIGURE after the first image
    if (t == 1 & xy > xy_range_loop(1) & SaveIm == 1) | (xy == total_xy+1 & SaveIm == 1);
       % save figure (high resolution)
       print(f(fig_count),strcat(saveplace,'/','maskoverlay_',num2str(fig_count),'.png'),'-dpng','-r300');    
            if xy == total_xy+1;
                close(h)
                return
            end
    end
    
    % every multi-plot image, create a new figure
    if t == 1;
        i = 1; % re-initialize i for subplot
        fig_count = fig_count + 1; % advance figure count   
        f(fig_count) = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        if SaveIm == 1;
           set(f(fig_count), 'Visible', 'off'); 
        end  
    end
       
    % load phase image
    phasename = bfImgs(t,xy).name;
    ph_adj_temp = imread(phasename);
    ph_adj_temp_vis = imadjust(ph_adj_temp);
   
    % create 4 x 6 matrix of images
    subplot(subplot_row,subplot_col,i);
        imshowpair(ph_adj_temp_vis, bwf{t,xy}); hold on;
        title(phasename) 
    % tweak the positioning and size of each subplot
    original_pos = get(gca, 'Position'); % Position = [left bottom width height]
    %new_pos = original_pos + [-0.05 -0.05 0.06 0.06]; % good for 4x6
    new_pos = original_pos + [-0.04 -0.04 0.03 0.03]; % good for 4x8
    %new_pos = original_pos + [-0.04 -0.04 0.045 0.045]; % good for 5x7
    set(gca, 'Position', new_pos);
    
    % advance i for subplot position
    i = i + 1;
    
    end
end

% clean workspace
if SaveIm == 1 | redo == 1;
    clearvars redo h t xy xy_range xy_range_loop q_range ans thresh_divide q q1 q2 q3
    clearvars threshold threshold1 threshold2 threshold3 phasename  
    clearvars f fig_count i new_pos original_pos ph_adj_temp ph_adj_temp_vis
    clearvars saveplace SaveIm subplot_col subplot_row total_im_per_subplot xy_figure
end


%% 04) [PART 2] optional; enter manually images to exclude, if any (e.g. due to evaporation)
t_exclude = [1:total_time];
xy_exclude = [68];

for t = t_exclude;
    for xy = xy_exclude;
        bfImgs(t,xy).name = [];
        bwf{t,xy} = {}; % make mask container empty
        bwf_5{t,xy} = {}; % make mask container empty
    end
end
clearvars t_exclude xy_exclude t xy

% save
save(strcat(analysis,'/bwf_final.mat'),'bwf'); % save the final masks variable
save(strcat(analysis,'/bwf_5_final.mat'),'bwf_5'); % save the final masks variable


%% 05) segment cells and calculate area and UQ

% initialize containers
Q75 = {}; area = {};

h = waitbar(0,'out of total time points'); steps = total_time;
for t = 1:total_time;
    waitbar(t / total_time);
    for xy = analyzed_xy;  
               
        % if there is no image corresponding to the t,xy
        if isempty(bfImgs(t,xy).name) == 1;
            bwf_5{t,xy} = [];
            bwf{t,xy} = [];
            continue            
        end
        
    phasename = bfImgs(t,xy).name;
    bwf_txy = bwf{t,xy};

    [cells{t,xy} area{t,xy} Q75{t,xy}] = CellPropertyCalculator(bwf_txy, phasename);

    end
end
close(h)
clearvars t xy phasename bwf_txy h steps ans saveplace threshold_new ans

% save workspace (without mask variables which take a long time)
save(strcat(analysis,'/03_CellPropertyCalculator.mat'), '-regexp', '^(?!(bwf|bwf_5)$).');


%% 06) [PART 1] Cell Property Collector & Histogram Generator
% collect cell area of all time for each xy position 
% visualize with histogram -> helps with determination of thresholds

% saveplace = strcat(analysis,'/mask plots'); mkdir(saveplace); % where updated masks will be saved

% collect all cell property values into structures for generating histograms
expt_cond_range = analyzed_q;
xy_range = analyzed_xy;

[Q75_cum_xy, area_cum_xy area_cum Q75_cum] = ...
        CellPropertyCollector(xy_range, total_time, expt_cond_range, expt_conditions_xy, area, Q75);

% generate histograms
Q75_hist = 1; % 1 = make histogram of UQ
a_hist = 1; % 1 = make histogram of area
bar_on = 0; % 1 = bar on histogram on; 0 = bar on histogram off (line only)
nbins_q = 100;
nbins_a = 100;

[q a] = CellPropertyHistogramGenerator(expt_cond_range,nbins_q,nbins_a,area,Q75,area_cum,Q75_cum,...
                                       expt_conditions_string, Q75_hist, a_hist, bar_on);

                                   
%% 06) [PART 2] adjust axis
xlim(a,[0 200]);
xlim(q,[600 1200]);


%% 06) [PART 3] save figures
saveas(a,strcat(saveplace,'/area_cells_histogram.png'));
saveas(q,strcat(saveplace,'/Q75_cells_histogram.png'));

% save workspace (without mask variables which take a long time)
clearvars saveplace Q75_hist a_hist bar_on q a ans xy_range expt_cond_range; close all;
save(strcat(analysis,'/04_CellPropertyCollector.mat'), '-regexp', '^(?!(bwf|bwf_5)$).');


%% 07) thresholds on area and UQ
    % I) iterate through thresholds on area and UQ - COARSE GRAIN ITERATION
    % II) real analysis (takes long)

% saveplace = strcat(analysis,'/threshold_1'); mkdir(saveplace); % where updated masks will be saved
% cd(wrkDir); % set directory to where TIFF images are located

close all;
RealAnalysis = 0;
TimeDivide = 0; % 1 if different thresh used for different time points (0 for coarse grain analysis)
redo = 0; % 1 if re-evaluating failed masks

% default save flags (dont change here; flags are changed in the if statements below)
FourPanelFigure = 0; % 4-panel mask image
hist_on = 0; % make histogram to visualize threshold or not
SaveFlag_hist = 0; % 1 = histogram figure suppressed and saved; 0 = figure pops up but not saved;
SingleMaskFigure = 0; % 1 for real analysis; 2-panel mask image
SaveFlag_mask = 0; % 1 for real analysis; 1 = mask figure suppressed and saved; 0 = figure pops up but not saved;

if RealAnalysis == 0;
    % set save flags
    FourPanelFigure = 1; % 4-panel mask image
    hist_on = 1; % make histogram to visualize threshold or not
    SaveFlag_hist = 1; % 1 = histogram figure suppressed and saved; 0 = figure pops up but not saved;
    
    q = 1; % MANUAL CHANGE for test; need to run through all q
    q_range = q; % for for-loop
    xy_range = expt_conditions_xy{q}(3) % an example xy position corresponding to the q
    t_range = [1 15]%[15:total_time]; % representative time points for threshold testing

    % determine thresholds
    percent_max_Q = 0; % percentage of cells to eliminate from the top
    thresh_min_area = 10;
    thresh_max_area = 200;
    thresh_max_Q(q) = 20000; % to specify specific value
    
    % calculate max UQ threshold for each expt condition (q) -> save final value in thresh_max_Q vector
    num_thresh_Q_max = round((percent_max_Q * length(Q75_cum{q})),0); % number of cells that are in threshold % of total 
    Q75_sorted = sort(Q75_cum{q},'descend'); % order Q75 from large to small
    
    % record
    thresh_percent_max_Q(q) = percent_max_Q; % record the percent
%    thresh_max_Q(q) = Q75_sorted(num_thresh_Q_max); % extract the Q75 value after the top % threshold    
    thresh_min_a(q) = thresh_min_area;
    thresh_max_a(q) = thresh_max_area; 
end
clearvars num_thresh_Q_max Q75_sorted

if RealAnalysis == 1 & redo == 0; % first time real analysis
    % set save flags
    SingleMaskFigure = 1; % 1 for real analysis; 2-panel mask image
    SaveFlag_mask = 1; % 1 for real analysis; 1 = mask figure suppressed and saved; 0 = figure pops up but not saved;

    if TimeDivide == 1;
        t_range1 = 1; % earlier time points
            thr_min_a1 = 10;
            thr_max_a1 = 180;
            thr_max_Q1 = 510;
            thr_percent_max_Q1 = 0.03;
        t_range2 = 2:5; % earlier time points
            thr_min_a2 = 10;
            thr_max_a2 = 180;
            thr_max_Q2 = 481;
            thr_percent_max_Q2 = 0.18;
    end
    
    t_range = 1:total_time;
    q_range = analyzed_q;
    
    h = waitbar(0,'out of total xy positions'); % waitbar only for real analysis
    
elseif RealAnalysis == 1 & redo == 1;
    
    % set save flags
    SingleMaskFigure = 1; % 1 for real analysis; 2-panel mask image
    SaveFlag_mask = 1; % 1 for real analysis; 1 = mask figure suppressed and saved; 0 = figure pops up but not saved;

    % MANUAL input on which images to re-do
    q_range = [4:8]; % dummy variable to let the for-loop proceed; the real variable that matters is xy_range_failed
    t_range_failed = [1:total_time];
    xy_range_failed = [22:56];
end

for q = q_range;
    
    if RealAnalysis == 1 & redo == 0;
        xy_range = expt_conditions_xy{q};
    elseif RealAnalysis == 1 & redo == 1;
        xy_range = xy_range_failed; % manual entry from above
        t_range = t_range_failed;
    end
    
    for xy = xy_range; 
            if RealAnalysis == 1 & redo == 0;
                waitbar(xy / total_xy);
            elseif RealAnalysis == 1 & redo == 1; % print xy since no waitbar
                xy
            end
            
        for t = t_range;
            
        % if there is no image corresponding to the t,xy
        if isempty(bfImgs(t,xy).name) == 1;
            continue            
        end
        
        % load image & set input variables
        phasename = bfImgs(t,xy).name;
        bwf_txy = bwf{t,xy};
        cells_txy = cells{t,xy};
        area_cell_temp = cell2mat(struct2cell(area{t,xy}));
        Q75_cell_temp = Q75{t,xy}; % UQ75 values of all cells in that image
        
        % set threshold - if dividing images by time or not
        if RealAnalysis == 1 & TimeDivide == 1;
            if ismember(t,t_range1); % early time point
               thr_min_a = thr_min_a1;
               thr_max_a = thr_max_a1;
               thr_max_Q = thr_max_Q1;
               thr_percent_max_Q = thr_percent_max_Q1;
            elseif ismember(t,t_range2); % later time point
               thr_min_a = thr_min_a2;
               thr_max_a = thr_max_a2;
               thr_max_Q = thr_max_Q2;
               thr_percent_max_Q = thr_percent_max_Q2;          
            end
        elseif RealAnalysis == 0 | (RealAnalysis == 1 & TimeDivide == 0);
                thr_min_a = thresh_min_a(q);
                thr_max_a = thresh_max_a(q);
                thr_max_Q = thresh_max_Q(q);
                thr_percent_max_Q = thresh_percent_max_Q(q);
        end
        
                
        % run script
        [a_elim(t,xy) Q_elim(t,xy) cells_post_thresh{t,xy}...
         bwf_thresh_a{t,xy} bwf_thresh_a_elim{t,xy} ...
         bwf_thresh_Q{t,xy} bwf_thresh_Q_elim{t,xy} ...
         bwf_post_thresh{t,xy} bwf_post_thresh_elim{t,xy}] = ...
            ThresholdIterator(q, phasename, bwf_txy, cells_txy, thr_min_a, thr_max_a, thr_max_Q, thr_percent_max_Q,...
            area_cell_temp, Q75_cell_temp, ...
            SaveFlag_mask, SaveFlag_hist, RealAnalysis, FourPanelFigure, SingleMaskFigure, saveplace, hist_on, Q75_cum{q}, area_cum{q});

        end % end of t
    end % end of xy
end % end of q
clearvars ans q xy t xy_range t_range thresh_min_area percent_max_Q
clearvars phasename bwf_txy cells_txy area_cell_temp Q75_cell_temp
clearvars FourPanelFigure SingleMaskFigure SaveFlag_mask SaveFlag_hist hist_on TimeDivide
clearvars t_range1 t_range2 thr_min_a thr_max_Q thresh_max_area thr_percent_max_Q

% save workplace
if RealAnalysis == 1 & redo == 0;
   close(h)
   save(strcat(analysis,'/05_ThresholdIterator_1.mat'), '-regexp', '^(?!(bwf|bwf_5)$).');
   save(strcat(analysis,'/bwf_post_thresh.mat'),'bwf_post_thresh'); % save the final masks variable
   save(strcat(analysis,'/cells_post_thresh.mat'),'cells_post_thresh'); % save the final masks variable
   clearvars h RealAnalysis redo
elseif RealAnalysis == 1 & redo == 1;
   save(strcat(analysis,'/05_ThresholdIterator_2.mat'), '-regexp', '^(?!(bwf|bwf_5)$).');
   save(strcat(analysis,'/bwf_post_thresh.mat'),'bwf_post_thresh'); % save the final masks variable
   save(strcat(analysis,'/cells_post_thresh.mat'),'cells_post_thresh'); % save the final masks variable
   clearvars RealAnalysis redo xy_range_failed t_range_failed q_range
end


%% 08) make thresholded mask figures for visual inspection

% visually check ensure mask formation success
close all;
SaveIm = 1; % 1 = save image; 0 = don't save image, set to 0 for redo if only a few xy positions
redo = 1; % make sure that SaveIm = 0 if redo = 1 and only a few xy positions

cells_elim_color = colormap(copper(2)); % value 2 gives the highest contrast
subplot_row = 4;
subplot_col = 8; % usually 8
total_im_per_subplot = subplot_row * subplot_col; 
xy_figure = [1:total_im_per_subplot:total_xy];

% define xy range for the loop
if redo == 0;
    xy_range_loop = [1:total_xy total_xy+1]; % 1 larger than total_xy
elseif redo == 1;
    xy_range = [22:56]; % MANUAL ENTRY for redo
    xy_range_loop = [xy_range];
end


fig_count = 0; % initialize fig_count (always start with 0)
if SaveIm == 1 & redo == 0;
     h = waitbar(0,'out of total xy positions'); steps = total_xy; % waitbar (but not for redo)
end
for xy = xy_range_loop;
    for t = 1:total_time;
        if SaveIm == 1 & redo == 0;
            waitbar(xy / total_xy);
        elseif SaveIm == 1 & redo == 1;
            xy % simply print xy if redo
        end    
            
    % SAVE FIGURE after the first image
    if (t == 1 & xy > xy_range_loop(1) & SaveIm == 1) | (xy == total_xy+1 & SaveIm == 1);
       % save figure (high resolution)
        print(f(fig_count),strcat(saveplace,'/','maskoverlay_thresh_',num2str(fig_count),'.png'),'-dpng','-r300');           
            if xy == total_xy+1;
                close(h)
                return
            end
    end
      
    % every multi-plot images, create a new figure
    if t == 1;
        i = 1; % re-initialize i for subplot
        fig_count = fig_count + 1; % advance figure count   
        f(fig_count) = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        if SaveIm == 1;
           set(f(fig_count), 'Visible', 'off'); 
        end  
    end
       
    % load phase image
    phasename = bfImgs(t,xy).name;
    ph_adj_temp = imread(phasename);
    ph_adj_temp_vis = imadjust(ph_adj_temp);
    BWoutline_final = imdilate(double(bwperim(bwf_post_thresh{t,xy})),strel('disk',1)); % define the thickness of outline

    % create subplot_row x subplot_col matrix of images
    subplot(subplot_row,subplot_col,i);   
    
    imshowpair(ph_adj_temp_vis, BWoutline_final);
         hold on;
    bwf_a_Q_elim = imshow(bwf_post_thresh_elim{t,xy}, cells_elim_color);
         set(bwf_a_Q_elim, 'AlphaData', 0.45); % eliminated cells look like they are glowing
    title(phasename);
             
    % tweak the positioning and size of each subplot
    original_pos = get(gca, 'Position'); % Position = [left bottom width height]
    %new_pos = original_pos + [-0.05 -0.05 0.06 0.06]; % good for 4x6
    %new_pos = original_pos + [-0.04 -0.04 0.045 0.045]; % good for 5x7
    new_pos = original_pos + [-0.04 -0.04 0.03 0.03]; % good for 4x8
    set(gca, 'Position', new_pos);
    
    % advance i for subplot position
    i = i + 1;
    
    end
end

if SaveIm == 1;
    clearvars ans cells_elim_color f fig_count i new_pos original_pos
    clearvars ph_adj_temp ph_adj_temp_vis phasename
    clearvars SaveIm saveplace subplot_col subplot_row xy xy_figure
end


%% 09) visualize the number of eliminated cells

% compile by expt condition
for q = analyzed_q;    
    xy_temp = expt_conditions_xy{q};
    
    % first compute elimination ratio (normalized to total detected cells
    for t = 1:total_time;        
        for xy = xy_temp;
            if isempty(bfImgs(t,xy).name) == 1;
                continue
            end
            
            numCells = cells{t,xy}.NumObjects;
            a_elim_value = a_elim(t,xy);
            Q_elim_value = Q_elim(t,xy);
            a_elim_ratio(t,xy) = a_elim_value / numCells;
            Q_elim_ratio(t,xy) = Q_elim_value / numCells;
        end
    end
    
    % second compute mean and standard deviation for each expt cond
    for t = 1:total_time;   
            a_elim_ratio_cond(t,q) = mean(a_elim_ratio(t,xy_temp));
            a_elim_ratio_cond_std(t,q) = std(a_elim_ratio(t,xy_temp));
            Q_elim_ratio_cond(t,q) = mean(Q_elim_ratio(t,xy_temp));
            Q_elim_ratio_cond_std(t,q) = std(Q_elim_ratio(t,xy_temp));   
    end
end
clearvars q xy_temp t xy numCells a_elim_value Q_elim_value

% plot line colors
colors = hsv(length(analyzed_q));

% plot number of eliminated cells over time
a = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
for q = analyzed_q;
    i = find(q==analyzed_q); % get the index
    p(i) = errorbar(a_elim_ratio_cond(:,q),a_elim_ratio_cond_std(:,q),'.-',...
                    'MarkerSize',30,...
                    'LineWidth',2,...
                    'Color',colors(i,:));
    hold on
end
hold off; 
xlabel('time point','FontSize',16);
ylabel('number of eliminated cells normalized by total cells detected','FontSize',16);
title('number of eliminated cells due to min & max cell area threshold','FontSize',16);
legend(p(:),{expt_conditions_string{analyzed_q}},'FontSize',16,'Location','NorthEast'); % 1 for bars, 2 for line
% axis([0 22 -0.05 0.2]);

% Q elim
uq = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
for q = analyzed_q;
    i = find(q==analyzed_q); % get the index
    p(i) = errorbar(Q_elim_ratio_cond(:,q),Q_elim_ratio_cond_std(:,q),'.-',...
                    'MarkerSize',30,...
                    'LineWidth',2,...
                     'Color',colors(i,:));
    hold on
end
hold off
xlabel('time point','FontSize',16);
ylabel('number of eliminated cells normalized by total cells detected','FontSize',16);
title('number of eliminated cells due to maximum upper quartile pixel intensity threshold','FontSize',16);
legend(p(:),{expt_conditions_string{analyzed_q}},'FontSize',16,'Location','NorthEast'); % 1 for bars, 2 for line
% axis([0 22 0 0.4]);


%% 09) [PART 2] save figures
saveplace = strcat(analysis,'/threshold plots'); mkdir(saveplace)
saveas(a,strcat(saveplace,'/cells_elim_area_threshold.png'));
saveas(uq,strcat(saveplace,'/cells_elim_UQ_threshold.png'));

clearvars a uq i p q colors
close all


%% 10) [PART 1] create background mask
% STEP 1: test thresholds and dilution radius

% cd(wrkDir)
close all;
% STEP 1: assess pixel histogram
q = 9; % MANUAL CHANGE
xy = expt_conditions_xy{q}(3);
t = 5;

% record thresholds
bckg_thresh_min(q) = 800;
bckg_thresh_max(q) = 1000;
dil_radius(q) = 30; % set the radius of inflation (usually ~10 pixels)

% load images
ph_adj_temp = imadjust(imread(bfImgs(t,xy).name)); % load a dummy phase image
ph = imread(bfImgs(t,xy).name); % load a dummy phase image

yfp_adj_temp = imadjust(imread(YFP(t,xy).name));

% form test background mask
[bckg_mask_test{t,xy}] = BckgMaskGenerator(bfImgs(t,xy).name,bwf_5{t,xy},dil_radius(q),bckg_thresh_min(q),bckg_thresh_max(q));

% visualize histogram and background mask
figure('units', 'normalized', 'outerposition', [0 0 1 1]); 
   % subplot(221); histogram(ph_adj_temp);
    subplot(221); histogram(ph);
        hold on;
        plot(bckg_thresh_min(q),0,'r*','MarkerSize',24); 
        plot(bckg_thresh_max(q),0,'r*','MarkerSize',24);
    subplot(222); imshow(yfp_adj_temp.*uint16(bckg_mask_test{t,xy})); title('bckg_mask YFP (final)');
    subplot(223); imshow(ph_adj_temp.*uint16(bckg_mask_test{t,xy})); title('bckg_mask phase (final)');
    subplot(224); imshow(ph_adj_temp.*uint16(~bckg_mask_test{t,xy})); title('bckg_mask phase(inverted, final)');
clearvars yfp_adj_temp ph_adj_temp q xy t
    

%% 10) [PART 2] create background mask
% STEP 2: generate background mask (takes some time)

close all;
tdivide = 0; % 1 = different thresholds used for different t

h = waitbar(0,'out of total experimental conditions, q');

for q = analyzed_q;
    waitbar(q / total_cond);
    for xy = expt_conditions_xy{q};       
        for t = 1:total_time;
                   
            if isempty(bfImgs(t,xy).name) == 1;
                continue
            end
            
            % if different thresholds used for different t, update this section
            if tdivide == 1;            
                if t == 1; % manual change
                    dil_r = 20; % manual change
                    bckg_thr_min = bckg_thresh_min(q); % manual change
                    bckg_thr_max = bckg_thresh_max(q); % manual change
                else 
                    dil_r = 15;
                    bckg_thr_min = bckg_thresh_min(q);
                    bckg_thr_max = bckg_thresh_max(q);
                end
            else % not divided by time
                dil_r = dil_radius(q);
                bckg_thr_min = bckg_thresh_min(q);
                bckg_thr_max = bckg_thresh_max(q);
            end
                     
        % specify variables
        phasename = bfImgs(t,xy).name;
        mask_cell = bwf_5{t,xy}; % cell mask, before border clearance

        % run scr420ipt
        [bckg_mask{t,xy}] = BckgMaskGenerator(phasename, mask_cell, ...
                                dil_r, bckg_thr_min, bckg_thr_max);
                            
        end
    end
end
close(h)
clearvars ans h q xy t phasename mask_cell bckg_mask_test 
clearvars ph_adj_temp yfp_adj_temp
clearvars dil_r bckg_thr_min bckg_thr_max tdivide


%% 10) [PART 3] create background mask
% STEP 3: test background mask visualize

t = 2; xy = 3; % load a dummy image
phase = imread(bfImgs(t,xy).name);
% yimg = imadjust(imread(YFP(t,xy).name));
% timg = imadjust(imread(TFP(t,xy).name));
figure('units', 'normalized', 'outerposition', [0 0 1 1]); 
subplot(121); imshow(imadjust(phase).*uint16(bckg_mask{t,xy})); title('bckg_mask phase (final)');
subplot(122); imshow(imadjust(phase).*uint16(~bckg_mask{t,xy})); title('bckg_mask phase (inverted, final)');

clearvars ans t xy phase img timg


%% 10) [PART 4] create background mask
% step 4: save background masks

close all;
save(strcat(analysis,'/bckg_mask.mat'),'bckg_mask','dil_radius','bckg_thresh_min','bckg_thresh_max'); % overwrite previous version


%% 11) [PART 1] calculate fluorescence

% clearvars -except motherdir expt_list 
% expt = 6 % MANUAL change

% reset directories
% cd(motherdir) % containing all DMSPkinetics folders
% expt_folder_name = dir(strcat(expt_list{expt}(1:10),'_*')) % experiment date folder name
% data_folder = strcat(expt_folder_name.name,'/analysis/DMSPkinetics') % experiment date folder directory
% cd(data_folder)

% load variables
% load('01_ParsedFileNames.mat', 'bfImgs','RFP','YFP','TFP','analyzed_xy','total_time','wrkDir');
% load('bckg_mask.mat')
% load('wrkDir.mat'); % if wrkDir had changed between PC and Mac, load the updated wrkDir variable 
% cd('shifted_mask_fix'); % folder containing the updated masks
% load('bwf_cells_post_thresh_fixed_final.mat');
% analysis = pwd % set current folder as place to save analysis file and raw fluo plots

% redefine variables for for-loop
% bwf_post_thresh = bwf_post_thresh_fixed; 
% cells_post_thresh = cells_post_thresh_fixed;

% for 5/16/2017 expt only
bwf_post_thresh = bwf;
cells_post_thresh = cells;

% original 
close all; cd(wrkDir)

h = waitbar(0,'out of total xy positions');
for xy = analyzed_xy;
    waitbar(xy / length(analyzed_xy));
    for t = 1:total_time;
        
        if isempty(bfImgs(t,xy).name) == 1;
            continue
        end
        
        if isempty(bckg_mask{t,xy}) == 1;
            continue
        end
        
     % specify variables
     rfp = double(imread(RFP(t,xy).name)); % originally kept as uint16
     yfp = double(imread(YFP(t,xy).name)); % originally kept as uint16 
     tfp = double(imread(TFP(t,xy).name)); % originally kept as uint16 
     c_mask = bwf_post_thresh{t,xy}; % cell mask
     b_mask = bckg_mask{t,xy}; % background mask
     cellpixelID = cells_post_thresh{t,xy}.PixelIdxList; % PixelIdxList
        
    [rfp_bckg(t,xy),      yfp_bckg(t,xy),      tfp_bckg(t,xy), ...
     rfp_bckg_std(t,xy),  yfp_bckg_std(t,xy),  tfp_bckg_std(t,xy), ...
     rfp_fintcells{t,xy}, yfp_fintcells{t,xy}, tfp_fintcells{t,xy}, ...
     rfp_fintview(t,xy),  yfp_fintview(t,xy),  tfp_fintview(t,xy), ...
     rfp_std(t,xy),       yfp_std(t,xy),       tfp_std(t,xy), ...
     rfp_sem(t,xy),       yfp_sem(t,xy),       tfp_sem(t,xy), ...
     totalpixelnum_bckg(t,xy)] = FluoCellComputation(rfp, yfp, tfp, c_mask, b_mask, cellpixelID);

    end
end
close(h)

clearvars ans h t xy rfp yfp tfp c_mask b_mask

% save workspace
save(strcat(analysis,'/06_FluoCellComputation_all.mat'), '-regexp', ...
    '^(?!(bwf|bwf_5|backg_mask|bwf_post_thresh|bwf_post_thresh_elim|bwf_thresh_a|bwf_thresh_a_elim|bwf_thresh_Q|bwf_thresh_Q_elim)$).');


%% 11) [PART 2] quick visualization of fluorescence analysis results

% old script:
    % xy = 40; figure;
    % plot(1:total_time,rfp_fintview(:,xy),'r')
    % hold on
    % plot(1:total_time,tfp_fintview(:,xy),'b')

    
% ADDED 12/21/2018 (after aggdata fix)
    % visualize all individual xy positions
    % visually check png images to note any unusual kinks
        % if any unusual kinks, fix with "fix_imanalysis_in_aggdata.m"
saveplace = strcat(analysis,'/raw fluo kinetics check'); mkdir(saveplace)
for xy = analyzed_xy
    
    figure('visible','off');
    
    % plot all colors
    plot(1:total_time, rfp_fintview(:,xy), 'ro-'); hold on
    plot(1:total_time, yfp_fintview(:,xy), 'ko-'); 
    plot(1:total_time, tfp_fintview(:,xy), 'bo-');
    
    % figure titles
    title(strcat('xy position : ', num2str(xy), ';   raw fluorescence kinetics in each color'))
    xlabel('time point')
    ylabel('fluo intensity (black = YFP; red = RFP; blue = TFP')
    
    % save image (low quality)
    print(strcat(saveplace,'/','raw_fluo_check_xy_',num2str(xy),'.png'),'-dpng','-r100');    
   
end

soundsc(1)

%% 12) extract other information - number of cells

% number of cells by xy position
for t = 1:total_time;
    for xy = 1:total_xy;       
        cell_num(t,xy) = length(yfp_fintcells{t,xy});
    end
end
% average number of cells by expt condition
for q = 1:total_cond;
    xy_range = expt_conditions_xy{q};
    for t = 1:total_time;        
        cell_num_avg(t,q) = mean(cell_num(t,xy_range));
        cell_num_std(t,q) = std(cell_num(t,xy_range));
    end
end
clearvars ans t xy q xy_range
save(strcat(analysis,'/07_CellNumber_preYFPthresh.mat'),'cell_num')%,'cell_num_avg','cell_num_std');


%% 13) extract other information - cell size

cell_size = {}; % initialize
for t = 1:total_time;
    for xy = 1:total_xy;       
        
        if isempty(bfImgs(t,xy).name) == 1;
            continue
        end
        
        % pull out cell info of each image
        cells_image = cells_post_thresh{t,xy}.PixelIdxList; 
        % initialize temporary container        
        size_temp = []; 
            
            % extract number of pixels in each cells
            for cell = 1:length(cells_image);
                size_temp(cell) = length(cells_image{cell});
            end
        
        % store in permanent container
        cell_size{t,xy} = size_temp;
        
        clearvars cells_image % initialize
    end
end
clearvars size_temp cells_image t xy 
save(strcat(analysis,'/08_CellSize_preYFPthresh.mat'),'cell_size');


%% 14) final save

% save only working variables
clearvars -except   rfp_bckg      yfp_bckg       tfp_bckg ...
                    rfp_bckg_std  yfp_bckg_std   tfp_bckg_std ...
                    rfp_fintcells yfp_fintcells  tfp_fintcells ...
                    rfp_fintview  yfp_fintview   tfp_fintview ...
                    rfp_std       yfp_std        tfp_std ...
                    rfp_sem       yfp_sem        tfp_sem ...
                    totalpixelnum_bckg ...
                    total_channel total_cond total_time total_xy total_xy_per_cond ...
                    wrkDir analysis bfImgs RFP YFP TFP ...
                    expt_conditions_string expt_conditions_xy ...
                    analyzed_q analyzed_xy expt_date time_expt_hours ...
                    cell_num cell_num_avg cell_num_std cell_size 
                    
                
save(strcat(analysis,'/06_FluoCellComputation_analysisonly.mat'));


%% NEXT STEP: SPECTRAL LEAKAGE CORRECTION using aggregated b1

