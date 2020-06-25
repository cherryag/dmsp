%% Image Analysis Pipeline (no time component)
% version November 2018

% Recognize cells and quantify fluorescence values for each cell

% 01) parse file names
% 02) recognize cells by pixel intensity [also iterate for failed masks] (CellRecognizer.m)
% 03) create multi-panel phase-mask images for inspection
% 04) [optional] enter manually images to exclude, if any (e.g. due to evaporation)
% 04) iterate new thresholds for failed masks (CellRecognizer.m)
% 05) segment cells and calculate area and UQ (CellPropertyCalculator.m)
% 06) Cell Property Collector & Histogram Generator (CellPropertyCollector.m & CellPropertyHistogramGenerator.m)
% 07) thresholds on area and UQ [also iterate for failed masks] (ThresholdIterator.m)
% 08) make thresholded mask figures for visual inspection
% 09) visualize the number of eliminated cells
% 10) create background mask [PARTS 1-4] (BckgMaskGenerator.m)
% 11) fluorescence calculation (FluoCellComputation.m)
% 12) extract other information - number of cells
% 13) extract other information - cell size
% 14) FINAL SAVE
% NEXT STEP: 03_SpectralLeakageCorrector_FPctrl.m


%% define directories 
% reverse slash direction for Mac

wrkDir = uigetdir(); wrkDir = strrep(wrkDir,'\','/'); cd(wrkDir); % folder with TIFF images
analysis = uigetdir(); analysis = strrep(analysis,'\','/'); % where analysis files will be saved
expt_date = '8/3/2018 FP controls'; % for titles


%% 01) parse file names & expt conditions

imgs = dir('*.tif');
lastImg = imgs(end).name
parsedFileName = strsplit(lastImg,{'FPctrl_','xy','c','.tif'})
total_xy = str2num(parsedFileName{2}) % total number of positions
total_channel = str2num(parsedFileName{3}) % total number of channels

% matrices of TIFF image file names with t rows x p columns 
for xy = 1:total_xy;     
    bfImgs(xy) = dir(strcat('*xy',num2str(xy,'%02i'),'c1.tif'));
    RFP(xy) = dir(strcat('*xy',num2str(xy,'%02i'),'c2.tif'));
    YFP(xy) = dir(strcat('*xy',num2str(xy,'%02i'),'c3.tif'));
    TFP(xy) = dir(strcat('*xy',num2str(xy,'%02i'),'c4.tif'));
end

clearvars imgs lastImg parsedFileName xy ans

expt_conditions_string = {'PA-mKate',...
                          'PA-YFP',...
                          'PA-TFP'}; % for legend  
                      
expt_conditions_xy = {[1:21] [22:42] [43:63]}; % MANUAL CHANGE

total_cond = length(expt_conditions_xy);
analyzed_q = 1:total_cond;
analyzed_xy = 1:total_xy; 
                      
% save workplace in the same folder as TIFF files
save(strcat(analysis,'/01_ParsedFileNames.mat')); clearvars saveplace;

%% 02) recognize cells by pixel intensity (CellRecognizer.m)
% run script CellRecognizer.m
% this step takes a long time if save mask overlay images
% FAILED MASKS: also iterate for failed masks by setting SaveFlag == 2

% cd(wrkDir);
% saveplace = strcat(analysis,'/masks_1'); mkdir(saveplace);

close all;
threshold = 447; % on non-adjusted image
SaveFlag = 1; % 0 = for setting thresholds; 1 = first iteration; 2 = redo for failed masks
thresh_divide = 1; % 1 = if there are different thresholds for different q; 0 = same threshold for all

% determine which images to run through for-loop
if SaveFlag == 0; % for tests, only look at representative images    
    q = 2; % MANUAL CHANGE
    xy_range = expt_conditions_xy{q}(6);    
elseif SaveFlag == 1; % all images   
    if thresh_divide == 1; 
        q1 = [1 2]; % q to be paired with threshold1 
            threshold1 = 447; 
        q2 = [3]; % q to be paired with threshold2
            threshold2 = 442;
        q3 = [0]; % q to be paired with threshold3
            threshold3 = 0;
    end
    xy_range = analyzed_xy;
    bwf_5 = {}; % initialize mask container, BEFORE border clearing
    bwf = {}; % initialize mask container, FINAL mask
    
    h = waitbar(0,'out of total xy points');
elseif SaveFlag == 2; % redo
    xy_range = 36:47;
    threshold = 490; 
end

for xy = xy_range;
         
    if SaveFlag == 1;
       waitbar(xy / length(xy_range));  
       
       if thresh_divide == 1; 
            if ismember(xy,[expt_conditions_xy{q1}]) == 1;
                    threshold = threshold1;
                elseif ismember(xy,[expt_conditions_xy{q2}]) == 1;
                    threshold = threshold2;
                elseif ismember(xy,[expt_conditions_xy{q3}]) == 1;
                    threshold = threshold3;
                else threshold = threshold1;
            end
       end % end of thresh_divide
       
    end
       
    phasename = bfImgs(xy).name;
    [bwf_5{xy} bwf{xy}] = CellRecognizer(SaveFlag,saveplace,phasename,threshold);

end

% save workspace
if SaveFlag == 1;
    close(h)
    save(strcat(analysis,'/02_CellRecognizer_1.mat')); % includes bwf
    save(strcat(analysis,'/bwf_final.mat'),'bwf'); % save the final masks variable
    save(strcat(analysis,'/bwf_5_final.mat'),'bwf_5'); % save the final masks variable
    clearvars h SaveFlag q1 q2 q3 thresh_divide SaveFlag 
    clearvars threshold1 threshold2 threshold3 threshold
elseif SaveFlag == 2;
    save(strcat(analysis,'/02_CellRecognizer_2.mat')); % includes bwf
    save(strcat(analysis,'/bwf_final.mat'),'bwf'); % save the final masks variable
    save(strcat(analysis,'/bwf_5_final.mat'),'bwf_5'); % save the final masks variable
    clearvars SaveFlag q1 q2 q3 thresh_divide SaveFlag 
    clearvars threshold1 threshold2 threshold3 threshold
end



%% 03) create multi-panel phase-mask images for inspection
% visually check ensure mask formation success
cd(wrkDir)
close all;
SaveIm = 1; % 1 = save image

subplot_row = 4;
subplot_col = 6;
total_im_per_subplot = subplot_row * subplot_col; 
xy_figure = [1:total_im_per_subplot:total_xy];

fig_count = 0; % initialize fig_count (always start with 0)
for xy = [xy_range xy_range(end)+1];  % add one more xy value s.t. include the last image
    
    % SAVE FIGURE after the first image
    if ismember(xy,xy_figure(2:end)) | xy == analyzed_xy(end)+1 & SaveIm == 1;
       % save figure (high resolution)
       print(f(fig_count),strcat(saveplace,'/','maskoverlay_',num2str(fig_count),'.png'),'-dpng','-r300');
           if xy == analyzed_xy(end)+1;
                % clean up workspace at the end of analysis
                clearvars xy xy_range q_range ans thresh_divide q q1 q2 q3
                clearvars threshold threshold1 threshold2 threshold3 phasename  
                clearvars f fig_count i new_pos original_pos ph_adj_temp ph_adj_temp_vis
                clearvars saveplace SaveIm subplot_col subplot_row total_im_per_subplot xy_figure
               return % end of images; get out of the for-loop
           end
    end
    
    % every (subplot_row x subplot_col) images, create a new figure
    if ismember(xy,xy_figure) == 1;
        i = 1; % re-initialize i for subplot
        fig_count = fig_count + 1; % advance figure count   
        f(fig_count) = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        if SaveIm == 1;
           set(f(fig_count), 'Visible', 'off'); 
        end  
    end
       
    % load phase image
    phasename = bfImgs(xy).name;
    ph_adj_temp = imread(phasename);
    ph_adj_temp_vis = imadjust(ph_adj_temp);
   
    % create (subplot_row x subplot_col) matrix of images
    subplot(subplot_row,subplot_col,i);
        imshowpair(ph_adj_temp_vis, bwf{xy}); hold on;
        title(phasename) 
        
    % tweak the positioning and size of each subplot
    original_pos = get(gca, 'Position'); % Position = [left bottom width height]
    new_pos = original_pos + [-0.05 -0.05 0.06 0.06]; % good for 4x6
    %new_pos = original_pos + [-0.04 -0.04 0.045 0.045]; % good for 5x7
    set(gca, 'Position', new_pos);
    
    % advance i for subplot position
    i = i + 1;
       
end

%% 04) [optional] enter manually images to exclude, if any (e.g. due to evaporation)

xy_exclude = [12,30,31];

for xy = xy_exclude;
    bfImgs(xy).name = [];
    bwf{xy} = []; % make mask container empty
    bwf_5{xy} = []; % make mask container empty
end
clearvars xy_exclude t xy


%% 05) segment cells and calculate area and UQ (CellPropertyCalculator.m)

% initialize containers
Q75 = {}; area = {};

h = waitbar(0,'out of total xy positions');
for xy = 1:total_xy;     
    waitbar(xy / total_xy);   

    phasename = bfImgs(xy).name;
    bwf_txy = bwf{xy};
    
    % skip excluded images
    if isempty(bfImgs(xy).name) == 1;
        continue            
    end

    [cells{xy} area{xy} Q75{xy}] = CellPropertyCalculator(bwf_txy, phasename);

end
close(h)
clearvars xy phasename bwf_txy h steps ans saveplace threshold_new ans

% save workspace (without mask variables which take a long time)
save(strcat(analysis,'/03_CellPropertyCalculator.mat'), '-regexp', '^(?!(bwf|bwf_5)$).');


%% 06) [PART 1] Cell Property Collector & Histogram Generator (CellPropertyCollector.m & CellPropertyHistogramGenerator.m)
% collect cell area and UQ for all xy positions.
% visualize with histogram -> helps with determination of thresholds.
% CellPropertyCollector.m is not suitable for this without time points.

% saveplace = strcat(analysis,'/mask plots'); mkdir(saveplace); % where updated masks will be saved

% collect all cell property values into structures for generating histograms
for expt_cond = 1:total_cond; % number of conditions
    
    cum_temp_a = [];  cum_temp_Q = []; % initialize
    cond_xy = expt_conditions_xy{expt_cond};
    
    for p = 1:length(cond_xy); % cycle through the xy positions of the expt cond       
        temp_a = []; temp_Q = []; % initialize
        
        xy_index = cond_xy(p); % call up the xy index
            % if there is no image corresponding to the t,xy
            if isempty(bfImgs(xy_index).name) == 1;
                continue            
            end
        temp_a = [area{xy_index}.Area]; % for area     
        temp_Q = Q75{xy_index}; % for Q75 
        
        cum_temp_a = [cum_temp_a; temp_a'];
        cum_temp_Q = [cum_temp_Q; temp_Q]; % append
             
    end
    
    area_cum{expt_cond} = cum_temp_a; % for area
    Q75_cum{expt_cond} = cum_temp_Q; % for Q75
    
end
clearvars expt_cond p cond_xy cum_temp_a cum_temp_Q temp_a temp_Q xy_index ans

% run script to generate histograms
Q75_hist = 1; % 1 = make histogram of UQ
a_hist = 1; % 1 = make histogram of area
bar_on = 0; % 1 = bar on histogram on; 0 = bar on histogram off (line only)
expt_cond_range = 1:total_cond;
nbins_q = 50;
nbins_a = 50;

[q a] = CellPropertyHistogramGenerator(expt_cond_range,nbins_q,nbins_a,area,Q75,area_cum,Q75_cum,expt_conditions_string, Q75_hist, a_hist, bar_on);

%% 06) [PART 2] adjust axes of plots
xlim(a,[0 100]); 
xlim(q,[380 500]);

%% 06) [PART 3] save figures
saveas(a,strcat(saveplace,'/area_cells_histogram.png'));
saveas(q,strcat(saveplace,'/Q75_cells_histogram.png'));

% save workspace (without mask variables which take a long time)
clearvars saveplace Q75_hist a_hist bar_on q a ans expt_cond_range nbins_a nbins_q; close all;
save(strcat(analysis,'/04_CellPropertyCollector.mat'), '-regexp', '^(?!(bwf|bwf_5)$).');


%% 07) thresholds on area and UQ [also iterate for failed masks] (ThresholdIterator.m)
    % 1) iterate through thresholds on area and UQ - COARSE GRAIN ITERATION
    % 2) REAL analysis (takes long)
    % 3) redo analysis (if there are failed masks)

% saveplace = strcat(analysis,'/threshold_1'); mkdir(saveplace); % where updated masks will be saved
% cd(wrkDir); % set directory to where TIFF images are located
close all;
RealAnalysis = 1; % 1 = eliminate cells that do not pass

% % determine COARSE or REAL analysis

    % should be 1 for coarse grain analysis step    
    FourPanelFigure = 0; % 4-panel mask image
    hist_on = 0; % make histogram to visualize threshold or not
    SaveFlag_hist = 0; % 1 = histogram figure suppressed and saved; 0 = figure pops up but not saved;
    % should be 1 for real analysis step
    SingleMaskFigure = 0; % 1 for real analysis; 2-panel mask image
    SaveFlag_mask = 0; % 1 for real analysis; 1 = mask figure suppressed and saved; 0 = figure pops up but not saved;
    
if RealAnalysis == 0;
    
    % HANDLES: should be 1 for coarse grain analysis step    
    FourPanelFigure = 1; % 4-panel mask image
    hist_on = 1; % make histogram to visualize threshold or not
    SaveFlag_hist = 1; % 1 = histogram figure suppressed and saved; 0 = figure pops up but not saved;
    
    q = 3; % MANUAL CHANGE for test
    q_range = q;
    xy_range = expt_conditions_xy{q}(3) % an example xy position corresponding to the q
    
    % determine thresholds
    thresh_percent_max_Q(q) = 0; % a relic from old code; percentage of cells to eliminate from the top
    thresh_min_a(q) = 10;
    thresh_max_a(q) = 120;
    thresh_max_Q(q) = 450; % specify specific value (instead of calculating from percent)
    
end
clearvars num_thresh_Q_max Q75_sorted

if RealAnalysis == 1;
    % HANDLES: should be 1 for real analysis step    
    SingleMaskFigure = 1; % 1 for real analysis; 2-panel mask image
    SaveFlag_mask = 1; % 1 for real analysis; 1 = mask figure suppressed and saved; 0 = figure pops up but not saved;

    q_range = 1:total_cond;
    h = waitbar(0,'out of total xy positions'); % waitbar only for real analysis
elseif RealAnalysis == 2; % redo
    % HANDLES: should be 1 for real analysis step    
    SingleMaskFigure = 1; % 1 for real analysis; 2-panel mask image
    SaveFlag_mask = 1; % 1 for real analysis; 1 = mask figure suppressed and saved; 0 = figure pops up but not saved;
    
    % try to redo a whole q and not just a few xy
    q_range = 3; % MANUAL entry
end


for q = q_range;
    
    if [RealAnalysis == 1 | RealAnalysis == 2];
        xy_range = expt_conditions_xy{q};
    end
    
    for xy = xy_range;    
        if RealAnalysis == 1; % for waitbar
            waitbar(xy / total_xy);
        end
        % skip excluded images
        if isempty(bfImgs(xy).name) == 1;
            continue            
        end
        % load image & set input variables
        phasename = bfImgs(xy).name;
        bwf_txy = bwf{xy};
        cells_txy = cells{xy};
        area_cell_temp = cell2mat(struct2cell(area{xy}));
        Q75_cell_temp = Q75{xy}; % UQ75 values of all cells in that image
        
        % set thresholds
        thr_min_a = thresh_min_a(q);
        thr_max_a = thresh_max_a(q);
        thr_max_Q = thresh_max_Q(q);
        thr_percent_max_Q = thresh_percent_max_Q(q);            
                
        % run script
        [a_elim(xy) Q_elim(xy) cells_post_thresh{xy}...
         bwf_thresh_a{xy} bwf_thresh_a_elim{xy} ...
         bwf_thresh_Q{xy} bwf_thresh_Q_elim{xy} ...
         bwf_post_thresh{xy} bwf_post_thresh_elim{xy}] = ...
            ThresholdIterator(q,phasename, bwf_txy, cells_txy, thr_min_a, thr_max_a, thr_max_Q, thr_percent_max_Q,...
            area_cell_temp, Q75_cell_temp, ...
            SaveFlag_mask, SaveFlag_hist, RealAnalysis, FourPanelFigure, SingleMaskFigure, saveplace, hist_on, Q75_cum{q}, area_cum{q});
           
    end
end
clearvars ans q xy xy_range q_range thresh_min_area percent_max_Q
clearvars phasename bwf_txy cells_txy area_cell_temp Q75_cell_temp
clearvars FourPanelFigure SingleMaskFigure SaveFlag_mask SaveFlag_hist hist_on 
clearvars thr_min_a thr_max_a thr_max_Q thr_percent_max_Q thresh_max_area


% save workplace
if RealAnalysis == 1;
   close(h)
   save(strcat(analysis,'/05_ThresholdIterator_1.mat'), '-regexp', '^(?!(bwf|bwf_5|bwf_thresh_a|bwf_thresh_a_elim|bwf_thresh_Q|bwf_thresh_Q_elim|bwf_post_thresh|bwf_post_thresh_elim)$).');
   save(strcat(analysis,'/bwf_post_thresh.mat'),'bwf_post_thresh','bwf_thresh_a','bwf_thresh_a_elim','bwf_thresh_Q','bwf_thresh_Q_elim','bwf_post_thresh_elim'); % save the final masks variable
   save(strcat(analysis,'/cells_post_thresh.mat'),'cells_post_thresh'); % save the final masks variable
   clearvars h RealAnalysis
elseif RealAnalysis == 2;
    save(strcat(analysis,'/05_ThresholdIterator_2.mat'), '-regexp', '^(?!(bwf|bwf_5|bwf_thresh_a|bwf_thresh_a_elim|bwf_thresh_Q|bwf_thresh_Q_elim|bwf_post_thresh|bwf_post_thresh_elim)$).');
   save(strcat(analysis,'/bwf_post_thresh.mat'),'bwf_post_thresh','bwf_thresh_a','bwf_thresh_a_elim','bwf_thresh_Q','bwf_thresh_Q_elim','bwf_post_thresh_elim'); % save the final masks variable
   save(strcat(analysis,'/cells_post_thresh.mat'),'cells_post_thresh'); % save the final masks variable
   clearvars RealAnalysis
end

%% 08) make thresholded mask figures for visual inspection

% visually check ensure mask formation success
close all;
SaveIm = 1; % 1 = save image

subplot_row = 4;
subplot_col = 6;
total_im_per_subplot = subplot_row * subplot_col; 
xy_figure = [1:total_im_per_subplot:total_xy];

cells_elim_color = colormap(copper(2)); % value 2 gives the highest contrast
fig_count = 0; % initialize fig_count (always start with 0)
for xy = analyzed_xy;
    
    % SAVE FIGURE after the first image
    if (ismember(xy,xy_figure(2:end)) | xy == analyzed_xy(end)) & SaveIm == 1;
       % save figure (high resolution)
       print(f(fig_count),strcat(saveplace,'/','maskoverlay_thresh_',num2str(fig_count),'.png'),'-dpng','-r300');        
    end
    
    % every __ images, create a new figure
    if ismember(xy,xy_figure) == 1;
        i = 1; % re-initialize i for subplot
        fig_count = fig_count + 1; % advance figure count   
        f(fig_count) = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        if SaveIm == 1;
           set(f(fig_count), 'Visible', 'off'); 
        end  
    end
    
    % skip excluded images
    if isempty(bfImgs(xy).name) == 1;
        continue            
    end
    
    % load phase image
    phasename = bfImgs(xy).name;
    ph_adj_temp = imread(phasename);
    ph_adj_temp_vis = imadjust(ph_adj_temp);
    BWoutline_final = imdilate(double(bwperim(bwf_post_thresh{xy})),strel('disk',1)); % define the thickness of outline

    % create subplot_row x subplot_col matrix of images
    subplot(subplot_row,subplot_col,i);   
    
    imshowpair(ph_adj_temp_vis, BWoutline_final);
         hold on;
    bwf_a_Q_elim = imshow(bwf_post_thresh_elim{xy}, cells_elim_color);
         set(bwf_a_Q_elim, 'AlphaData', 0.45); % eliminated cells look like they are glowing
    title(phasename);
             
    % tweak the positioning and size of each subplot
    original_pos = get(gca, 'Position'); % Position = [left bottom width height]
    new_pos = original_pos + [-0.05 -0.05 0.06 0.06]; % good for 4x6
    %new_pos = original_pos + [-0.04 -0.04 0.045 0.045]; % good for 5x7
    set(gca, 'Position', new_pos);
    
    % advance i for subplot position
    i = i + 1;
       
end

clearvars ans cells_elim_color f fig_count i new_pos original_pos
clearvars ph_adj_temp ph_adj_temp_vis phasename BWoutline_final
clearvars SaveIm saveplace subplot_col subplot_row xy xy_figure



%% 09) [PART 1] visualize the number of eliminated cells
close all; 
% compile by expt condition
for q = 1:total_cond;
    
    xy_temp = expt_conditions_xy{q};
    
    % first compute elimination ratio (normalized to total detected cells    
        for xy = xy_temp;
            % skip excluded images
            if isempty(bfImgs(xy).name) == 1;
                continue            
            end
            numCells = cells{xy}.NumObjects;
            a_elim_value = a_elim(xy);
            Q_elim_value = Q_elim(xy);
            a_elim_ratio(xy) = a_elim_value / numCells;
            Q_elim_ratio(xy) = Q_elim_value / numCells;
        end
  
    % second compute mean and standard deviation for each expt cond
        a_elim_ratio_cond(q) = mean(a_elim_ratio(xy_temp));
        a_elim_ratio_cond_std(q) = std(a_elim_ratio(xy_temp));
        Q_elim_ratio_cond(q) = mean(Q_elim_ratio(xy_temp));
        Q_elim_ratio_cond_std(q) = std(Q_elim_ratio(xy_temp));

end
clearvars q xy_temp xy numCells a_elim_value Q_elim_value


% plot line colors
colors = hsv(total_cond);

% plot number of eliminated cells (a and Q elim)
a_uq = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
for q = analyzed_q;
    i = find(q==analyzed_q); % get the index
    
    % AREA THRESHOLD: plot number of eliminated cells
    p(i) = errorbar(1, a_elim_ratio_cond(:,q),a_elim_ratio_cond_std(:,q),'.-',...
                    'MarkerSize',30,...
                    'LineWidth',2,...
                    'Color',colors(i,:));
        hold on;
    % UQ THRESHOLD: plot number of eliminated cells         
    p(i) = errorbar(1.5, Q_elim_ratio_cond(:,q),Q_elim_ratio_cond_std(:,q),'.-',...
                'MarkerSize',30,...
                'LineWidth',2,...
                 'Color',colors(i,:));
end
hold off; 

set(gca,'xtick',[1,1.5],'xticklabel',{'min & max area','max upper-quartile intensity'},'FontSize',14)
ylabel('number of eliminated cells normalized by total cells detected in each image','FontSize',16);
suptitle({'number of eliminated cells due to minimum & maximum cell area threshold and maximum upper-quartile intensity',...
          strcat('errorbars = standard deviation across images;         experiment: ',expt_date)});
xlim([0.5 2]);

% legend
legend(p(:),expt_conditions_string,'FontSize',16); % 1 for bars, 2 for line

%% 09) [PART 2] save figures

saveplace = strcat(analysis,'/mask plots'); % mask plots folder
saveas(a_uq,strcat(saveplace,'/cells_elim_area_UQ_threshold.png'));

clearvars i p q ans a_uq colors
close all; 


%% 10) [STEP 1] create background mask (BckgMaskGenerator.m)
% step 1: test thresholds and dilution radius
% cd(wrkDir)
close all;

% STEP 1: assess pixel histogram
q = 3; % MANUAL CHANGE
xy = expt_conditions_xy{q}(11);

% record thresholds
bckg_thresh_min(q) = 440;
bckg_thresh_max(q) = 520;
dil_radius(q) = 25; % set the radius of inflation (usually ~10 pixels)

% load images
ph_adj_temp = imread(bfImgs(xy).name); % load a dummy phase image
ph_adj_temp_vis = imadjust(ph_adj_temp);
yfp_adj_temp = imadjust(imread(YFP(xy).name));

% form test background mask
[bckg_mask_test{xy}] = BckgMaskGenerator(bfImgs(xy).name,bwf_5{xy},dil_radius(q),bckg_thresh_min(q),bckg_thresh_max(q));

% visualize histogram and background mask
figure('units', 'normalized', 'outerposition', [0 0 1 1]); 
    subplot(221); histogram(ph_adj_temp);
        hold on;
        plot(bckg_thresh_min(q),0,'r*','MarkerSize',24); 
        plot(bckg_thresh_max(q),0,'r*','MarkerSize',24);
    subplot(222); imshow(yfp_adj_temp.*uint16(bckg_mask_test{xy})); title('bckg_mask YFP (final)');
    subplot(223); imshow(ph_adj_temp_vis.*uint16(bckg_mask_test{xy})); title('bckg_mask phase (final)');
    subplot(224); imshow(ph_adj_temp_vis.*uint16(~bckg_mask_test{xy})); title('bckg_mask phase(inverted, final)');

    
%% 10) [PART 2] create background mask
% STEP 2: generate background mask (takes some time)

close all;
h = waitbar(0,'out of total xy positions');
for q = analyzed_q;
    for xy = expt_conditions_xy{q};       
        waitbar(xy / total_xy);
        
        % skip excluded images
        if isempty(bfImgs(xy).name) == 1;
            continue            
        end
        
        % specify variables
        phasename = bfImgs(xy).name;
        mask_cell = bwf_5{xy}; % cell mask, before border clearance

        % run script
        [bckg_mask{xy}] = BckgMaskGenerator(phasename, mask_cell, ...
                                dil_radius(q), bckg_thresh_min(q), bckg_thresh_max(q));
                            
    end
end
close(h)
clearvars ans h q xy t phasename mask_cell bckg_mask_test 
clearvars ph_adj_temp yfp_adj_temp


%% 10) [PART 3] test background mask visualize

xy = 33; % load a dummy image
phase = imread(bfImgs(xy).name);
% yimg = imadjust(imread(YFP(xy).name));
% timg = imadjust(imread(TFP(xy).name));
figure('units', 'normalized', 'outerposition', [0 0 1 1]); 
subplot(121); imshow(phase.*uint16(bckg_mask{xy}),[300 600]); title('bckg_mask phase (final)');
subplot(122); imshow(phase.*uint16(~bckg_mask{xy}),[300 600]); title('bckg_mask phase (inverted, final)');

clearvars ans xy phase


%% 10) [PART 4] save background masks
close all;
save(strcat(analysis,'/bckg_mask.mat'),'bckg_mask','dil_radius','bckg_thresh_min','bckg_thresh_max'); % overwrite previous version


%% 11) fluorescence calculation (FluoCellComputation.m)
% MAJOR CHANGE in November 2018 version: convert images to double from
% uint16 so that negative numbers become possible. However, this conversion
% to double doesn't change anything at the level of this
% FluoCellComputation function because I'm subtracting at the cell level,
% and not at the pixel level.

h = waitbar(0,'out of total xy positions');
for xy = 1:total_xy;
    waitbar(xy / total_xy);
       
    % skip excluded images
    if isempty(bfImgs(xy).name) == 1;
        continue            
    end
                        
    % specify variables
    rfp = double(imread(RFP(xy).name)); % convert image from uint16 to double
    yfp = double(imread(YFP(xy).name)); % convert image from uint16 to double
    tfp = double(imread(TFP(xy).name)); % convert image from uint16 to double   
    c_mask = bwf_post_thresh{xy}; % cell mask
    b_mask = bckg_mask{xy}; % background mask
    cellpixelID = cells_post_thresh{xy}.PixelIdxList; % PixelIdxList

   [rfp_bckg(xy),      yfp_bckg(xy),      tfp_bckg(xy), ...
    rfp_bckg_std(xy),  yfp_bckg_std(xy),  tfp_bckg_std(xy), ...
    rfp_fintcells{xy}, yfp_fintcells{xy}, tfp_fintcells{xy}, ...
    rfp_fintview(xy),  yfp_fintview(xy),  tfp_fintview(xy), ...
    rfp_std(xy),       yfp_std(xy),       tfp_std(xy), ...
    rfp_sem(xy),       yfp_sem(xy),       tfp_sem(xy), ...
    totalpixelnum_bckg(xy)] = FluoCellComputation(rfp, yfp, tfp, c_mask, b_mask, cellpixelID);

end
close(h)
clearvars ans h xy rfp yfp tfp c_mask b_mask cellpixelID

% save workspace
save(strcat(analysis,'/06_FluoCellComputation_all.mat'));


%% 12) extract other information - number of cells

% number of cells by xy position
for xy = 1:total_xy;       
    % skip excluded images
    if isempty(bfImgs(xy).name) == 1;
        continue            
    end
    
     cell_num(xy) = length(yfp_fintcells{xy});
end

% average number of cells by expt condition
for q = 1:total_cond;
    xy_range = expt_conditions_xy{q}; 
        cell_num_avg(q) = mean(cell_num(xy_range));
        cell_num_std(q) = std(cell_num(xy_range));
end
clearvars ans xy q xy_range

save(strcat(analysis,'/07_CellNumber_preThresh.mat'),'cell_num','cell_num_avg','cell_num_std');


%% 13) extract other information - cell size

cell_size = {}; % initialize

for xy = 1:total_xy;     
    
    % skip excluded images
    if isempty(bfImgs(xy).name) == 1;
        continue            
    end
    
    % pull out cell info of each image
    cells_image = cells_post_thresh{xy}.PixelIdxList; 
    % initialize temporary container        
    size_temp = []; 

        % extract number of pixels in each cells
        for cell = 1:length(cells_image);
            size_temp(cell) = length(cells_image{cell});
        end

    % store in permanent container
    cell_size{xy} = size_temp;

    clearvars cells_image % initialize
end
clearvars size_temp cells_image xy cell

save(strcat(analysis,'/08_CellSize_preThresh.mat'),'cell_size');

%% 15) FINAL SAVE

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
                    expt_conditions_string expt_conditions_xy expt_date ...
                    analyzed_q analyzed_xy ...
                    cell_num cell_num_avg cell_num_std cell_size
                
save(strcat(analysis,'/06_FluoCellComputation_analysisonly.mat'));


%% NEXT STEP: 03_SpectralLeakageCorrector_FPctrl.m (calculate b1)
    % THEN
% 03_SpectralLeakageCorrector_FPctrl.m (to correct FP control signal with b1)
    % THEN
% 04_runFPpositiveFunc_FPcontrol.m (threshold FP signals and calculate RT
%                                   conversion factor -- but this is an old
%                                   relic step)


