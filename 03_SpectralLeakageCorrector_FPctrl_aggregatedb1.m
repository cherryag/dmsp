%% aggregated b1 calculation for spectral leakage correction, using all FPctrl datasets (aggregated)
%    && spectral leakage correction of FPctrl data

% Version November 2018 (updated 1/1/2019)

% In this version, b1 is calculated according to aggregated data from all
% FPctrl datasets of all experiments. 

% runFPpositiveFunc_FPctrl.m comes AFTER this SpectralLeakageCorrector_FPctrl.m 
% to be consistent with the DMSPkinetics analysis (where I also do spectral 
% leakage correction before YFP thresholding).


% PART I:
% 01) specify directories
% 02) load & aggregate analysis results of FPctrl images from different experiments
% 03) calculate b1 for EACH EXPERIMENT and store in fpagg (no figures)
% 04) calculate b1 for ALL cells in ALL experiments (aggregated) (with regression plots)
% 05) SAVE b1 and rsq 

% PART II (updated 1/1/2019):
% 06) load aggregated analysis results of FPctrl & aggregated b1 for spectral leakage correction of FPctrl data
% 07) spectral leakage correction of FPctrl data
% 08) collect amount of correction of each cell by expt condition
% 09) visualize leakage correction magnitude for each color
% 10) calculate mean, std, and sem with spectral leakage corrected values
% 11) SAVE spectral leakage corrected variables 

% % NEXT: runFPpositiveFunc_FPctrl.m  [[ calculate RT conversion factor ]]



%% 01) specify directories

analysis = uigetdir; % where analysis results are stored
motherdir = '/Volumes/Latte 4TB/_DMSPAvailabilityHypothesis_Datasets'; % where all experimental folders are located


%% 02) load & aggregate analysis results of FPctrl images from different experiments

% dates of all expts to be included in the b1 calculation
expts_desired = {'2018-02-24',...
                 '2018-02-25',...
                 '2018-02-27',...
                 '2018-02-28',...
                 '2018-03-01',...
                 '2018-03-02',...
                 '2018-03-03',...
                 '2018-05-09',...
                 '2018-05-10',...
                 '2018-08-03'}; 
% initialize output container
fpagg = {}; 

% change to mother directory
cd(motherdir); 

for expt = 1:length(expts_desired);
    
    % visit each experiment FPctrl folder
    expt_folder_temp = dir(strcat(expts_desired{expt},'*')); % '*' is to indicate whatever text
    data_folder_temp = strcat(expt_folder_temp.name,'/analysis/FPctrl/');
    cd(data_folder_temp); 
    
    % load variables of interest
    load('01_ParsedFileNames.mat','expt_date','expt_conditions_xy','expt_conditions_string');
    load('06_FluoCellComputation_analysisonly.mat','rfp_fintcells', 'yfp_fintcells', 'tfp_fintcells');
    
    % assign variables of interest
    fpagg(expt).expt_date = expt_date;
    fpagg(expt).expt_conditions_xy = expt_conditions_xy;
    fpagg(expt).expt_conditions_string = expt_conditions_string;
    fpagg(expt).rfp_fintcells = rfp_fintcells;
    fpagg(expt).yfp_fintcells = yfp_fintcells;
    fpagg(expt).tfp_fintcells = tfp_fintcells;
    
    % initialize for next iteration
    expt_folder_temp = []; % initialize temp folder
    clearvars -except fpagg expt expts_desired analysis motherdir % clear the variables loaded from previous iteration
    cd('/Volumes/Latte 4TB/_DMSPAvailabilityHypothesis_Datasets'); % change to mother directory
    
end

% save workspace
clearvars expt
save(strcat(analysis,'/01_AggregatedData_AllExpts.mat')); 


%% 03) calculate b1 for EACH EXPERIMENT and store in fpagg (no figures)

colors = [{'RFP','YFP','TFP'}];

for expt = 1:length(expts_desired); % for each experimental date
    
    % define variables
    rfp_fintcells = fpagg(expt).rfp_fintcells;
    yfp_fintcells = fpagg(expt).yfp_fintcells;
    tfp_fintcells = fpagg(expt).tfp_fintcells;
    expt_conditions_xy = fpagg(expt).expt_conditions_xy;
    
    % initialize temp variables
    b1 = ones(3); % b1 = slope or regression coefficient (relationship is y = b1*x)
    rsq = ones(3); % rsq = R^2, coefficient of determination
    fpagg(expt).b1 = []; fpagg(expt).rsq = []; % initialize containers
    
    for i = 1:length(colors); % true FP
       for  j = 1:length(colors); % leak FP

        % specify FP colors of interest
        FP_color = colors{i};
        leak_color = colors{j};

        if FP_color == leak_color % don't assess the diagonals of the matrix
            continue
        end

    % DEFINE LEAK & TRUE DATA FOR b1 CALCULATION      
        % specify TRUE color (MANUAL DEFINITION of q)
        if FP_color == 'RFP'
            a = 3; % matrix position for linear regression matrix
            q = 1;
            FP_data_temp = rfp_fintcells;       
        elseif FP_color == 'YFP'
            a = 2; % matrix position
            q = 2;
            FP_data_temp = yfp_fintcells;
        elseif FP_color == 'TFP'
            a = 1; % matrix position
            q = 3;
            FP_data_temp = tfp_fintcells;       
        end

        % specify LEAKAGE color
        if leak_color == 'RFP'
           b = 3; % matrix position
           leak_data_temp = rfp_fintcells;      
        elseif leak_color == 'YFP'
           b = 2; % matrix position
           leak_data_temp = yfp_fintcells;
        elseif leak_color == 'TFP'
           b = 1; % matrix position
           leak_data_temp = tfp_fintcells;
        end

    % SIMPLE LINEAR REGRESSION
        % https://ch.mathworks.com/help/matlab/data_analysis/linear-regression.html#buva8q5
        % here, I'm assuming that the regression goes through 0,0; in the end, I want a matrix of slopes b1

        x = []; y = []; % initialize data containers on which linear regression will be performed
        for xy = expt_conditions_xy{q}; % q is defined manually according to color                  
            x = [x FP_data_temp{xy}(~isnan(FP_data_temp{xy}))]; % only include non-Nan values
            y = [y leak_data_temp{xy}(~isnan(leak_data_temp{xy}))]; % only include non-NaN values  
        end

        % perform linear regression    
        x = x'; y = y'; % must be vertical for mldivide \ to work
        format long; % displays many digits of output
        b1(b,a) = x\y;  % b1 is the slope or regression coefficient (relationship is y = b1*x)
                        % the \ operator performs a least-squares regression.
        rsq(b,a) = 1 - sum((y - b1(b,a)*x).^2)/sum((y - mean(y)).^2); % coefficient of determination (R^2)
        
       end
    end
    
    % store calculated values in permanent container
    fpagg(expt).b1 = b1;
    fpagg(expt).rsq = rsq;
    
end
clearvars a b b1 colors expt FP_color FP_data_temp i j leak_color leak_data_temp q rsq x xy y

% save workspace
save(strcat(analysis,'/02_b1_each_expt.mat')); 


%% 04) calculate b1 for ALL cells in ALL experiments (aggregated) (with regression plots)
% [updated on 6/3/2019 to generate eps images; added if linregress_on == 1 line]
    % 1) perform linear regression on each experiment separately -> plot as different colors
    % 2) perform linear regression on aggregated data -> plot as black dotted line

% cd('/Volumes/Latte 4TB/_DMSPAvailabilityHypothesis_Datasets/aggregated_b1/all expts_b1 final');  
% saveplace = uigetdir; % folder to save figures, "regression"
 
close all;
colors = [{'RFP','YFP','TFP'}];
plot_colors = repmat(0,10,3); % parula(length(expts_desired)); % how to differentiate each experiment

hist_on = 0; % 1 = pop up and save histogram with threshold; 0 = no histogram generated
SaveFig = 0; % 1 = save figures
linregress_on = 0; % 1 = do linear regression on scatter plot and plot individual experiment
common_axis_on = 1; % 1 = set common axis
common_axis = [0 10000 -20 180]; % set common axis to be able to compare plots
true_fp_thresh = 0; % thresholding on 0 because otherwise, too many negative numbers and the b1 calculation does not proceed properly
                    % not a real threshold
                    
                   
% initialize data aggregation containers
if linregress_on == 1
    b1_agg = ones(3); % initialize ; b1 = slope or regression coefficient (relationship is y = b1*x)
    rsq_agg = ones(3); % initialize ; rsq = R^2, coefficient of determination
end

for i = 1:length(colors); % true FP
   for  j = 1:length(colors); % leak FP
       
         % specify FP colors of interest
            FP_color = colors{i};
            leak_color = colors{j};

            if FP_color == leak_color % don't assess the diagonals of the matrix
                continue
            end
                                       
            
       % POPUP A FIGURE for each color combination
            FluoLeakScatter = figure('units','normalized','outerposition',[0 0 1 1]); 
       
       % INITIALIZE data aggregation containers
            x_agg = []; y_agg = []; index_agg = []; b1_string_for_legend = {}; p = []; % aggregated for each color combination
                              
       % ITERATE through each experiment for this color combination
       for expt = 1:length(expts_desired); 
                
            % DEFINE DATA TO INCLUDE AND PLOT        
                % specify expt condition for each FP color control (MANUAL CHANGE of q)
                if FP_color == 'RFP'
                    a = 3; % matrix position for linear regression matrix
                    q = 1;
                    FP_data_temp = fpagg(expt).rfp_fintcells; 
                elseif FP_color == 'YFP'
                    a = 2; % matrix position
                    q = 2;
                    FP_data_temp = fpagg(expt).yfp_fintcells;
                elseif FP_color == 'TFP'
                    a = 1; % matrix position
                    q = 3;
                    FP_data_temp = fpagg(expt).tfp_fintcells;   
                end

                % specify the leakage channel data to plot
                if leak_color == 'RFP'
                   b = 3; % matrix position
                   leak_data_temp = fpagg(expt).rfp_fintcells; 
                elseif leak_color == 'YFP'
                   b = 2; % matrix position
                   leak_data_temp = fpagg(expt).yfp_fintcells;
                elseif leak_color == 'TFP'
                   b = 1; % matrix position
                   leak_data_temp = fpagg(expt).tfp_fintcells;
                end     
                  
            % AGGREGATE all data points from each experiment
                x_expt = []; y_expt = []; index_expt = []; % initialize container for each experiment
                for xy = fpagg(expt).expt_conditions_xy{q}; % for all xy positions (i.e. images) of all colors) 
                   % aggregate for each experiment (just for plot)
                   x_expt = [x_expt FP_data_temp{xy}];
                   y_expt = [y_expt leak_data_temp{xy}];
                   % aggregate across all experiments (for aggregated b1 calculation)
                   x_agg = [x_agg FP_data_temp{xy}];
                   y_agg = [y_agg leak_data_temp{xy}];
                end

            % THRESHOLD at 0 on true FP color (histogram below)
                index_expt = find(x_expt > true_fp_thresh);
                    x_expt = x_expt(index_expt); y_expt = y_expt(index_expt);
                    num_cells_color_expt(i,expt) = length(x_expt);
                index_agg = find(x_agg > true_fp_thresh);
                    x_agg = x_agg(index_agg); y_agg = y_agg(index_agg);
                    num_cells_color_agg(i) = length(x_agg)
         
            % SCATTER PLOT figure (additive after each expt; each expt is different color)    
                % scatter(x_expt(1:10:end), y_expt(1:10:end)... if every 10 cells
                p(expt) = scatter(x_expt, y_expt, 15, 'filled',...
                                                      'MarkerFaceColor', plot_colors(expt,:), ...
                                                      'MarkerFaceAlpha', 1); 
                hold on;
                     
            % LINEAR REGRESSION PLOT (for each expt)
            if linregress_on == 1
                rsq_expt_plot = []; b1_expt_plot = [];
                    b1_expt_plot = fpagg(expt).b1; % define b1 (just for plot)
                    rsq_expt_plot = fpagg(expt).rsq; % define rsq (just for plot
                plot(x_expt, b1_expt_plot(b,a)*x_expt, 'LineWidth',1.5,'Color',plot_colors(expt,:),'HandleVisibility','off'); % HandleVisibility is for legend
                    % build the text for the legend
                    b1_string_for_legend{expt} = num2str(b1_expt_plot(b,a)); 
            end

       end % END of experiment iteration
            
       % LEGEND
            if linregress_on == 1 
            leg = legend(p(1:expt), strcat(expts_desired(1:expt),';  b1 = ',cellfun(@num2str,b1_string_for_legend,'un',0)), 'Location','NorthWest');         
            leg.FontSize = 14;
            set(gca,'FontSize',24); % tic size
            end
            
       % LINEAR REGRESSION on aggregated data (b1 generation for each color combination)
            % https://ch.mathworks.com/help/matlab/data_analysis/linear-regression.html#buva8q5
            % here, I'm assuming that the regression goes through 0,0; in the end, I want a matrix of slopes b1.
            if linregress_on == 1
                x_agg = x_agg'; y_agg = y_agg'; % must be vertical for mldivide \ to work
                format long; % displays many digits of output
                b1_agg(b,a) = x_agg(1:1000:end)\y_agg(1:1000:end);  % b1 is the slope or regression coefficient (relationship is y = b1*x)
                b1_agg(b,a) = x_agg\y_agg;  % b1 is the slope or regression coefficient (relationship is y = b1*x)
                    % the \ operator performs a least-squares regression.
                rsq_agg(b,a) = 1 - sum((y_agg - b1_agg(b,a)*x_agg).^2)/sum((y_agg - mean(y_agg)).^2); % coefficient of determination (R^2)
            end

%        % PLOT aggregated linear regression            
%             agg_txt = strcat('linear regression (slope=',num2str(b1_agg(b,a)),'; R2=',num2str(rsq_agg(b,a)),')');
%             x_interval = (max(x_agg)-min(x_agg)) / 500; % want 500 dots for each plot
%             x_range = min(x_agg):x_interval:max(x_agg); % define range of x for line
%             y_range = b1_agg(b,a) * x_range;
%             plot(x_range, y_range,':r', 'LineWidth',3); % plot the line
%                        
%       % PLOT PROPERTIES                            
%             xlabel(strcat(FP_color,' (true FP channel)'),'FontSize',24);
%             ylabel(strcat(leak_color,' (leak channel)'),'FontSize',24);
%             title(strcat('aggregated spectral leakage assessment;    black dotted line = ', agg_txt),'FontSize',14);          
            
       % set common axis (if desired)         
            if common_axis_on == 1;
                axis(common_axis); % common axis for all figures
            end

       % save figure
            if SaveFig == 1;
                if common_axis_on == 0;
                    print(FluoLeakScatter,strcat(saveplace,'/aggregated_threshat0_fp',FP_color,'_leakin',leak_color),'-dpng');
                    close all
                elseif common_axis_on == 1;
%                    print(FluoLeakScatter,strcat(saveplace,'/aggregated_threshat0_fp',FP_color,'_leakin',leak_color,'_sameaxis'),'-dpng');
%                    print(FluoLeakScatter,strcat(saveplace,'/aggregated_every10cell_threshat0_fp',FP_color,'_leakin',leak_color,'_sameaxis'),'-deps');
                    print(FluoLeakScatter,strcat(saveplace,'/aggregated_threshat0_fp',FP_color,'_leakin',leak_color,'_sameaxis'),'-deps');
                    close all
                end
            end
                 
       % quick histogram to sanity check the threshold at 0
       if hist_on == 1;
            f = figure; h = histogram(x_agg,500); hold on; plot(true_fp_thresh,0:10:max(h.Values),'r.','MarkerSize',10);
            xlabel(strcat(FP_color,' (true FP channel)'),'FontSize',24); ylabel('counts'); set(gca, 'fontsize',18);   
                print(f,strcat(saveplace,'/hist_thresh_',num2str(true_fp_thresh),'_fp',FP_color,'_leakin',leak_color),'-dpng'); close all;
       end
       
                
   end % END of leak color
end % END of true FP color

% number of cells included in aggregated b1 regression (every 1000)
num_cells_color_agg

% number of cells from every experiment
num_cells_color_expt
sum(num_cells_color_expt(1,:)) == num_cells_color_agg(1) % RFP should be same length
sum(num_cells_color_expt(2,:)) == num_cells_color_agg(2) % YFP
sum(num_cells_color_expt(3,:)) == num_cells_color_agg(3) % TFP

% mean and std number of cells
mean(num_cells_color_expt(1,:))
std(num_cells_color_expt(1,:))

mean(num_cells_color_expt(2,:))
std(num_cells_color_expt(2,:))

mean(num_cells_color_expt(3,:))
std(num_cells_color_expt(3,:))


% output b1
b1_agg


%% 05) SAVE b1 and rsq 

close all;
clearvars a b b1 colors expt FP_color FP_data_temp i j leak_color leak_data_temp q rsq x xy y ans f h
clearvars leg p SaveFig common_axis_on FluoLeakScatter hist_on index_expt index_agg
clearvars x_agg y_agg x_expt y_expt x_interval x_range y_range agg_text rsq_expt_plot b1_expt_plot

% save workspace
save(strcat(analysis,'/02_b1_each_expt.mat')); 

% save linear regression data
save(strcat(analysis,'/RegressionSlope_b1_aggregated.mat'),'b1_agg'); % overwrite previous version
save(strcat(analysis,'/Rsquared_rsq_aggregated.mat'),'rsq_agg'); % overwrite previous version






%% PART II (updated 1/1/2019 for RT conversion calculation)
% 06) load aggregated analysis results of FPctrl & aggregated b1 for spectral leakage correction of FPctrl data

% load aggregated FPctrl data (fintcells)
load('/Volumes/Latte 4TB/_DMSPAvailabilityHypothesis_Datasets/aggregated_b1/all expts_b1 final/01_AggregatedData_AllExpts.mat')

% load aggregated b1
load('/Volumes/Latte 4TB/_DMSPAvailabilityHypothesis_Datasets/aggregated_b1/all expts_b1 final/RegressionSlope_b1_aggregated.mat')

% directory
cd('/Volumes/Latte 4TB/_DMSPAvailabilityHypothesis_Datasets/aggregated_b1/all expts_b1 final')


%% 07) spectral leakage correction of FPctrl data
% using aggregated b1
% measured signal = (matrix of slopes) x (actual signal)
% Y = b1*X  --> solve for X
% X = Y*b1 or X = b1\Y

total_time = 1; % for FPctrl, no time range

for expt = 1:length(expts_desired)
    
    % initialize before processing each experiment
    allc_fintcells = {}; leakcorrect_fintcells = {};
    tfp_leakcorr_fintcells = {}; yfp_leakcorr_fintcells = {}; rfp_leakcorr_fintcells = {};
    tfp_leakcorrdiff_fintcells = {}; yfp_leakcorrdiff_fintcells = {}; rfp_leakcorrdiff_fintcells = {};
    
    clearvars rfp_fintcells yfp_fintcells tfp_fintcells
    
    % print progress
    strcat('the experiment being processed is : ', fpagg(expt).expt_date, '(expt number : ',num2str(expt),')')
    
    % define datasets for experiment
    total_xy = length(fpagg(expt).rfp_fintcells);
    rfp_fintcells = fpagg(expt).rfp_fintcells;
    yfp_fintcells = fpagg(expt).yfp_fintcells;
    tfp_fintcells = fpagg(expt).tfp_fintcells;
    
    for t = 1:total_time;
        for xy = 1:total_xy;
            c = length(tfp_fintcells{t,xy});
            
            for cell = 1:c;

                % 3x1 vector containing the fluo value of a cell in each color
                allc_fintcells{t,xy}{1,cell}(1,1) = tfp_fintcells{t,xy}(1,cell);
                allc_fintcells{t,xy}{1,cell}(2,1) = yfp_fintcells{t,xy}(1,cell);
                allc_fintcells{t,xy}{1,cell}(3,1) = rfp_fintcells{t,xy}(1,cell);

                % solve the matrix
                leakcorrect_fintcells{t,xy}{1,cell} = b1_agg \ allc_fintcells{t,xy}{1,cell}; 

                % separate
                tfp_leakcorr_fintcells{t,xy}(1,cell) = leakcorrect_fintcells{t,xy}{1,cell}(1,1);
                yfp_leakcorr_fintcells{t,xy}(1,cell) = leakcorrect_fintcells{t,xy}{1,cell}(2,1);
                rfp_leakcorr_fintcells{t,xy}(1,cell) = leakcorrect_fintcells{t,xy}{1,cell}(3,1);

                % difference between pre-correction and post-correction
                tfp_leakcorrdiff_fintcells{t,xy}{1,cell} = leakcorrect_fintcells{t,xy}{1,cell}(1,1) - allc_fintcells{t,xy}{1,cell}(1,1);
                yfp_leakcorrdiff_fintcells{t,xy}{1,cell} = leakcorrect_fintcells{t,xy}{1,cell}(2,1) - allc_fintcells{t,xy}{1,cell}(2,1);
                rfp_leakcorrdiff_fintcells{t,xy}{1,cell} = leakcorrect_fintcells{t,xy}{1,cell}(3,1) - allc_fintcells{t,xy}{1,cell}(3,1);

            end % end of cell
            
        end % end of xy
    end % end of time
    
    % store in permanent container
    fpagg(expt).tfp_leakcorr_fintcells = tfp_leakcorr_fintcells;
    fpagg(expt).yfp_leakcorr_fintcells = yfp_leakcorr_fintcells;
    fpagg(expt).rfp_leakcorr_fintcells = rfp_leakcorr_fintcells;
    fpagg(expt).tfp_leakcorrdiff_fintcells = tfp_leakcorrdiff_fintcells;  
    fpagg(expt).yfp_leakcorrdiff_fintcells = yfp_leakcorrdiff_fintcells; 
    fpagg(expt).rfp_leakcorrdiff_fintcells = rfp_leakcorrdiff_fintcells;
    
end % end of expt

clearvars ans c t xy total_xy cell expt allc_fintcells leakcorrect_fintcells tfp_leakcorr_fintcells yfp_leakcorr_fintcells
clearvars rfp_leakcorr_fintcells tfp_leakcorrdiff_fintcells yfp_leakcorrdiff_fintcells rfp_leakcorrdiff_fintcells
clearvars rfp_fintcells yfp_fintcells tfp_fintcells
    

%% 08) collect amount of correction of each cell by expt condition

for expt = 1:length(expts_desired)
    
    % define data
    tfp_leakcorrdiff_fintcells = fpagg(expt).tfp_leakcorrdiff_fintcells;
    yfp_leakcorrdiff_fintcells = fpagg(expt).yfp_leakcorrdiff_fintcells;
    rfp_leakcorrdiff_fintcells = fpagg(expt).rfp_leakcorrdiff_fintcells;
    
    for q = 1:length(fpagg(expt).expt_conditions_xy) % for each FP color
        cond = fpagg(expt).expt_conditions_xy{q};
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

            end % end of t
        end % end of xy

        tfp_leakcorrdiff_fintcells_q{q} = tfp_temp; 
        yfp_leakcorrdiff_fintcells_q{q} = yfp_temp;
        rfp_leakcorrdiff_fintcells_q{q} = rfp_temp;

        clearvars t_temp y_temp r_temp 
        clearvars tfp_temp yfp_temp rfp_temp 
        
    end % end of q (FP color)
    
    % store in fpagg
    fpagg(expt).tfp_leakcorrdiff_fintcells_q = tfp_leakcorrdiff_fintcells_q;
    fpagg(expt).yfp_leakcorrdiff_fintcells_q = yfp_leakcorrdiff_fintcells_q;    
    fpagg(expt).rfp_leakcorrdiff_fintcells_q = rfp_leakcorrdiff_fintcells_q;
    
end % end of expt


%% 09) visualize leakage correction magnitude for each color
% visualize the amount of correction with histogram, separated by color and
% for all experimental conditions. Save histogram figures.

% saveplace = strcat(pwd,'/spectral leakage correction plots'); mkdir(saveplace); 
close all;
cmap = hsv(3); % for 3 FP colors
SaveFig = 1; % 1 = save figure


for expt = 1:length(expts_desired)
    
    % define data
    rfp_leakcorrdiff_fintcells_q = fpagg(expt).rfp_leakcorrdiff_fintcells_q;
    yfp_leakcorrdiff_fintcells_q = fpagg(expt).yfp_leakcorrdiff_fintcells_q;
    tfp_leakcorrdiff_fintcells_q = fpagg(expt).tfp_leakcorrdiff_fintcells_q;
 
    for FPcolors = 1:3;
        if FPcolors == 1;
            data = rfp_leakcorrdiff_fintcells_q;
          %  xaxislim = [-1 0.5];
            color = 'RFP';
        elseif FPcolors == 2;
            data = yfp_leakcorrdiff_fintcells_q;
          %  xaxislim = [-0.8 0.1];
            color = 'YFP';
        elseif FPcolors == 3;
            data = tfp_leakcorrdiff_fintcells_q;
         %   xaxislim = [-80 5]; % should be biggest correction
            color = 'TFP';
        end
        
        figure('units','normalized','outerposition',[0 0 1 1]); 

        for q = 1:3; % for each FP color, different histogram line

            histdata_q(q) = histogram(data{q},'FaceColor',cmap(q,:)); hold on;
                centers = histdata_q(q).BinEdges + histdata_q(q).BinWidth/2; % center of each bar in histogram
                heights = [histdata_q(q).Values,0];        
                    hold on
                line_q(q) = plot(centers,heights,'LineWidth',3,'Color',cmap(q,:));    

        end
        
        title(strcat('magnitude of spectral leakage correction (in ',color,' channel); final value - original ; ', fpagg(expt).expt_date),'FontSize',20)
        xlabel('magnitude of spectral leakage correction (final value minus original)','FontSize',20);
        ylabel('frequency','FontSize',20)
        set(histdata_q(:),'facealpha',0,'edgecolor','none'); % turn OFF bars
        legend(line_q(:),{fpagg(expt).expt_conditions_string{:}},'FontSize',16,'Location','NorthWest'); 
       % xlim(xaxislim);
       
        % save figure
        if SaveFig == 1;
            saveas(gcf,strcat(saveplace,'/',color,'_LeakCorrectionDeltaHistogram_expt_',num2str(expt),'.png'));
            close all;
        end
    end

end % end of expt

clearvars cmap FPcolors SaveFig data xaxislim color q 
clearvars ans histdata_q centers heights line_q
clearvars expt cond t xy total_time
clearvars rfp_leakcorrdiff_fintcells_q yfp_leakcorrdiff_fintcells_q tfp_leakcorrdiff_fintcells_q
 

%% 10) calculate mean, std, and sem with spectral leakage corrected values

for expt = 1:length(expts_desired)
    
    % define data for each expt
    rfp_leakcorr_fintcells = fpagg(expt).rfp_leakcorr_fintcells;
    yfp_leakcorr_fintcells = fpagg(expt).yfp_leakcorr_fintcells;
    tfp_leakcorr_fintcells = fpagg(expt).tfp_leakcorr_fintcells;
    total_xy = length(rfp_leakcorr_fintcells);
    
    for xy = 1:total_xy;
        for t = 1; % for FPctrl, no time component
            
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

    % store in permanent container
    fpagg(expt).rfp_leakcorr_fintview = rfp_leakcorr_fintview;
    fpagg(expt).rfp_leakcorr_sem = rfp_leakcorr_sem;
    fpagg(expt).rfp_leakcorr_std = rfp_leakcorr_std;
    fpagg(expt).yfp_leakcorr_fintview = yfp_leakcorr_fintview;
    fpagg(expt).yfp_leakcorr_sem = yfp_leakcorr_sem;
    fpagg(expt).yfp_leakcorr_std = yfp_leakcorr_std;
    fpagg(expt).tfp_leakcorr_fintview = tfp_leakcorr_fintview;
    fpagg(expt).tfp_leakcorr_sem = tfp_leakcorr_sem;
    fpagg(expt).tfp_leakcorr_std = tfp_leakcorr_std;    
end % end of expt

% clear variables before save 
clearvars -except motherdir b1_agg expts_desired fpagg saveplace analysis


%% 11) SAVE spectral leakage corrected variables 

% cleanup fpagg
fpagg_leakcorr_details = fpagg; % contains leakcorrdiff variables
fpagg = rmfield(fpagg, {'tfp_leakcorrdiff_fintcells','rfp_leakcorrdiff_fintcells','yfp_leakcorrdiff_fintcells'}) % remove some fields

save('03_AggregatedData_SpectralLeakageCorrected.mat'); 



%% NEXT: runFPpositiveFunc_FPctrl.m
