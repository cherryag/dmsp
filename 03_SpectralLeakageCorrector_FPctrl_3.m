%% Spectral leakage corrector for PA-FP control cell data
% Version July 2018

% In this version, order of spectral leakage and FP thresholding was
% changed. With this change, runFPpositiveFunc_FPctrl.m comes AFTER this 
% SpectralLeakageCorrector_FPctrl.m to be consistent with the DMSPkinetics analysis
% (where I also do spectral leakage correction before YFP thresholding).


% PART I:
% 1) load analysis results of FPctrl & load plot styles
% 2) LINEAR REGRESSION: fluorescence true vs. leak channel scatter plots 
% 3) save b1 and rsq 

% PART II: 
% 4) load analysis results of FPctrl & b1 for spectral leakage correction
% 5) spectral leakage correction
% 6) visualize leakage correction magnitude for each color
% 7) calculate mean, std, and sem with spectral leakage corrected values
% 8) SAVE spectral leakage corrected variables (09_SpectralLeakageCorrector_FPctrl.mat)
% % NEXT: runFPpositiveFunc_FPctrl.m



%% 1) load analysis results of FPctrl & load plot styles
analysis_fpctrl = uigetdir; cd(analysis_fpctrl); % choose folder with analysis results data of FPctrl
load(strcat(analysis_fpctrl,'/06_FluoCellComputation_analysisonly.mat')); % of FPctrl
saveplace = strcat(analysis_fpctrl,'/leakage figures'); mkdir(saveplace);

% specify colors (RGB) and plot styles
% http://shirt-ediss.me/matlab-octave-more-colours/
% load H:\PlotAppearance_Variables % for PC
load('/Volumes/Latte 4TB/PlotAppearance_Variables.mat') % for Mac

%% 2) LINEAR REGRESSION: fluorescence true vs. leak channel scatter plots 
% Change in v3 of this code: feed raw cell fluorescence data, not
% thresholded for positive fluorescence yet (changed order of the
% analysis such that the fluo signal threshold comes after spectral leakage correction).

% same code as in SpectralLeakageCorrector.m (for DMSPkinetics) 
    % -- if b1 is already generated, simply load the variable.
% outputs 6 single images for summary PowerPoint
% output matrix (of regression values only), for PowerPoint:
    % columns = a (FP); rows = b (leak)
    
close all;
colors = [{'RFP','YFP','TFP'}];
b1 = ones(3); % initialize ; b1 = slope or regression coefficient (relationship is y = b1*x)
rsq = ones(3); % initialize ; rsq = R^2, coefficient of determination

SaveFig = 0; % 1 = save figures
common_axis_on = 0; % 1 = set common axis
common_axis = [-20 10000 -20 160]; % set common axis to be able to compare plots

for i = 1:length(colors); % true FP
   for  j = 1:length(colors); % leak FP
   
    % specify FP colors of interest
    FP_color = colors{i};
    leak_color = colors{j};
    
    if FP_color == leak_color % don't assess the diagonals of the matrix
        continue
    end
    
    % specify plot style
    marker = '.';
    linewidth = 1.5;
    markersize = 10;
    ManualColor = 0; % 0 = set color automatically; 1 = set color manually
        if ManualColor == 1;
            color_in = color_tfp_in; 
            color_out = color_rfp_out; % when ManualColor setting desired
        end

% DEFINE DATA TO PLOT        
    % specify expt condition for each FP color control (MANUAL CHANGE of q)
    if FP_color == 'RFP'
        a = 3; % matrix position for linear regression matrix
        q = 1;
        FP_data_temp = rfp_fintcells; 
            % plot style
            if ManualColor == 0;
               color_in = color_rfp_in;
               color_out = color_rfp_out;
            end            
    elseif FP_color == 'YFP'
        a = 2; % matrix position
        q = 2;
        FP_data_temp = yfp_fintcells;
            % plot style
            if ManualColor == 0;
               color_in = color_yfp_in;
               color_out = color_yfp_out;
            end
    elseif FP_color == 'TFP'
        a = 1; % matrix position
        q = 3;
        FP_data_temp = tfp_fintcells;       
           % plot style
           if ManualColor == 0;
              color_in = color_tfp_in;
              color_out = color_tfp_out;
           end
    end

    % specify the leakage channel data to plot
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
    % here, I'm assuming that the regression goes through 0,0
    % in the end, I want a matrix of slopes b1

    x = []; y = []; % initialize
    for xy = expt_conditions_xy{q};                  
        x = [x FP_data_temp{xy}(~isnan(FP_data_temp{xy}))]; % only include non-Nan values
        y = [y leak_data_temp{xy}(~isnan(leak_data_temp{xy}))]; % only include non-NaN values                        
    end

    % perform linear regression    
    x = x'; y = y'; % must be vertical for mldivide \ to work
    format long; % displays many digits of output
    b1(b,a) = x\y;  % b1 is the slope or regression coefficient (relationship is y = b1*x)
                    % the \ operator performs a least-squares regression.
    rsq(b,a) = 1 - sum((y - b1(b,a)*x).^2)/sum((y - mean(y)).^2); % coefficient of determination (R^2)

    
% SCATTER PLOT figure
FluoLeakScatter = figure('units','normalized','outerposition',[0 0 1 1]); 
                   
for xy = expt_conditions_xy{q};  
    hold on;
    p = plot(FP_data_temp{xy}, leak_data_temp{xy},...
                 marker, 'LineWidth', linewidth, ...
                 'Color', color_out, ...
                 'MarkerSize', 10);              
end

xlabel(strcat(FP_color,' (true FP channel)'),'FontSize',24);
ylabel(strcat(leak_color,' (leak channel)'),'FontSize',24);
title('spectral leakage assessment','FontSize',20);
set(gca,'FontSize',24);
leg = legend(p,expt_conditions_string{q});
leg.FontSize = 16;

% plot linear regression
plot(x, b1(b,a)*x, 'LineWidth',1.5,'Color','k',...
                   'DisplayName',strcat('linear regression (slope=',num2str(b1(b,a)),'; R2=',num2str(rsq(b,a)),')'));
hold off

if common_axis_on == 1;
    axis(common_axis); % common axis for all figures
end

% save figure
if SaveFig == 1;
%    print(FluoLeakScatter,strcat(saveplace,'\fp',FP_color,'_leakin',leak_color),'-depsc');
    print(FluoLeakScatter,strcat(saveplace,'\fp',FP_color,'_leakin',leak_color),'-dpng');
end

   end
end


%% 3) SAVE b1 and rsq 
close all;

% save linear regression data
save(strcat(analysis_fpctrl,'\RegressionSlope_b1_final.mat'),'b1'); % overwrite previous version
save(strcat(analysis_fpctrl,'\Rsquared_rsq_final.mat'),'rsq'); % overwrite previous version

clearvars FP_data_temp leak_data_temp q p i xy expt_cond leg
clearvars FP_color leak_color FluoLeakScatter
clearvars x y a b j ans





%% PART II
% 4) load analysis results of FPctrl & b1 for spectral leakage correction

% load 06_FluoCellComputation_analysisonly of FPctrl and b1 and reg analysis
clearvars -except b1 saveplace analysis_fpctrl
load(strcat(analysis_fpctrl,'\06_FluoCellComputation_analysisonly.mat')); % load from FPctrl folder



%% 5) spectral leakage correction

% measured signal = (matrix of slopes) x (actual signal)
% Y = b1*X  --> solve for X
% X = Y*b1 or X = b1\Y
total_time = 1; % for FPctrl, no time range
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


%% 6) visualize leakage correction magnitude for each color
% visualize the amount of correction with histogram, separated by color and
% for all experimental conditions
close all;
cmap = hsv(total_cond);
SaveFig = 0; % 1 = save figure
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
   % xlim(xaxislim);
    % save figure
    if SaveFig == 1;
        saveas(fig_hist(FPcolors),strcat(saveplace,'\',color,'_LeakCorrectionDeltaHistogram.png'));
        close all;
    end
end

clearvars cmap FPcolors SaveFig data xaxislim color fig_hist q 
clearvars ans histdata_q centers heights line_q
 

%% 7) calculate mean, std, and sem with spectral leakage corrected values

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

% clear variables before save 
clearvars ans c cell f final_timepoint h heights LeakCorr q 
clearvars start_timepoint t xy 


%% 8) SAVE spectral leakage corrected variables 

save(strcat(analysis_fpctrl,'\09_SpectralLeakageCorrector_FPctrl.mat'))



%% NEXT: runFPpositiveFunc_FPctrl.m

