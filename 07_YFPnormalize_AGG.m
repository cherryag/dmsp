%% Normalize each cell fluorescence by YFP
% version July 2018

% All analyses should involve YFP normalization.
% Run this code after runYFPpositiveFunc_4.m.

% This allows comparison between cells and leads to a more fair averaging process.

% 0) load data files
% 1) collect all RFP and TFP values of all cells
% 2) scatter plot RFP vs. YFP and TFP vs. YFP
% 3) divide RFP and TFP value of each cell by its YFP value
% 4) collect all YFP normalized values of all cells by expt condition
% 5) average fluo value of each image (each xy position) after YFP normalization
% 6) average fluo value of each expt condition (each q) after YFP normalization
% 7) collect all fluorescence intensities AFTER YFP normalized
% 8) visualize RFP vs. YFP and TFP vs. YFP of each cell AFTER YFP normalization 
% 9) SAVE variables



%% 0) load data files

% load relevant data
clearvars -except analysis_kinetics
thresh_YFP = 50; % MANUAL ENTRY; indicate which YFP threshold data to load
cd(analysis_kinetics);
load(strcat('11_YFPon_thresh_',num2str(thresh_YFP),'.mat'));
load 01_ParsedFileNames.mat

% specify colors (RGB) and plot styles
% http://shirt-ediss.me/matlab-octave-more-colours/
load H:\PlotAppearance_Variables

%% 1) collect all fluorescence intensities and do a scatter plot
% Collect all fluo intensities of all cells, spectral leakage-corrected AND YFP threshold. 

inputdata_r = rfp_fintcells_ypos;
inputdata_y = yfp_fintcells_ypos;
inputdata_t = tfp_fintcells_ypos;
rfp_allcells_leakcorr = [];
yfp_allcells_leakcorr = [];
tfp_allcells_leakcorr = [];
for xy = analyzed_xy;
    for t = 1:total_time;       
        rfp_allcells_leakcorr = [rfp_allcells_leakcorr inputdata_r{t,xy}];
        yfp_allcells_leakcorr = [yfp_allcells_leakcorr inputdata_y{t,xy}];
        tfp_allcells_leakcorr = [tfp_allcells_leakcorr inputdata_t{t,xy}];  
    end
end
clearvars inputdata_r inputdata_t inputdata_y
savename = strcat('\12_FPallcells_afterLeakCorr_Ythresh_',num2str(thresh_YFP),'.mat')
save(strcat(analysis_kinetics,savename),'rfp_allcells_leakcorr','tfp_allcells_leakcorr','yfp_allcells_leakcorr'); % save variable


%% 2) visualize RFP vs. YFP and TFP vs. YFP of each cell before normalization
% OPTIONAL.
% not very useful because not separated by q.

scatterPlot = figure('units','normalized','outerposition',[0 0 1 1]);

x_data = yfp_allcells_leakcorr;    x_text = 'YFP';
y_data = rfp_allcells_leakcorr;    y_text = 'RFP';

color = color_rfp_out;
markersize = 4;
title_text = {'Scatter Plot: each dot = one cell; all cells of all time',expt_date};
fig_savename = '\YFPvsRFP_scatter_leakcorr.png';

leg_text = expt_date;

scatter(x_data,y_data,markersize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.1);


% plot label
suptitle(title_text);
xlabel(x_text,'FontSize',20)
ylabel(y_text,'FontSize',20);

% legend
leg = leg_text; % legend
legend({leg},'FontSize',20);


%% 2) visualize RFP vs. YFP and TFP vs. YFP of each cell before normalization 
% OPTIONAL.
% scatter plot with different color dots for time. 

time_range = [1 5 10 30];
fitPerform = 0;
x_data = rfp_fintcells_q_ypos;    x_text = 'RFP';
y_data = tfp_fintcells_q_ypos;    y_text = 'TFP';
color = hsv(4);
markersize = 4;
leg_text = {'t1','t5','t10','t30'};
q_range = 2;

[scatterPlot] = scatter_2D(q_range, time_range, x_data, y_data, color, markersize, fitPerform,...
                           title_text, x_text, y_text, leg_text);

% axis range
% axis([-20 2500 -20 120]);

% print(scatterPlot,strcat(saveplace,fig_savename),'-dpng','-r300'); % resolution 300


%% 3) divide RFP and TFP value of each cell by its YFP value
close all

% initialize
rfp_fintcells_leakcorr_ypos_ynorm = {};
tfp_fintcells_leakcorr_ypos_ynorm = {};

for xy = 1:total_xy;
    for t = 1:total_time;
        
    % exit if the image is evaporated
%     if ismember([t,xy],evaporated_txy,'rows');
%         strcat('[t, xy] of [',num2str(t),', ',num2str(xy),'] is an evaporated image and is skipped')
%         continue % skip to the next cycle 
%     end
        
    yfp_cells_temp = yfp_fintcells_ypos{t,xy};
    rfp_cells_temp = rfp_fintcells_ypos{t,xy};
    tfp_cells_temp = tfp_fintcells_ypos{t,xy};
    
    % initialize before each image
    rfp_ynorm_temp = []; tfp_ynorm_temp = [];
    
    for i = 1:length(yfp_cells_temp); % refer to each cell in the image   
        
        % normalize each cell by its YFP value
        if isnan(yfp_cells_temp(i)) == 0; % if the cell fluo value is not NaN
           rfp_ynorm_temp(i) = rfp_cells_temp(i)./yfp_cells_temp(i);
           tfp_ynorm_temp(i) = tfp_cells_temp(i)./yfp_cells_temp(i);
        else
            continue % assigned 0 if NaN
        end      
        
    end
    
    rfp_fintcells_leakcorr_ypos_ynorm{t,xy} = rfp_ynorm_temp;
    tfp_fintcells_leakcorr_ypos_ynorm{t,xy} = tfp_ynorm_temp;        
        
    end
end
clearvars yfp_cells_temp rfp_cells_temp tfp_cells_temp
clearvars t xy i rfp_ynorm_temp tfp_ynorm_temp


% 4) collect all YFP normalized values of all cells by expt condition

for t = 1:total_time;
    for q = q_range_analyzed;
        r_data_cum = []; t_data_cum = []; % initialize
        
        for xy = expt_conditions_xy{q};
            r_data_temp = []; t_data_temp = []; % initialize
            
            % collect all cells by expt conditions
            r_data_temp = rfp_fintcells_leakcorr_ypos_ynorm{t,xy};
            t_data_temp = tfp_fintcells_leakcorr_ypos_ynorm{t,xy};

            % append
            r_data_cum = [r_data_cum r_data_temp];
            t_data_cum = [t_data_cum t_data_temp];
        
        end
        
        rfp_fintcells_q_leakcorr_ypos_ynorm{t,q} = r_data_cum;
        tfp_fintcells_q_leakcorr_ypos_ynorm{t,q} = t_data_cum;
               
    end
end
clearvars r_data_temp t_data_temp r_data_cum t_data_cum t q xy


% 5) average fluo value of each image (each xy position) after YFP normalization

for xy = 1:total_xy;
    for t = 1:total_time;
        
        % exit if the image is evaporated
%         if ismember([t,xy],evaporated_txy,'rows');
%             strcat('[t, xy] of [',num2str(t),', ',num2str(xy),'] is an evaporated image and is skipped')
%             continue % skip to the next cycle 
%         end        

        % eliminate 0-value cells (i.e. cells that are YFP negative)
        data_r = rfp_fintcells_leakcorr_ypos_ynorm{t,xy}(~rfp_fintcells_leakcorr_ypos_ynorm{t,xy}==0);
        data_t = tfp_fintcells_leakcorr_ypos_ynorm{t,xy}(~tfp_fintcells_leakcorr_ypos_ynorm{t,xy}==0);

        rfp_fintview_leakcorr_ypos_ynorm(t,xy) = mean(data_r);
        tfp_fintview_leakcorr_ypos_ynorm(t,xy) = mean(data_t);

        rfp_std_leakcorr_ypos_ynorm(t,xy) = std(data_r);
        tfp_std_leakcorr_ypos_ynorm(t,xy) = std(data_t);

        rfp_sem_leakcorr_ypos_ynorm(t,xy) = std(data_r) ./ sqrt(length(data_r));
        tfp_sem_leakcorr_ypos_ynorm(t,xy) = std(data_t) ./ sqrt(length(data_t));

        clearvars data_r data_t

    end
end
clearvars xy t ans data_r data_t


% 6) average fluo value of each expt condition (each q) after YFP normalization

for t = 1:total_time;
    for q = q_range_analyzed;
        
        % eliminate 0-value cells (i.e. cells that are YFP negative)
        r_data_temp = rfp_fintcells_q_leakcorr_ypos_ynorm{t,q}(~rfp_fintcells_q_leakcorr_ypos_ynorm{t,q}==0);
        t_data_temp = tfp_fintcells_q_leakcorr_ypos_ynorm{t,q}(~tfp_fintcells_q_leakcorr_ypos_ynorm{t,q}==0);
        
        rfp_fintview_q_leakcorr_ypos_ynorm(t,q) = mean(r_data_temp);
        tfp_fintview_q_leakcorr_ypos_ynorm(t,q) = mean(t_data_temp);
        
        rfp_std_q_leakcorr_ypos_ynorm(t,q) = std(r_data_temp);
        tfp_std_q_leakcorr_ypos_ynorm(t,q) = std(t_data_temp);
        
        rfp_sem_q_leakcorr_ypos_ynorm(t,q) = std(r_data_temp) ./ sqrt(length(r_data_temp));
        tfp_sem_q_leakcorr_ypos_ynorm(t,q) = std(t_data_temp) ./ sqrt(length(t_data_temp));
        
    end
end
    
clearvars r_data_temp t_data_temp t q


%% 7) collect all fluorescence intensities AFTER YFP normalized
% Collect all fluo intensities of all cells, spectral leakage-corrected AND YFP thresholded AND YFP normalized. 
% For scatter plot in the next section.

inputdata_r = rfp_fintcells_leakcorr_ypos_ynorm;
inputdata_t = tfp_fintcells_leakcorr_ypos_ynorm;
rfp_allcells_leakcorr_Ynorm = [];
tfp_allcells_leakcorr_Ynorm = [];
for xy = analyzed_xy;
    for t = 1:total_time;       
        rfp_allcells_leakcorr_Ynorm = [rfp_allcells_leakcorr_Ynorm inputdata_r{t,xy}];
        tfp_allcells_leakcorr_Ynorm = [tfp_allcells_leakcorr_Ynorm inputdata_t{t,xy}];
    end
end
clearvars inputdata_r inputdata_t xy t
save(strcat(analysis_kinetics,'\12_FPallcells_afterLeakCorr_Ythresh_',num2str(thresh_YFP),'_Ynorm.mat'),...
    'rfp_allcells_leakcorr_Ynorm','tfp_allcells_leakcorr_Ynorm'); % save variable


%% 8) visualize RFP vs. YFP and TFP vs. YFP of each cell AFTER YFP normalization
% OPTIONAL.

scatterPlot = figure('units','normalized','outerposition',[0 0 1 1]);

x_data = rfp_allcells_leakcorr_Ynorm;    x_text = 'RFP';
y_data = tfp_allcells_leakcorr_Ynorm;    y_text = 'TFP';

color = color_yfp_out;
markersize = 4;
title_text = {'Scatter Plot: each dot = one cell; all cells of all time AFTER YFP normalization',expt_date};
fig_savename = '\RFPvsTFP_scatter_leakcorr_afterYnorm.png';

leg_text = expt_date;

scatter(x_data,y_data,markersize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5);

% plot label
suptitle(title_text);
xlabel(x_text,'FontSize',20)
ylabel(y_text,'FontSize',20);

% legend
leg = leg_text; % legend
legend({leg},'FontSize',20);

% axis range
% axis([-0.4 0.6 -0.5 3]);

% save image
% print(scatterPlot,strcat(saveplace,fig_savename),'-dpng','-r300'); % resolution 300

clearvars scatterPlot x_data y_data x_text y_text color markersize
clearvars title-text fig_savename leg_text

%% 9) SAVE variables

savevars = {'rfp_fintcells_leakcorr_ypos_ynorm','tfp_fintcells_leakcorr_ypos_ynorm',...
            'rfp_fintcells_q_leakcorr_ypos_ynorm','tfp_fintcells_q_leakcorr_ypos_ynorm',...
            'rfp_fintview_leakcorr_ypos_ynorm', 'tfp_fintview_leakcorr_ypos_ynorm',...
            'rfp_std_leakcorr_ypos_ynorm','tfp_std_leakcorr_ypos_ynorm',...
            'rfp_sem_leakcorr_ypos_ynorm', 'tfp_sem_leakcorr_ypos_ynorm',...
            'rfp_fintview_q_leakcorr_ypos_ynorm','tfp_fintview_q_leakcorr_ypos_ynorm'...
            'rfp_std_q_leakcorr_ypos_ynorm','tfp_std_q_leakcorr_ypos_ynorm',...
            'rfp_sem_q_leakcorr_ypos_ynorm','tfp_sem_q_leakcorr_ypos_ynorm'};
        
% run at the end
save(strcat(analysis_kinetics,'/13_YFP_normalized_Ythresh_',num2str(thresh_YFP),'.mat'),savevars{:});
 

%% next: PlotGenerator
