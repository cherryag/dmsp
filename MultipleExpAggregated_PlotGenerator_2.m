%% Plot generator for multiple experiments aggregated
% Version July 2018

% combine v3 and v4
% 1) generate histograms of YFP (raw fluo, spectral leakage corrected, background subtracted)
% 2) normalize by 250 uM signal (not sure if I need to)
% 3) calculate average across experiments
% 4) visualize average across experiments (time vs. signal)

    % .date = date of experiment
    % .time = time points (hours since treatment)
    % .Ythresh = YFP threshold used
    
    % .Rc_raw = raw RFP fluorescence of all cells (background subtracted, YFP threshold applied, by q)
    % .Tc_raw = raw TFP fluorescence of all cells (background subtracted, YFP threshold applied, by q)
    % .Yc_raw = raw YFP fluorescence of all cells (background subtracted, YFP threshold applied, by q) -> this way, can do a histogram across all days
    % .Rm_raw
    % .Tm_raw
    % .Ym_raw
    % .Rstd_raw
    % .Tstd_raw
    % .Ystd_raw
    % .Rsem_raw
    % .Tsem_raw
    % .Ysem_raw
        
    % .Rc_Ynorm = RFP fluorescence of each cell normalized by its own YFP (background subtracted, YFP threshold applied)
    % .Tc_Ynorm = TFP fluorescence of each cell normalized by its own YFP (background subtracted, YFP threshold applied)
    % .Rm_Ynorm
    % .Tm_Ynorm
    % .Rstd_Ynorm
    % .Tstd_Ynorm
    % .Rsem_Ynorm
    % .Tsem_Ynorm
    
    % added:
    
    % .Yc_raw_agg = all YFP values of all cells, aggregated per experiment (for histogram)
    % .Rc_Ynorm_250norm
    

expt_range_processed = [5 10]; % list of experiments to process ; normally 1:length(expt_list)

% load publication colormap (Parula, 12/2018)
load('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/Matlab/plot color map/colormap_publication.mat')


%% 01) 2-panel figures for each experiment (dmdA and dddW) 
% added 2/14/2019: show mid-exponential time points on each panel
% figure 2
% all data are background-subtracted, spectral leakage-corrected, YFP-thresholded, YFP-normalized

% 1) 3cv3-TFP
% 2) 3cv3-RFP
% 3) 3cv4-TFP
% 4) 3cv4-RFP, all expts, all colors, all mean for each q

% 5) all of the above, mean across all expts
% 6) overlay with Pohnert 
% 7) normalize by initial cell number? 

% cd(motherdir); saveplace = 'aggregated_analysis/figure 2/expt_separated'; mkdir(saveplace);
close all;
SaveIm = 0; 
midpoint_show = 0; % show mid-exponential time point (load ratios_time_points.mat)

% figure for each experiment
for expt = expt_range_processed
        
    % initialize
    s = []; fp_color = {};
    
    % define labels and plot style
    q_range = 1:length(aggdata.xy_q{expt});
    cmap = parula(length(q_range));
    leg_string = aggdata.string_q{expt};
    expt_date = aggdata.date{expt};
    
    % define fluorescence colors according to strain
    if sum(leg_string{1}(1:4) == '3cv3') == 4
        fp_color{1} = 'TFP'; % color of demethylation
        fp_color{2} = 'RFP'; % color of cleavage
        data_plot_A = aggdata.Tm_Ynorm_t2{expt};
        data_plot_W = aggdata.Rm_Ynorm_t2{expt};
        error_plot_A = aggdata.Tsem_Ynorm_t2{expt};
        error_plot_W = aggdata.Rsem_Ynorm_t2{expt};
        if midpoint_show == 1
        midexpo_t_a = startexpo_t.tfp{expt};
        midexpo_t_w = startexpo_t.rfp{expt};
        end
    elseif sum(leg_string{1}(1:4) == '3cv4') == 4
        fp_color{1} = 'RFP'; % color of demethylation
        fp_color{2} = 'TFP'; % color of cleavage
        data_plot_A = aggdata.Rm_Ynorm_t2{expt};
        data_plot_W = aggdata.Tm_Ynorm_t2{expt};
        error_plot_A = aggdata.Rsem_Ynorm_t2{expt};
        error_plot_W = aggdata.Tsem_Ynorm_t2{expt};
        if midpoint_show == 1
        midexpo_t_a = startexpo_t.rfp{expt};
        midexpo_t_w = startexpo_t.tfp{expt};
        end
    else
        'experiment doesn''t match either of the strains'
    end
    
    % define plot variables
    t_range = aggdata.time{expt};

    
    % generate a figure
    if SaveIm == 0
        fig{expt} = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'on'); % figure pops up
    elseif SaveIm == 1
        fig{expt} = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off'); % no figure popup
       % fig{expt} = figure('units','normalized','outerposition',[0 0 0.5 1], 'visible', 'on'); % no figure popup       
    end
    
    % plot (each subplot = one pathway)
    for q = 1:9 %q_range % iterate through DMSP concentrations
        s(1) = subplot(121); % demethylation
            errorbar(t_range, data_plot_A(:,q), error_plot_A(:,q), '.-', 'markersize', 20, 'color', cmap_final(q,:)); hold on;
            
            if midpoint_show == 1
            plot(t_range(midexpo_t_a(q)) , data_plot_A(midexpo_t_a(q),q), 'r.', 'markersize',20);
            mid_y = [0:0.01:data_plot_A(midexpo_t_a(q),q)];
            t_plot = repmat(t_range(midexpo_t_a(q)),length(mid_y));
            plot(t_plot(:,1),mid_y, ':', 'linewidth', 5, 'color', cmap_final(q,:));
            end
             
        s(2) = subplot(122); % cleavage
            errorbar(t_range, data_plot_W(:,q), error_plot_W(:,q), '.-', 'markersize', 20, 'color', cmap_final(q,:)); hold on; 
            
            if midpoint_show == 1
            plot(t_range(midexpo_t_w(q)) , data_plot_W(midexpo_t_w(q),q), 'r.', 'markersize',20);
            mid_y = [0:0.01:data_plot_W(midexpo_t_w(q),q)];
            t_plot = repmat(t_range(midexpo_t_w(q)),length(mid_y));
            plot(t_plot(:,1),mid_y, ':', 'linewidth', 5, 'color', cmap_final(q,:));
            end
    end
    
    % labels
    xlabel(s(1), 'time since treatment (hours)'); 
    xlabel(s(2), 'time since treatment (hours)'); 
    ylabel(s(1), strcat('fluorescence (',fp_color{1},')'))
    ylabel(s(2), strcat('fluorescence (',fp_color{2},')'))
    title(s(1),'dmdA (demethylation)')
    title(s(2),'dddW (cleavage)')
    sgtitle({expt_date,'fluorescence is background-subtracted, spectral leakage-corrected, YFP-thresholded, and YFP-normalized (each cell)'})
 %   legend(s(1),leg_string,'FontSize',15,'location','northwest')
    ax = findall(gcf,'type','axes'); set(ax,'fontsize', 20); % set all axis labels the same size
    
    if SaveIm == 1
        savename = strcat('/startexpo_t_vs_fluo_2panel',expt_date(1:10));
        print(fig{expt},strcat(saveplace,savename),'-dpng','-r300'); % resolution 300
    end
    
end


%     xlabel('time since treatment (hours)'); 
%     ylabel(strcat('fluorescence (',fp_color{2},')'))
%   %  title('dmdA (demethylation)')
%  %   title(s(2),'dddW (cleavage)')
%     sgtitle({strcat(expt_date,';   dddW (cleavage)'),'fluorescence is background-subtracted, spectral leakage-corrected, YFP-thresholded, and YFP-normalized (each cell)'})
%     legend(s(1),leg_string,'FontSize',15,'location','northwest')
%     ax = findall(gcf,'type','axes'); set(ax,'fontsize', 16); % set all axis labels the same size


%% 01) [VERSION 2] separate figure for each color FP
% figure 2
% all data are background-subtracted, spectral leakage-corrected, YFP-thresholded, YFP-normalized

% 1) 3cv3-TFP
% 2) 3cv3-RFP
% 3) 3cv4-TFP
% 4) 3cv4-RFP, all expts, all colors, all mean for each q

% 5) all of the above, mean across all expts
% 6) overlay with Pohnert 
% 7) normalize by initial cell number? 

% cd(motherdir); saveplace = 'aggregated_analysis/figure 2/expt_separated'; mkdir(saveplace);
close all;
SaveIm = 0; 
pathway_string = {'demethylation (dmdA)','cleavage (dddW)'};
axis_lim_same = 1; % make y axis limits the same for each color-pathway

% figure for each experiment
for expt = 10% expt_range_processed
        
    % initialize
    s = []; fp_color = {};
    
    % define labels and plot style
    q_range = 1:9; % length(aggdata.xy_q{expt});
    leg_string = aggdata.string_q{expt};
    expt_date = aggdata.date{expt};
    
    % define fluorescence colors according to strain
    if sum(leg_string{1}(1:4) == '3cv3') == 4
        fp_color{1} = 'TFP'; % color of demethylation
        fp_color{2} = 'RFP'; % color of cleavage
        data_plot_A = aggdata.Tm_Ynorm_eacht{expt};
        data_plot_W = aggdata.Rm_Ynorm_eacht{expt};
        error_plot_A = aggdata.Tsem_Ynorm_eacht{expt};
        error_plot_W = aggdata.Rsem_Ynorm_eacht{expt};
%         data_plot_A = aggdata.Tm_Ynorm_t2{expt};
%         data_plot_W = aggdata.Rm_Ynorm_t2{expt};
%         error_plot_A = aggdata.Tsem_Ynorm_t2{expt};
%         error_plot_W = aggdata.Rsem_Ynorm_t2{expt};
    elseif sum(leg_string{1}(1:4) == '3cv4') == 4
        fp_color{1} = 'RFP'; % color of demethylation
        fp_color{2} = 'TFP'; % color of cleavage
        data_plot_A = aggdata.Rm_Ynorm_eacht{expt};
        data_plot_W = aggdata.Tm_Ynorm_eacht{expt};
        error_plot_A = aggdata.Rsem_Ynorm_eacht{expt};
        error_plot_W = aggdata.Tsem_Ynorm_eacht{expt};
%         data_plot_A = aggdata.Rm_Ynorm_t2{expt};
%         data_plot_W = aggdata.Tm_Ynorm_t2{expt};
%         error_plot_A = aggdata.Rsem_Ynorm_t2{expt};
%         error_plot_W = aggdata.Tsem_Ynorm_t2{expt};
    else
        'experiment doesn''t match either of the strains'
    end
    
    % define plot variables
    t_range = aggdata.time{expt};

        for pathway = 1:2 % a figure for each pathway
            % 1 = demethylation; 2 = cleavage

        % generate a figure
        if SaveIm == 0
            fig{expt} = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'on'); % figure pops up
        elseif SaveIm == 1
       %     fig{expt} = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off'); % no figure popup
            fig{expt} = figure('units','normalized','outerposition',[0 0 0.5 1], 'visible', 'off'); % no figure popup       
        end

        % plot (each subplot = one pathway)
        for q = q_range % iterate through DMSP concentrations
            if pathway == 1 % demethylation
                errorbar(t_range, data_plot_A(:,q), error_plot_A(:,q), '.-', 'markersize', 20, 'color', cmap_final(q,:)); hold on;
            elseif pathway == 2 % cleavage
                errorbar(t_range, data_plot_W(:,q), error_plot_W(:,q), '.-', 'markersize', 20, 'color', cmap_final(q,:)); hold on; 
            end
        end

        % labels [NEEDS TO FIX]
        xlabel('time since treatment (hours)'); 
        ylabel(strcat('fluorescence (',fp_color{pathway},')'))
        sgtitle({strcat(expt_date,';   ', pathway_string{pathway}),...
                 'fluorescence is background-subtracted, spectral leakage-corrected, YFP-thresholded, and YFP-normalized (average YFP at each timepoint)'}) % ('t2 value for each q')
    %    legend(leg_string,'FontSize',15,'location','northwest')
        ax = findall(gcf,'type','axes'); set(ax,'fontsize', 20); % set all axis labels the same size
        xlim([0 24])
        
        if axis_lim_same == 1
            if pathway == 1 && strcmp(fp_color{pathway}, 'TFP') % dmdA, TFP
                ylim([0 0.3])
            elseif pathway == 1 && strcmp(fp_color{pathway},'RFP') % dmdA, RFP
                ylim([0 0.065]) 
                % ylim([0 0.06]) % for fig. 2
            elseif pathway == 2 && strcmp(fp_color{pathway},'TFP') % dddW, TFP
                ylim([0 0.55])
                % ylim([0 0.45]) % for fig. 2
            elseif pathway == 2 && strcmp(fp_color{pathway},'RFP') % dddW, RFP
                ylim([0 3.2])
            end
        end
        
        if SaveIm == 1
%            savename = strcat('/t_vs_fluo_Ynormt2_',pathway_string{pathway}((end-4):(end-1)),'_',expt_date(1:10));
            savename = strcat('/t_vs_fluo_Ynorm_Yavg_each_t_',pathway_string{pathway}((end-4):(end-1)),'_',expt_date(1:10));
           % print(fig{expt},strcat(saveplace,savename),'-dpng','-r300'); % resolution 300
            print(fig{expt},strcat(saveplace,savename), '-depsc'); % for publication
           saveas(fig{expt},strcat(saveplace,savename)); % as .fig 
        end
    
    
    end
    
end


%% 02) [PART 1] all experiments on the same plot (4 plots, RFP-dddW, TFP-dddW, RFP-dmdA, TFP-dmdA)
% each figure: one strain, one color, all expts

% cd(motherdir); saveplace = 'aggregated_analysis/figure 2/expt_together'; mkdir(saveplace);
close all;
SaveIm = 0; 
strain_list = {'3cv3','3cv4'};
color_list = {'TFP','RFP'};

% make strain_dates (separate dates for 3cv3 and 3cv4)
j = 1; jj = 1; % initialize
strain_dates = {{},{}}; strain_index = {[],[]}; % 3cv3 and 3cv4 dates
for i = 1:length(expt_list)
    
    if strcmp('(3cv3)',expt_list{i}(end-5:end)) == 1
        strain_dates{1}{i} = expt_list{i};
        strain_index{1}(j) = i;
            j = j + 1;
        
    elseif strcmp('(3cv4)',expt_list{i}(end-5:end)) == 1
        strain_dates{2}{jj} = expt_list{i};
        strain_index{2}(jj) = i;
            jj = jj + 1;
        
    end
    
end


%% 02) [PART 3] all experiments on the same plot (4 plots, RFP-dddW, TFP-dddW, RFP-dmdA, TFP-dmdA)
% different options for normalization -> none worked
close all
marker_expt = {'.-','o-','*-','d-','s-','^-','>-','+-','p-','h-'}; % different marker for each experiment

% strain_index = {[5,6,7],[8,9,10]}; % after looking at the cell numbers

cellnum_plot = 0; 
SaveIm = 0; 

% loop through for plot
for s = 1:length(strain_list) % for each strain
    strain = strain_list{s}; % specify 3cv3 or 3cv4

    for fp = 1%:length(color_list) % for each FP color
        fp_color = color_list{fp}; % tfp or rfp
        
        % generate a figure for each strain and FP color
        if SaveIm == 0
            figure('units','normalized','outerposition',[0 0 0.5 1], 'visible', 'on'); % figure pops up
        elseif SaveIm == 1
            figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off'); % no figure popup
        end
               
        
        for expt = 8%strain_index{s} % for each experiment
     %   figure
            % determine strain, FP color, pathway to plot according to strain & FP color
            if fp_color == 'TFP'
                data_plot = aggdata.Tm_Ynorm{expt};%./ mean(slope_growth{expt});% ./ mean(slope_growth{expt}) ./ mean(slope_growth{expt}); %./ aggdata.Tm_Ynorm{expt}(32,2)./ slope_growth(expt); %aggdata.cell_num_ypos_avg_q{expt}(2,4);
                error_plot = aggdata.Tsem_Ynorm{expt};%./ mean(slope_growth{expt});% ./ mean(slope_growth{expt})./ mean(slope_growth{expt});% ./slope_growth(expt)  ./ aggdata.Tm_Ynorm{expt}(32,2)./ slope_growth(expt);%./ aggdata.cell_num_ypos_avg_q{expt}(2,4);
                if cellnum_plot == 1
                    data_plot = aggdata.cell_num_ypos_avg_q{expt};
                    error_plot = aggdata.cell_num_ypos_std_q{expt};
                end

                if strain == '3cv3'
                    pathway = 'dmdA';
                elseif strain == '3cv4'
                    pathway = 'dddW';
                end

            elseif fp_color == 'RFP'
                data_plot = aggdata.Rm_Ynorm{expt};%./ mean(slope_growth{expt}); % ./slope_growth(expt) %./ aggdata.Rm_Ynorm{expt}(32,2)./ slope_growth(expt);%./ aggdata.cell_num_ypos_avg_q{expt}(2,4);
                error_plot = aggdata.Rsem_Ynorm{expt};%./ mean(slope_growth{expt}); %./slope_growth(expt) ./ aggdata.Rm_Ynorm{expt}(32,2)./ slope_growth(expt);% aggdata.cell_num_ypos_avg_q{expt}(2,4);
                if cellnum_plot == 1
                    data_plot = aggdata.cell_num_ypos_avg_q{expt};
                    error_plot = aggdata.cell_num_ypos_std_q{expt};
                end
                if strain == '3cv3'
                    pathway = 'dddW';
                elseif strain == '3cv4'
                    pathway = 'dmdA';
                end

            end 
                                
            % define labels and plot style
            q_range = 1:9;%[1 2 4 7 9]; %1:9; % 1:length(aggdata.xy_q{expt});
            cmap = parula(9); %parula(length(q_range));
            leg_string = aggdata.string_q{expt};
            expt_date = aggdata.date{expt};
            t_range = aggdata.time{expt};
               
            % plot (each subplot = one pathway)
         
            for q = q_range % iterate through DMSP concentrations
             
                
        % if expt_separate == 1  % not the mean
              if ismember(q,[2:9])
                    errorbar(t_range, data_plot(:,q), error_plot(:,q), marker_expt{expt} , 'markersize', 20, 'color', cmap(q,:),...
                        'HandleVisibility','off'); hold on;
              elseif ismember(q,1) % for legend
                    errorbar(t_range, data_plot(:,q), error_plot(:,q), marker_expt{expt} , 'markersize', 20, 'color', cmap(q,:),...
                        'HandleVisibility','on','DisplayName',aggdata.date{expt}); hold on;
              end
              
              sgtitle({strcat('strain: ',strain,';       FP color: ',fp_color,';       pathway: ',pathway),...
                    'each marker is a different experiment'}, 'fontsize', 16)
            xlabel('time since treatment (hours)');
            ylabel({'fluorescence intensity', '(bckg subtracted; spectral leakage corrected; YFP thresholded at 30; YFP normalized)'})
            legend('Fontsize',16,'location','northwest')
            set(gca,'Fontsize',16)             
             
       %  elseif expt_separate == 0
      %       data_plot_mean(
         end
            
            
          end
            

            
        end % end of experiment
        
      
        
            % save image (each strain-color combination)
            if SaveIm == 1
                savename = strcat('/t_vs_fluo_',strain,'_',fp_color);
                print(strcat(saveplace,savename),'-dpng','-r300'); % resolution 300
            end
            
    end % end of fp color
%end % end of strain


%% 03) Mean fluorescence across experiments (4 plots, RFP-dddW, TFP-dddW, RFP-dmdA, TFP-dmdA)
% updated: 2/1/2019
% to make a plot with mean across all experiments (need to be finalized)
% average across experiments
close all
strain_index = {[5,6,7],[8,9,10]}; % after looking at the cell numbers
strain_list = {'3cv3','3cv4'};
color_list = {'RFP','TFP'};
% saveplace = uigetdir;
cellnum_plot = 0; 
% loop through for plot
SaveIm = 0;

for s = 1:2 % length(strain_list) % for each strain
    strain = strain_list{s}; % specify 3cv3 or 3cv4

    for fp = 1:length(color_list) % for each FP color
        fp_color = color_list{fp}; % tfp or rfp
        
        % generate a figure for each strain and FP color
        if SaveIm == 0
            figure('units','normalized','outerposition',[0 0 0.5 1], 'visible', 'on'); % figure pops up
        elseif SaveIm == 1
%            figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off'); % no figure popup
            figure('units','normalized','outerposition',[0 0 0.5 1], 'visible', 'off'); % no figure popup

        end               
        
        % determine strain, FP color, pathway to plot according to strain & FP color
        if fp_color == 'TFP'
            
            data_plot = aggdata.Tm_Ynorm_t2;
            error_plot = aggdata.Tsem_Ynorm_t2;
            if cellnum_plot == 1
                data_plot = aggdata.cell_num_ypos_avg_q;
                error_plot = aggdata.cell_num_ypos_std_q;
            end

            if strain == '3cv3'
                pathway = 'dmdA';
            elseif strain == '3cv4'
                pathway = 'dddW';
            end

        elseif fp_color == 'RFP'
            data_plot = aggdata.Rm_Ynorm_t2; 
            error_plot = aggdata.Rsem_Ynorm_t2; 
            if cellnum_plot == 1
                data_plot = aggdata.cell_num_ypos_avg_q;
                error_plot = aggdata.cell_num_ypos_std_q;
            end
            if strain == '3cv3'
                pathway = 'dddW';
            elseif strain == '3cv4'
                pathway = 'dmdA';
            end

        end 
                                
        % define labels and plot style
        q_range = 1:9; % 1:length(aggdata.xy_q{expt});
       % cmap = parula(9); %parula(length(q_range));
        leg_string = aggdata.string_q{5};
 %       expt_date = aggdata.date{expt};
        t_range = aggdata.time{5}; % need to FIX

        % plot (each subplot = one pathway)
        data_plot_mean = []; error_plot_mean = [];
        for q = q_range % iterate through DMSP concentrations
            for tt = 1:length(t_range) 
                data_temp = [];
                for date_iterate = strain_index{s}

                    data_temp = [data_temp data_plot{date_iterate}(tt,q)];
                end
               data_plot_mean(tt,q) =  mean(data_temp);
               error_plot_mean(tt,q) = std(data_temp)/sqrt(length(data_temp));
            end
            errorbar(t_range, data_plot_mean(:,q), error_plot_mean(:,q), '.-', 'markersize', 20, 'color', cmap_final(q,:)); hold on
        end
                           
        sgtitle({strcat('strain: ',strain,';       FP color: ',fp_color,';       pathway: ',pathway),...
                'each marker an average across 3 experiments (errorbar = SEM)'}, 'fontsize', 16)
        xlabel('time since treatment (hours)');
        xlim([0 24])
        ylabel({'fluorescence intensity', '(bckg subtracted; spectral leakage corrected; YFP thresholded at 50; YFP normalized using t2 value)'})
        legend(leg_string,'Fontsize',16,'location','northwest')
        set(gca,'Fontsize',16)             
        
        % save image (each strain-color combination)
        if SaveIm == 1
            savename = strcat('/t_vs_fluo_',strain,'_',fp_color,'_avgexpt_resized');
            print(strcat(saveplace,savename),'-dpng','-r300'); % resolution 300
        end
            
    end
end % end of experiment
        
      
%% 04) growth curve [2/7/2019]

% collect cell numbers and calculate normalized cell numbers (by t1 cell number)
cell_num_ypos_q = {}; 
cell_num_ypos_q_norm_by_t1 = {};
for expt = 5:10
    for q = 1:9
        for t = 1:32
                       
            cell_num_ypos_q{expt}(t,q) = sum(~isnan(aggdata.Rc_Ynorm_t2{expt}{t,q}));
            cell_num_ypos_q_norm_by_t1{expt}(t,q) = cell_num_ypos_q{expt}(t,q) / cell_num_ypos_q{expt}(1,q);
                        
        end
    end
end

save('cell_number_ypos_q_for_growth_curve.mat','cell_num_ypos_q','cell_num_ypos_q_norm_by_t1');

%% growth curve plots [polished 5/21/2019]
close all
% saveplace = 'growth curves figures/mean FOV_no_t1/same axis'; mkdir(saveplace)
SaveIm = 1;

% indicate what to plot
plot_norm_by_t1 = 0; 
plot_sum_cell_num = 0; 
plot_mean_cell_num = 1;

% define data to plot
if plot_norm_by_t1 == 1
    plot_data = cell_num_ypos_q;  
elseif plot_sum_cell_num == 1
    plot_data = cell_num_ypos_q_norm_by_t1;
elseif plot_mean_cell_num == 1
    plot_data = aggdata.cell_num_ypos_avg_q;
    plot_std = aggdata.cell_num_ypos_std_q;
end

% define time range to plot
plot_trange = 2:32; % all time ranges have length 32

% generate a plot for each experiment
for expt = 5:10
    
    t_range = aggdata.time{expt}(plot_trange); % aggdata.time{expt}
    figure('units','normalized','outerposition',[0 0 0.5 1], 'visible', 'on'); 
    
    for q = 1:9
        
        if plot_norm_by_t1 == 1 || plot_sum_cell_num == 1 % no error bars
            plot(t_range, plot_data{expt}(plot_trange,q), '.-', 'markersize', 20, 'color', cmap_final(q,:)); hold on;        
            ytext = 'number of cells (cumulative for each q, across 7 xy positions)'; % y label
            
            % determine y lim
            if plot_norm_by_t1 == 1
                y_lim_range = [0 22]; % for normalized, including outlier
  %              y_lim_range = [0 7]; % for normalized, excluding outlier
            elseif plot_sum_cell_num ==1
                y_lim_range = [100 8000]; % for absolute cell counts
            end
            
        elseif plot_mean_cell_num == 1 % for average FOV; error bars
            errorbar(t_range, plot_data{expt}(plot_trange,q), plot_std{expt}(plot_trange,q),...
                '.-', 'markersize', 20, 'color', cmap_final(q,:)); hold on;        
            ytext = 'mean number of cells (across 7 xy positions after Ythresh; error bars are std)'; % y label
            y_lim_range = [50 1200]; % for average FOV
        end
        
    end % end of q
    
        % label plots
        xlabel('time since treatment (hours)'); 
        ylabel(ytext)
        title(strcat('growth curves;  ', aggdata.date{expt}));
        set(gca,'fontsize', 16);

        % axis limits
        xlim([0 24])
        ylim(y_lim_range)


        % save figure
        if SaveIm == 1
            print(gcf,strcat(saveplace,'/cell_num_normalized_',aggdata.date{expt}(1:10)),'-depsc')
            close all
        end
    
end

%% 02) mean number of cells in a field of view across replicates [6/4/2019 added]
% calculation for figure caption

t = 2; 

mean( [mean(aggdata.cell_num_ypos_avg_q{5}(t,:)), mean(aggdata.cell_num_ypos_avg_q{6}(t,:)), mean(aggdata.cell_num_ypos_avg_q{7}(t,:)), ...
       mean(aggdata.cell_num_ypos_avg_q{8}(t,:)), mean(aggdata.cell_num_ypos_avg_q{9}(t,:)), mean(aggdata.cell_num_ypos_avg_q{10}(t,:)) ] )
   
std(  [mean(aggdata.cell_num_ypos_avg_q{5}(t,:)), mean(aggdata.cell_num_ypos_avg_q{6}(t,:)), mean(aggdata.cell_num_ypos_avg_q{7}(t,:)), ...
       mean(aggdata.cell_num_ypos_avg_q{8}(t,:)), mean(aggdata.cell_num_ypos_avg_q{9}(t,:)), mean(aggdata.cell_num_ypos_avg_q{10}(t,:)) ] )


%% 02) [PART 2] growth curve (quantify slope of growth curve)
    % max value - min value (not t1) / time(max) - time(min)
% future: quantify slope properly 
% q = 5,6,7
% (square, ^, >)
slope_growth = {};
for expt = 5:10
    for q = 1:9
        data = [];
        data = aggdata.cell_num_ypos_avg_q{expt}(:,q);
        t_range = aggdata.time{expt};
        
        [value_max,index_max] = max(data);
        [value_min,index_min] = min(data(7:end)) ;  % not including t = 1
        slope_growth{expt}(q) = ( value_max - value_min ) / ( t_range(index_max) - t_range(index_min+1) ); % index_min+1 because not including t =1 

    end
end

figure
for expt = 5:10
    for q = 1:9
        
        plot(q,slope_growth{expt}(q), marker_expt{expt},'markersize',15,'color',cmap_final(q,:)); 
        hold on;
    
    end
end

figure
for expt = 5:7
    
        
        plot(1:9,slope_growth{expt}, marker_expt{expt},'markersize',15)%,'color',cmap_final(q,:)); 
        hold on;
    
    
end



%% 04) slope of growth vs. slope of fluo or max fluo


marker_expt = {'.-','o-','*-','d-','s-','^-','>-','+-','p-','h-'};

% MAX FLUO VALUES
figure;
for expt = 5:7% 8:10 % 5:7
   % figure
   % for q = 2:9
    data = []; max_fluo = [];
    data = aggdata.Rm_Ynorm{expt}(:,1:9);
    max_fluo = max(data);
   

   %  end
   
   for q = 2:8
        plot(slope_growth{expt}(q),max_fluo(q),marker_expt{expt},'markersize',15,'color',cmap_final(q,:)); hold on;
   end

end

xlabel('slope of cell number increase')
ylabel('max fluo (RFP, dmdA)')
%title('slope of cell number increase vs. maximum fluorescence value (TFP, dmdA, 3cv3); square = 8/10/2018; ^ = 8/14/2018; > = 8/16/2018')
title('slope of cell number increase vs. maximum fluorescence value (RFP, dmdA, 3cv4); + = 2/28/2018; star = 3/1/2018; hexa = 3/3/2018')

set(gca,'fontsize',16)


% SLOPE OF FLUO INCREASE

% calculate slope of fluo increase
slope_fluo_increase = {};
for expt = 5:10
    for q = 1:9
        data = [];
        data = aggdata.Rm_Ynorm{expt}(:,q);
        t_range = aggdata.time{expt};
        
        [value_max,index_max] = max(data);
        [value_min,index_min] = min(data) ;  % not including t = 1
        slope_fluo_increase{expt}(q) = ( value_max - value_min ) / ( t_range(index_max) - t_range(2) ); % index_min+1 because not including t =1 

    end
end

% plot slope of fluo increase
   figure
for expt = 5:7

    for q = 2:9
  %  data = []; max_fluo = [];
   % data = aggdata.Rm_Ynorm{expt}(:,1:9);
   % max_fluo = max(data);
  
    
  
    plot(slope_growth{expt}(q),slope_fluo_increase{expt}(q),marker_expt{expt},'markersize',15,'color',cmap_final(q,:)); hold on;


    end
end

xlabel('slope of cell number increase')
ylabel('slope of fluo increase (RFP, dddW)')
title('slope of cell number increase vs. maximum fluorescence value (RFP, dddW, 3cv3); square = 8/10/2018; ^ = 8/14/2018; > = 8/16/2018')
%title('slope of cell number increase vs. maximum fluorescence value (RFP, dmdA, 3cv4); + = 2/28/2018; star = 3/1/2018; hexa = 3/3/2018')

set(gca,'fontsize',16)




%% 2) normalize by average fluo value of 3 time points of 250 uM condition
% DO I NEED TO DO THIS?!
% 250 uM condition is chosen because it always reliably saturates.
% normalize each experiment by its own 250 uM signal.
% normalize each color (i.e. each pathway) by its own 250 uM signal. -> if so, cannot compare between colors! 
% normalize each time point
% .Rc_Ynorm_250norm

% for 3cv3 and 3cv4 (use aggdata)
for i = 1:4;
    Rdata = aggdata.Rm_Ynorm{i} % replace with Rm_Ynorm
    Tdata = aggdata.Tm_Ynorm{i} % replace with Tm_Ynorm
        
    norm_value_R = mean(Rdata((end-2):end,3)) % last 3 time points (t30-32), 250uM condition (q = 3)
    norm_value_T = mean(Tdata((end-2):end,3)) % last 3 time points (t30-32), 250uM condition (q = 3)
    
    Rdata_250norm_temp = [];
    Tdata_250norm_temp = [];
    
    for q = 1:9
        for t = 1:32;
            Rdata_250norm_temp(t,q) = Rdata(t,q)./norm_value_R;  
            Tdata_250norm_temp(t,q) = Tdata(t,q)./norm_value_T;         
        end
    end
      
    aggdata.Rc_Ynorm_250norm{i} = Rdata_250norm_temp;
    aggdata.Tc_Ynorm_250norm{i} = Tdata_250norm_temp;

end
clearvars i t q Rdata_250norm_temp 


%% 3) calculate average across experiments 
datasets = {aggdata_v3 aggdata_v4};

for data = 1:2; % v3 and v4
for t = 1:32
    for q = 1:9;
        Rdata_temp = []; Tdata_temp = [];
        for i = 1:3;
            Rdata_temp = [Rdata_temp datasets{data}.Rm{i}(t,q)];
            Tdata_temp = [Tdata_temp datasets{data}.Tm{i}(t,q)];
        end
        
        % calculate mean for each t and q
        if data == 1;
        aggdata_v3.Rm_expt(t,q) = mean(Rdata_temp)
        aggdata_v3.Rstd_expt(t,q) = std(Rdata_temp)      
        aggdata_v3.Tm_expt(t,q) = mean(Tdata_temp)
        aggdata_v3.Tstd_expt(t,q) = std(Tdata_temp)      
        end
        if data == 2;
        aggdata_v4.Rm_expt(t,q) = mean(Rdata_temp)
        aggdata_v4.Rstd_expt(t,q) = std(Rdata_temp)      
        aggdata_v4.Tm_expt(t,q) = mean(Tdata_temp)
        aggdata_v4.Tstd_expt(t,q) = std(Tdata_temp)      
        end
    end
end
end



%% 4) visualize (need to make better)

cmap = parula(17);
figure('units','normalized','outerposition',[0 0 0.25 0.8]);
%for t_range = 1:32;
    t = 1:32; %aggdata_v4{1}.time(:,1);
    for q = [1 2 4 6];
        % dddW = v3 RFP
        y1 = aggdata_v3.Rm_expt(t,q);
        y1_sem = aggdata_v3.Rstd_expt(t,q);
        
                
        % dmdA = v4 RFP
      %  for i = 1:3;
       y2 = aggdata_v4.Rm_expt(t,q);       
        y2_sem = aggdata_v4.Rstd_expt(t,q);
        
        if q == 1;
            c = cmap(10+3,:);
        else
            c = cmap(q+3,:);
        end
        
      %  yyaxis left
        errorbar(t,y1,y1_sem,'.','MarkerSize',50,'color',c,'linewidth',3); hold on;
     %   yyaxis right
     %   errorbar(t,y2,y2_sem,'o','MarkerSize',15,'color',c,'linewidth',3); hold on;
      %  end
    end
%end
xlim([0 24]);
ylim([0 3])
legend({'1 mM glucose','1 mM DMSP','100 \muM DMSP','50 \muM DMSP',},'FontSize',20,'location','northwest')
%title('3cv3 RFP (dddW), 2/25/2018 experiment, YFP-normalized, YFP threshold = 50','FontSize',10)
suptitle({'3cv4 RFP (dmdA), 2/28/2018 experiment, YFP-normalized,',... 
        'YFP threshold = 50'})
    
suptitle({'3cv3 RFP (dddW), 2/25/2018 experiment, YFP-normalized,',... 
        'YFP threshold = 50'})    
xlabel('time (hours)')
ylabel('dmdA expression (A.U.)')
set(gca, 'FontSize',20)
grid on

plot(DMSP(q),aggdata_v3.raw_ratio{1,1}(t,q),'.','MarkerSize',75,'color',cmap(q+3,:),'linewidth',3); hold on;
plot(DMSP(q),aggdata_v4.raw_ratio{1,1}(t,q),'o','MarkerSize',20,'color',cmap(q+3,:),'linewidth',3); hold on;

%% visualize normalized values (fluo vs. time for all expts)
% use parula with filled vs. not filled colors (like in poster)

plot_dmdA_v3 = 1;
plot_dddW_v3 = 0;
plot_dmdA_v4 = 0;
plot_dddW_v4 = 1;

notnorm = 1; % fp_fintview_q_leakcorr_ypos_ynorm
sub = 0; % subtracted
norm = 0; % subtracted and normalized

% manual 
title_text = ['3cv3-dmdA (TFP) and 3cv4-dmdA (RFP) of all experiments; each data set subtracted by its own t1 value',...
            'each data set subtracted by its own t1 value; then, all values in each experiment is normalized by average fluo intensity of last 5 time points of the 1 mM DMSP condition'];

marker_s = {'*-','*--','*:'};
marker_o = {'o-','o--','o:'};

figure('units','normalized','outerposition',[0 0 1 1]);
cmap = hsv(9);
for i = 1:3
    for q = [1 2 4 6];
        
        if plot_dmdA_v4 == 1;
            time_data = aggdata_v4.time{i}(:,1);
            if notnorm == 1;
                plot_data = aggdata_v4.Rm{i}(:,q);
            elseif sub == 1;               
                plot_data = aggdata_v4.Rm_sub{i}(:,q);
            elseif norm == 1;
                plot_data = aggdata_v4.Rc_Ynorm_250norm{i}(:,q);
            end
            plot(time_data,plot_data,marker_s{i},'MarkerSize',15,'Color',cmap(q,:),'LineWidth',2); hold on;               
        end
        
        if plot_dddW_v4 == 1;
            time_data = aggdata_v4.time{i}(:,1);
            if notnorm == 1;
                plot_data = aggdata_v4.Tm{i}(:,q);
            elseif sub == 1;  
                plot_data = aggdata_v4.Tm_sub{i}(:,q);
            elseif norm == 1;
                plot_data = aggdata_v4.Tc_Ynorm_250norm{i}(:,q);
         end
            plot(time_data,plot_data,marker_s{i},'MarkerSize',15,'Color',cmap(q,:),'LineWidth',2); hold on;          
        end        
        
        if plot_dmdA_v3 == 1;
            if i == 2 | i == 3;
                continue
            end
            time_data = aggdata_v3.time{i}(:,1);
            if notnorm == 1;
                 plot_data = aggdata_v3.Tm{i}(:,q);
            elseif sub == 1;  
                 plot_data = aggdata_v3.Tm_sub{i}(:,q);
            elseif norm == 1;
                 plot_data = aggdata_v3.Tc_Ynorm_250norm{i}(:,q);
           end
            plot(time_data,plot_data,marker_o{i},'MarkerSize',15,'Color',cmap(q,:),'LineWidth',2); hold on;       
        end        
        
        if plot_dddW_v3 == 1;
            if i == 2 | i == 3;
                continue
            end
            time_data = aggdata_v3.time{i}(:,1);
            if notnorm == 1;
                plot_data = aggdata_v3.Rm{i}(:,q);
            elseif sub == 1;  
                plot_data = aggdata_v3.Rm_sub{i}(:,q);
            elseif norm == 1;
                plot_data = aggdata_v3.Rc_Ynorm_250norm{i}(:,q);
               end
            plot(time_data,plot_data,marker_o{i},'MarkerSize',15,'Color',cmap(q,:),'LineWidth',2); hold on;            
        end
        
    end     
end

suptitle({title_text,'o = 3cv3;   * = 3cv4','glucose, DMSP 1mM, DMSP 250 uM, DMSP 75 uM, DMSP 1 uM conditions only'})
xlabel('time since treatment (hours)');
ylabel('fluorescence intensity, YFP-normalized');
legend({'glucose 1 mM',...
       'DMSP 1 mM',...
       'DMSP 250 \muM',...
       'DMSP 100 \muM',...
       'DMSP 75 \muM',...
       'DMSP 50 \muM',...
       'DMSP 25 \muM',...
       'DMSP 10 \muM',...
       'DMSP 1 \muM'},...
       'location','northwest');
set(gca,'FontSize',20);
xlim([0 24])
%%
print(strcat(saveplace,'/dmdA_3cv3_3cv4_normalized_sparse'),'-dpng','-r300'); % resolution 300

        


%% [added 8/15/2018, from PlotGenerator_4] plot xy vs. time - generic
% For scenario where there are multiple curves (for multiple xy positions) per expt condition.
    % i.e. 7 curves per color per subplot

PlotYFP = 1; % 1 = plot YFP data; 0 = no YFP data plotted but still need to give dummy data
ErrorBarOn = 1; % 1 = with errorbar; 0 = no errorbar, but still need to give dummy error data
SaveIm = 0; 

q_range = 1:9;%q_range_analyzed; % will determine xy data from expt_conditions_xy
i = 1;

% define data to plot (multiple data array per expt condition)
plot_rfp =  aggdata.Rc_Ynorm{i};
plot_tfp = aggdata.Tm_Ynorm{i};
plot_yfp = 0%aggdata.Yc_Ynorm{i};
error_rfp = aggdata.Rsem_Ynorm{i};
error_tfp = aggdata.Tsem_Ynorm{i};
error_yfp = 0;
time_expt_hours = aggdata.time{i};
expt_conditions_string = {'3cv3 - glucose 1 mM',...
                          '3cv3 - DMSP 1 mM',...
                          '3cv3 - DMSP 250 \muM',...
                          '3cv3 - DMSP 100 \muM',...
                          '3cv3 - DMSP 75 \muM',...
                          '3cv3 - DMSP 50 \muM',...
                          '3cv3 - DMSP 25 \muM',...
                          '3cv3 - DMSP 10 \muM',...
                          '3cv3 - DMSP 1 \muM'}; % for legend


% plot labels
title_text = {sprintf('3cv3 DMSP kinetics (leakage corrected; YFP threshold = %.f); each point = avg fluo of cells in image',aggdata.Ythresh{i})...
              strcat('error bars = SEM;    expt date:   ',aggdata.date{i},')')};
x_text = 'time from treatment (hours)';
y_text = 'fluorescence (A.U.)';
leg_text = {'RFP (dddW)',...
            'TFP (dmdA)',...
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


%% 4) plot q vs. time - generic
% For scenario where there is one curve per expt condition.
    % i.e. 1 curve per color per subplot
% h = {}
PlotYFP = 0; % 1 = plot YFP data; 0 = no YFP data plotted but still need to give dummy data
ErrorBarOn = 1; % 1 = with errorbar; 0 = no errorbar, but still need to give dummy error data

q_range = 1:total_cond; % will determine xy data from expt_conditions_xy

for i = 1:3;

% define data to plot (1 data array per expt condition)
plot_rfp =  aggdata.Rm_Ynorm{i} ./ aggdata.Rm_Ynorm{i}(30,3);
plot_tfp = aggdata.Tm_Ynorm{i} ./ aggdata.Rm_Ynorm{i}(30,3);
plot_yfp = 0%aggdata.Yc_Ynorm{i};
error_rfp = aggdata.Rsem_Ynorm{i};
error_tfp = aggdata.Tsem_Ynorm{i};
error_yfp = 0;
time_expt_hours = aggdata.time{i};


% plot labels
x_text = 'time from treatment (hours)';
y_text = 'fluorescence (A.U.)';

% y_text = 'RATIO RFP/TFP (dddW/dmdA)';
if ismember(i,[1:3]) == 1; % 3cv3
             leg_text = {'RFP (dddW)',...
            'TFP (dmdA)',...
            'YFP (constitutive)'};
        title_text = {sprintf('3cv3 DMSP kinetics (leakage corrected; YFP threshold = %.f); each point = avg fluo of cells in image',aggdata.Ythresh{i})...
              strcat('both RFP and TFP values are divided by t30 RFP value of DMSP 250 uM (q3); error bars = SEM;    expt date:   ',aggdata.date{i},')')};

          elseif ismember(i,[4:6]) == 1; % 3cv4
             leg_text = {'RFP (dmdA)',...
            'TFP (dddW)',...
            'YFP (constitutive)'};
         title_text = {sprintf('3cv4 DMSP kinetics (leakage corrected; YFP threshold = %.f); each point = avg fluo of cells in image',aggdata.Ythresh{i})...
              strcat('both RFP and TFP values are divided by t30 TFP value of DMSP 250 uM (q3); error bars = SEM;    expt date:   ',aggdata.date{i},')')};

          end


[h{i} q_t] = plot_q_t(time_expt_hours, PlotYFP, ErrorBarOn, q_range,...
                   expt_conditions_string, expt_conditions_xy,...
                   plot_rfp, plot_tfp, plot_yfp, error_rfp, error_tfp, error_yfp,...
                   title_text, x_text, y_text, leg_text);
%%
axis(h{i}(:),[0 24 0 3.5]);   

% save figure
fig_savename = strcat('\AvgFluo_q_Ynorm_D250RFPnorm_i_',num2str(i),'.png');
print(q_t,strcat(saveplace,fig_savename),'-dpng','-r150'); % resolution 300


end

                            
clearvars PlotYFP ErrorBarOn plot_rfp plot_tfp plot_yfp error_rfp error_tfp error_yfp
clearvars title_text x_text y_text leg_text fig_savename
 
      
%% two-sample t-test
% at each time point of glucose condition, compare the same color in different strains (after YFP normalization):
    % plot 1) TFP: dddW in 3cv4 and dmdA in 3cv3
    % plot 2) RFP: dmdA in 3cv4 and dddW in 3cv3


%%
saveplace = uigetdir; % data_figure folder

% plot fluorescence over time (glucose condition only)
% TFP
fluo_time_TFP = figure('units','normalized','outerposition',[0 0 1 1]);
marker_s = {'g*-','g*--','g*:'}
marker_o = {'go-','go--','go:'}

q = 2; % for glucose condition
for i = 1:3

    errorbar(aggdata_v4.time{i}(:,q),...
             aggdata_v4.Tm{i}(:,q),...
             aggdata_v4.Tsem{i}(:,q),...
             marker_s{i},'MarkerSize',15,'LineWidth',2);
    hold on;

    errorbar(aggdata_v3.time{i}(:,q),...
             aggdata_v3.Tm{i}(:,q),...
             aggdata_v3.Tsem{i}(:,q),...
             marker_o{i},'MarkerSize',10,'LineWidth',2);
     
end
leg = {strcat('3cv4-TFP (dddW)__',expt_date_list_v4{1}),...
       strcat('3cv3-TFP (dmdA)__',expt_date_list_v3{1}),...
       strcat('3cv4-TFP (dddW)__',expt_date_list_v4{2}),...
       strcat('3cv3-TFP (dmdA)__',expt_date_list_v3{2}),...
       strcat('3cv4-TFP (dddW)__',expt_date_list_v4{3})};
legend({leg{:}},'FontSize',13)

suptitle({'TFP: comparison between 2 strains; 1 mM glucose condition (blue) and 1 mM DMSP condition (green); each cell normalized to its own YFP before averaging',...
          strcat('error bars = SEM; YFP threshold = ',num2str(Ythresh)),...
          'o = 3cv3 TFP;  * = 3cv4 TFP'})
xlabel('time since treatment (hours)');
ylabel('fluorescence intensity (YFP normalized)');
set(gca, 'FontSize', 20); xlim([0 24]);
print(fluo_time_TFP,strcat(saveplace,'/TFP_glucose_3cv3_vs_3cv4'),'-dpng','-r300'); % resolution 300



% RFP    
fluo_time_RFP = figure('units','normalized','outerposition',[0 0 1 1]);
marker_s = {'m*-','m*--','m*:'}
marker_o = {'mo-','mo--','mo:'}

for i = 1:3

errorbar(aggdata_v4.time{i}(:,q),...
         aggdata_v4.Rm{i}(:,q),...
         aggdata_v4.Rsem{i}(:,q),...
         marker_s{i},'MarkerSize',15,'LineWidth',2);
hold on;

errorbar(aggdata_v3.time{i}(:,q),...
         aggdata_v3.Rm{i}(:,q),...
         aggdata_v3.Rsem{i}(:,q),...
         marker_o{i},'MarkerSize',10,'LineWidth',2);
     
end
leg = {strcat('3cv4-RFP (dmdA)__',expt_date_list_v4{1}),...
       strcat('3cv3-RFP (dddW)__',expt_date_list_v3{1}),...
       strcat('3cv4-RFP (dmdA)__',expt_date_list_v4{2}),...
       strcat('3cv3-RFP (dddW)__',expt_date_list_v3{2}),...
       strcat('3cv4-RFP (dmdA)__',expt_date_list_v4{3})};
h = legend({leg{:}},'FontSize',13)
     
suptitle({'RFP: comparison between 2 strains; 1 mM glucose condition (red) & 1 mM DMSP condition (magenta); each cell normalized to its own YFP before averaging',...
          strcat('error bars = SEM; YFP threshold = ',num2str(Ythresh)),...
          'o = 3cv3 RFP;  * = 3cv4 RFP'})
xlabel('time since treatment (hours)');
ylabel('fluorescence intensity (YFP normalized)');
set(gca, 'FontSize', 20); xlim([0 24]); ylim([-0.01 0.1]);
print(fluo_time_RFP,strcat(saveplace,'/RFP_glucose_3cv3_vs_3cv4_3'),'-dpng','-r300'); % resolution 300



%% statistical test at each time point
    % at each time point in 1 mM glucose condition, compare:
        % 1) TFP of 2 strains
        % 2) RFP of 2 strains
        
q = 1; % glucose condition      
h_R = []; h_T = [];
q_R = []; q_T = [];

for i = 1:3; 
    for t = 1:32;
    
%     x = fp_fintcells_q_leakcorr_ypos_ynorm_s.v3_Rc{t,q};
%     y = fp_fintcells_q_leakcorr_ypos_ynorm_s.v4_Rc{t,q};
    
%    x = aggdata_v3.Tc(i){t,q};
%    y = aggdata_v4.Tc(i){t,q};
    
    % get rid of zeros
    x = nonzeros(x);
    y = nonzeros(y);
        
    [h,p] = ttest2(x,y,'Alpha',0.05);% set the alpha  (0.05 is default)
    
%     h_R(t) = h;
%     p_R(t) = p;
    h_T(t) = h;
    p_T(t) = p;
        
    clearvars x y h p
    end
end
clearvars q t x y h p

save('SameColor_3cv3_vs_3cv4_gluc_allt.mat')

%% two-sample t-test #2
% for each color, compare with the OFF condition (PA-YFP of each day, after spectral correction)
    % box plot 1) dmdA: RFP in 3cv4 and TFP in 3cv3
    % box plot 2) dddW: TFP in 3cv4 and RFP in 3cv3
   
    
    
% SPECTRAL LEAKAGE CORRECTION    
 % measured signal = (matrix of slopes) x (actual signal)
% Y = b1*X  --> solve for X
% X = Y*b1 or X = b1\Y
% TFP should be the most corrected channel (YFP leaks into TFP a lot)

allc_fintcells = {}; % initialize
q = 2; % for PA-YFP cells only
c = length(tfp_fintcells_q_FPpos{q});
for cell = 1:c;

    % PA-YFP cells only
    allc_fintcells{1,cell}(1,1) = tfp_fintcells_q_FPpos{q}(1,cell);
    allc_fintcells{1,cell}(2,1) = yfp_fintcells_q_FPpos{q}(1,cell);
    allc_fintcells{1,cell}(3,1) = rfp_fintcells_q_FPpos{q}(1,cell);

    % solve the matrix
    leakcorrect_fintcells{1,cell} = b1 \ allc_fintcells{1,cell}; 

    % separate
    tfp_leakcorr_fintcells(1,cell) = leakcorrect_fintcells{1,cell}(1,1);
    yfp_leakcorr_fintcells(1,cell) = leakcorrect_fintcells{1,cell}(2,1);
    rfp_leakcorr_fintcells(1,cell) = leakcorrect_fintcells{1,cell}(3,1);

    % difference between pre-correction and post-correction
    tfp_leakcorrdiff_fintcells{1,cell} = leakcorrect_fintcells{1,cell}(1,1) - allc_fintcells{1,cell}(1,1);
    yfp_leakcorrdiff_fintcells{1,cell} = leakcorrect_fintcells{1,cell}(2,1) - allc_fintcells{1,cell}(2,1);
    rfp_leakcorrdiff_fintcells{1,cell} = leakcorrect_fintcells{1,cell}(3,1) - allc_fintcells{1,cell}(3,1);

end

save('FPctrl_LeakCorrected_YFP_2018-02-25.mat')
        
% visualize the amount of correction with histogram, separated by color and
% for all experimental conditions
close all;
SaveFig = 0; % 1 = save figure
cmap = hsv(3);

for FPcolors = 1:3;
    if FPcolors == 1;
        data = cell2mat(rfp_leakcorrdiff_fintcells);
        xaxislim = [-1.8 0];
        color = 'RFP';
    elseif FPcolors == 2;
        data = cell2mat(yfp_leakcorrdiff_fintcells);
      %  xaxislim = [-0.8 0.1];
        color = 'YFP';
    elseif FPcolors == 3;
        data = cell2mat(tfp_leakcorrdiff_fintcells);
     %   xaxislim = [-80 5]; % should be biggest correction
        color = 'TFP';
    end
    
    fig_hist(FPcolors) = figure('units','normalized','outerposition',[0 0 1 1]);

    histdata_q(FPcolors) = histogram(data,'FaceColor',cmap(FPcolors,:)); hold on;
      %      centers = histdata_q(q).BinEdges + histdata_q(q).BinWidth/2; % center of each bar in histogram
      %      heights = [histdata_q(q).Values,0];        
     %           hold on
       %     line_q(q) = plot(centers,heights,'LineWidth',3,'Color',cmap(q,:));    

  %  end
    title(strcat('magnitude of spectral leakage correction (',color,'); final value - original'),'FontSize',20)
    xlabel('magnitude of spectral leakage correction (final value minus original)','FontSize',20);
    ylabel('frequency','FontSize',20)
    set(histdata_q(:),'facealpha',0,'edgecolor','none'); % turn OFF bars
    %legend(line_q(:),{expt_conditions_string{:}},'FontSize',16,'Location','NorthWest'); 
   % xlim(xaxislim);
    % save figure
    if SaveFig == 1;
        saveas(fig_hist(FPcolors),strcat(saveplace,'\',color,'_LeakCorrectionDeltaHistogram.png'));
        close all;
    end
end

clearvars cmap FPcolors SaveFig data xaxislim color fig_hist q 
clearvars ans histdata_q centers heights line_q

%% box plots
% for each color, compare with the OFF condition (PA-YFP of each day, after spectral correction)
    % box plot 1) dmdA: RFP in 3cv4 and TFP in 3cv3
    % box plot 2) dddW: TFP in 3cv4 and RFP in 3cv3
% On each box, the central mark indicates the median, and the bottom and top 
% edges of the box indicate the 25th and 75th percentiles, respectively. 
% The whiskers extend to the most extreme data points not considered outliers, and the 
% outliers are plotted individually using the '+' symbol.
wrkdir = uigetdir; cd(wrkdir);




% dmdA
fig = figure('units','normalized','outerposition',[0 0 1 1]);
clearvars x1 x2 x3 x4 x5 x6 x g

data_source_v3 = aggdata_v3; % indicate appropriate structure
data_source_v4 = aggdata_v4;

data_source_off_v3 = aggdata_PAyfp_v3; % indicate appropriate structure
data_source_off_v4 = aggdata_PAyfp_v4;

q = 1;
t = 1;

% TFP OFF aggregated (3cv3 and 3cv4)
x1_t = [data_source_off_v3.Tc_Ynorm{:} data_source_off_v4.Tc_Ynorm{:}];
    x1_t = x1_t(~isnan(x1_t))'; x1_t = x1_t(x1_t~=0); % get rid of nan and 0
    x1_t = x1_t(x1_t < 0.1); x1_t = x1_t(x1_t > -0.1);

% RFP OFF aggregated (3cv3 and 3cv4)
x1_r = [data_source_off_v3.Rc_Ynorm{:} data_source_off_v4.Rc_Ynorm{:}];
    x1_r = x1_r(~isnan(x1_r))'; x1_r = x1_r(x1_r~=0); % get rid of nan and 0
    x1_r = x1_r(x1_r < 0.1); x1_r = x1_r(x1_r > -0.1);
    
% 3cv3 RFP #1
x2 = [data_source_v3.Rc{1}{t,q}];
    x2 = x2(~isnan(x2))'; x2 = x2(x2~=0);

% 3cv3 RFP #2
x3 = [data_source_v3.Rc{2}{t,q}];
    x3 = x3(~isnan(x3))';   x3 = x3(x3~=0); 

% 3cv4 RFP #1
x4 = [data_source_v4.Rc{1}{t,q}];
    x4 = x4(~isnan(x4))';    x4 = x4(x4~=0);

% 3cv4 RFP #2
x5 = [data_source_v4.Rc{2}{t,q}];
    x5 = x5(~isnan(x5))';   x5 = x5(x5~=0);

% 3cv4 RFP #3
x6 = [data_source_v4.Rc{3}{t,q}];
    x6 = x6(~isnan(x6))';    x6 = x6(x6~=0);
    
    
% 3cv3 TFP #1
x2_t = [data_source_v3.Tc{1}{t,q}];
    x2_t = x2_t(~isnan(x2_t))'; x2_t = x2_t(x2_t~=0);

% 3cv3 TFP #2
x3_t = [data_source_v3.Tc{2}{t,q}];
    x3_t = x3_t(~isnan(x3_t))';   x3_t = x3_t(x3_t~=0); 

% 3cv4 TFP #1
x4_t = [data_source_v4.Tc{1}{t,q}];
    x4_t = x4_t(~isnan(x4_t))';    x4_t = x4_t(x4_t~=0);

% 3cv4 TFP #2
x5_t = [data_source_v4.Tc{2}{t,q}];
    x5_t = x5_t(~isnan(x5_t))';   xx5_t5 = x5_t(x5_t~=0);

% 3cv4 TFP #3
x6_t = [data_source_v4.Tc{3}{t,q}];
    x6_t = x6_t(~isnan(x6_t))';    x6_t = x6_t(x6_t~=0);

    
    
    
    
    
% pool data
x_2_r = [x2; x3];
x_3_r = [x4; x5; x6];
x_2_t = [x2_t; x3_t];
x_3_t = [x4_t; x5_t; x6_t];
x = [x1_r;x_2_r;x_3_r; x1_t; x_2_t; x_3_t];
figure;
g = [zeros(length(x1_r), 1); ones(length(x_2_r), 1); 2*ones(length(x_3_r), 1); 3*ones(length(x1_t), 1); 4*ones(length(x_2_t), 1); 5*ones(length(x_3_t), 1)];
boxplot(x,g,'Notch','on','Labels',{'RFP (OFF) aggregated','3cv3 RFP (dddW) aggregated',...
                                   '3cv4 RFP (dmdA) aggregated',...
                                   'TFP (OFF) aggregated','3cv3 TFP (dmdA) aggregated',...
                                   '3cv4 TFP (dddW) aggregated'});

    
     
    
figure;
x = [x1;x2;x3;x4;x5;x6];
g = [zeros(length(x1), 1); ones(length(x2), 1); 2*ones(length(x3), 1); 3*ones(length(x4), 1); 4*ones(length(x5), 1); 5*ones(length(x6), 1)];
boxplot(x,g,'Notch','on','Labels',{'TFP (OFF) aggregated','3cv3 TFP (dmdA) 2018-02-25','3cv3 TFP (dmdA) 2018-03-02',...
                                   '3cv4 TFP (dddW) 2018-02-28','3cv4 TFP (dddW) 2018-03-01',...
                                   '3cv4 TFP (dddW) 2018-03-03'});
ylim([-20000 20000]);
suptitle({'TFP & RFP in PA-YFP (OFF), 3cv3 and 3cv4, 1mM glucose condition, t1, each cell divided by its own YFP signal',...
          'TFP (OFF) and RFP (OFF) cutoff at -0.1 and 0.1'});


% 
%         'RFP (3cv4) = dmdA signal from 2018-02-28 experiment, t1 glucose condition only (YFP thresh = 50; YFP-normalized)',...
%         'TFP (OFF) = PA-YFP from 2018-02-25 experiment in TFP channel, each cell divided by its own YFP signal',...
%      'TFP (3cv3) = dmdA signal from 2018-02-25 experiment, t1 glucose condition only (YFP thresh = 50; YFP-normalized)',...
%      'COMPARISON OF dmdA SIGNAL WITH OFF SIGNAL (line = median)'});
ylabel('fluorescence intensity (YFP-normalized)','FontSize',20);
set(gca,'FontSize',17)

[h,p] = ttest2(x_2_t,x_3_t,'Alpha',0.01)


%% aggregate 





 %%   
q_range = 1:9;
% h_samecolor_diffstrains = [];
h_3cv3_diffq_R = ones(length(q_range)); % n by n matrix
h_3cv3_diffq_T = ones(length(q_range)); % n by n matrix
h_3cv4_diffq_R = ones(length(q_range)); % n by n matrix
h_3cv4_diffq_T = ones(length(q_range)); % n by n matrix

add = 0;

for t = 1;
    for q = q_range
        for m = q_range;
        
            if q == m | q < m; 
                continue
            end
            
     %   x = fp_fintcells_q_leakcorr_ypos_ynorm_s.v3_T{t,q};
      %  y = fp_fintcells_q_leakcorr_ypos_ynorm_s.v4_T{t,q};
    %    x = fp_fintcells_q_leakcorr_ypos_ynorm_s.v3_R{t,q};
     %   y = fp_fintcells_q_leakcorr_ypos_ynorm_s.v4_R{t,q};
     
           x = fp_fintcells_q_leakcorr_ypos_ynorm_s.v3_R{t,q};
       y = fp_fintcells_q_leakcorr_ypos_ynorm_s.v3_R{t,m};
     
        % get rid of zeros
        x = nonzeros(x);
        y = nonzeros(y);
        
        [h,p] = ttest2(x,y,'Alpha',0.005);% set the alpha  (0.05 is default)
        % h = 1 if rejects the null hypothesis at default 5% significance level (i.e. if 1, the 2 samples are significantly different)

    %    h_samecolor_diffstrains(q+add) = h;
    %   p_samecolor_diffstrains(q+add) = p;
    

        
        h_3cv3_diffq_R(q,m) = h;
        p_3cv3_diffq_R(q,m) = p;
       
%        h_3cv3_diffq_T(q,m) = 
        
        clearvars h p x y
        
        end 
    end
end
clearvars t q

save('SameColor_3cv3_vs_3cv4.mat','h_samecolor_diffstrains','p_samecolor_diffstrains');


%% normalize [DMSP] vs fluorescence data

data_temp = aggdata_v4.Tm; 
data_subtracted = {}; data_normalized = {}; % initialize
for i = 1:length(data_temp);
   % large_average = mean(data_temp{i}((end-5:end),3)); % 1mM DMSP condition, last 5 time points
    for q = 1:9;     
     %   data_subtracted{i}(:,q) = data_temp{i}(:,q) - data_temp{i}(1,q); % subtract t1 of that q from each data point of that q
        data_normalized{i}(:,q) = data_temp{i}(:,q) ./ large_average;
    end    
end

% store in permanent variable
if isequal([data_temp{1}],[aggdata_v4.Rm{1}]);
    aggdata_v4.Rm_sub = data_subtracted;
    aggdata_v4.Rm_norm = data_normalized;
elseif isequal([data_temp{1}],[aggdata_v4.Tm{1}]); 
    aggdata_v4.Tm_sub = data_subtracted;
    aggdata_v4.Tm_norm = data_normalized;
elseif isequal([data_temp{1}],[aggdata_v3.Rm{1}]); 
    aggdata_v3.Rm_sub = data_subtracted;
    aggdata_v3.Rm_norm = data_normalized;
elseif isequal([data_temp{1}],[aggdata_v3.Tm{1}]); 
    aggdata_v3.Tm_sub = data_subtracted;
    aggdata_v3.Tm_norm = data_normalized;
end
clearvars data_temp data_subtracted data_normalized large_average i q

% %% [DMSP] fluo vs. time for all experiments (copied to above)
% % normalize
% 
% plot_dmdA_v3 = 1;
% plot_dddW_v3 = 0;
% plot_dmdA_v4 = 1;
% plot_dddW_v4 = 0;
% 
% notnorm = 0; % fp_fintview_q_leakcorr_ypos_ynorm
% sub = 0; % subtracted
% norm = 1; % subtracted and normalized
% 
% % manual 
% title_text = ['3cv3-dmdA (TFP) and 3cv4-dmdA (RFP) of all experiments; each data set subtracted by its own t1 value',...
%             'each data set subtracted by its own t1 value; then, all values in each experiment is normalized by average fluo intensity of last 5 time points of the 1 mM DMSP condition'];
% 
% 
% marker_s = {'*-','*--','*:'};
% marker_o = {'o-','o--','o:'};
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% cmap = hsv(9);
% for i = 1:3
%     for q = 1:9;
%         
%         if plot_dmdA_v4 == 1;
%             time_data = aggdata_v4.time{i}(:,1);
%             if notnorm == 1;
%                 plot_data = aggdata_v4.Rm{i}(:,q);
%             elseif sub == 1;               
%                 plot_data = aggdata_v4.Rm_sub{i}(:,q);
%             elseif norm == 1;
%                 plot_data = aggdata_v4.Rm_norm{i}(:,q);
%             end
%             plot(time_data,plot_data,marker_s{i},'MarkerSize',15,'Color',cmap(q,:),'LineWidth',2); hold on;               
%         end
%         
%         if plot_dddW_v4 == 1;
%             time_data = aggdata_v4.time{i}(:,1);
%             if notnorm == 1;
%                 plot_data = aggdata_v4.Tm{i}(:,q);
%             elseif sub == 1;  
%                 plot_data = aggdata_v4.Tm_sub{i}(:,q);
%             elseif norm == 1;
%                 plot_data = aggdata_v4.Tm_norm{i}(:,q);
%          end
%             plot(time_data,plot_data,marker_s{i},'MarkerSize',15,'Color',cmap(q,:),'LineWidth',2); hold on;          
%         end        
%         
%         if plot_dmdA_v3 == 1;
%             if i == 3;
%                 continue
%             end
%             time_data = aggdata_v3.time{i}(:,1);
%             if notnorm == 1;
%                  plot_data = aggdata_v3.Tm{i}(:,q);
%             elseif sub == 1;  
%                  plot_data = aggdata_v3.Tm_sub{i}(:,q);
%             elseif norm == 1;
%                  plot_data = aggdata_v3.Tm_norm{i}(:,q);
%            end
%             plot(time_data,plot_data,marker_o{i},'MarkerSize',15,'Color',cmap(q,:),'LineWidth',2); hold on;       
%         end        
%         
%         if plot_dddW_v3 == 1;
%             if i == 3;
%                 continue
%             end
%             time_data = aggdata_v3.time{i}(:,1);
%             if notnorm == 1;
%                 plot_data = aggdata_v3.Rm{i}(:,q);
%             elseif sub == 1;  
%                 plot_data = aggdata_v3.Rm_sub{i}(:,q);
%             elseif norm == 1;
%                 plot_data = aggdata_v3.Rm_norm{i}(:,q);
%                end
%             plot(time_data,plot_data,marker_o{i},'MarkerSize',15,'Color',cmap(q,:),'LineWidth',2); hold on;            
%         end
%         
%     end     
% end
% 
% suptitle({title_text,'o = 3cv3;   * = 3cv4','glucose, DMSP 1mM, DMSP 250 uM, DMSP 75 uM, DMSP 1 uM conditions only'})
% xlabel('time since treatment (hours)');
% ylabel('fluorescence intensity, YFP-normalized');
% legend({'glucose 1 mM',...
%        'DMSP 1 mM',...
%        'DMSP 250 \muM',...
%        'DMSP 100 \muM',...
%        'DMSP 75 \muM',...
%        'DMSP 50 \muM',...
%        'DMSP 25 \muM',...
%        'DMSP 10 \muM',...
%        'DMSP 1 \muM'},...
%        'location','northwest');
% set(gca,'FontSize',20);
% xlim([0 24])
% %%
% print(strcat(saveplace,'/dmdA_3cv3_3cv4_normalized_sparse'),'-dpng','-r300'); % resolution 300
% 
%         


%% derivatives
marker_s = {'*-','*--','*:'}; % 3cv4
marker_o = {'o-','o--','o:'}; % 3cv3

%figure
for i = 1:2;
    for q = [1 2 3 5 7];
        x = 1:32;%aggdata_v3.time{i}(:,1);
       % y = smooth(aggdata_v3.Tm_norm{i}(:,q)); % smoothed raw data
        y = aggdata_Am(:,q); 
      %  y = aggdata_v3.Tm_norm{i}(:,q);
        dydx = diff(y(:))./diff(x(:));

        %ddydx = diff(dydx(:))./diff(x(1:31));
        plot(x(2:end),dydx,marker_o{i},'MarkerSize',10,'Color',cmap(q,:),'LineWidth',2); hold on;
    %  plot(x,y,'o-');hold on;
    end
end
leg = {expt_conditions_string{[1 2 3 5 7]}}
legend(leg,'FontSize',20)
suptitle({'o = dddW (mean signal of all expts) derivatives, after smoothing raw data',...
       '* = dmdA (mean signal of all expts) derivatives, after smoothing raw data'})
ylabel('derivative')
xlabel('time point (every 45 minutes)')
set(gca,'FontSize',20)



%% make a plot with averaged normalized signal vs. 
% average across all expts (3cv3 and 3cv4)
% copied from PlotGenerator_4.m
avgFluo_shaded = figure('units','normalized','outerposition',[0 0 1 1]);

q_range = 1:9;

% take the average at each time point
% dmdA
for i = 1:3;
    for t = 1:32;
        for q = 1:9
                      
            if t <= 22
            data = [aggdata_v3.Tm_norm{1}(t,q),aggdata_v3.Tm_norm{2}(t,q),...
                    aggdata_v4.Rm_norm{1}(t,q),aggdata_v4.Rm_norm{2}(t,q),aggdata_v4.Rm_norm{3}(t,q)];
            elseif t > 22;
                  data = [aggdata_v3.Tm_norm{1}(t,q),...
                    aggdata_v4.Rm_norm{1}(t,q),aggdata_v4.Rm_norm{2}(t,q),aggdata_v4.Rm_norm{3}(t,q)];
             end
                             
            aggdata_Am(t,q) = mean(data);
            aggdata_Astd(t,q) = std(data);

        end
    end
end

% dddW
for i = 1:3;
    for t = 1:32;
        for q = 1:9
            
            
            if t <= 22
            data = [aggdata_v3.Rm_norm{1}(t,q),aggdata_v3.Rm_norm{2}(t,q),...
                    aggdata_v4.Tm_norm{1}(t,q),aggdata_v4.Tm_norm{2}(t,q),aggdata_v4.Tm_norm{3}(t,q)];
            elseif t > 22;
                  data = [aggdata_v3.Rm_norm{1}(t,q),...
                    aggdata_v4.Tm_norm{1}(t,q),aggdata_v4.Tm_norm{2}(t,q),aggdata_v4.Tm_norm{3}(t,q)];
             end
             
                
            aggdata_Wm(t,q) = mean(data);
            aggdata_Wstd(t,q) = std(data);

        end
    end
end

plot_rfp =aggdata_Wm;
plot_tfp = 0% tfp_mean_q_ypos_ynorm_t;

error_rfp = aggdata_Wstd;
error_tfp =0% tfp_sem_q_ypos_ynorm_t;

cmap = hsv(9);
%figure;
for q = 1:9;
y = plot_rfp(:,q);
x = [1:32]';
dy = error_rfp(:,q);  % made-up error values
%fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.9 .9 .9],'color',cmap(q,:)); hold on;
errorbar(x,y,dy,'LineWidth',2,'Color',cmap(q,:)); hold on;
%line(x,y)
end
leg = {'1 mM glucose','1 mM DMSP','250 \muM DMSP','100 \muM DMSP',...
    '75 \muM DMSP','50 \muM DMSP','25 \muM DMSP','10 \muM DMSP','1 \muM DMSP'}
legend({leg{:}},'FontSize',16)
suptitle({'aggregated gene expression signal over time; includes all expts, 3cv3 and 3cv4 normalized by YFP, normalized by average of last 5 time points of 250 \muM DMSP condition',...
    'one expt ends at t22'});
xlabel('time point')
ylabel('gene expression signal (RFP or TFP of each cell normalized by its own YFP)')
set(gca,'FontSize',20);

% the real plot
for q = 1:9;
  %  h(q) = subplot(1,9,q); hold on;
    
    % specify data
    x = time_expt_hours(:,1);
    yr = plot_rfp(:,q);
    dyr = error_rfp(:,q);    
 %   yt = plot_tfp(:,q);
  %  dyt = error_tfp(:,q);
    
hold on
    % shaded errorbars
        % flipud -> can make continuous line that encapsulates a shape 
        % in this case, shape = -ve and +ve errorbars to be filled
    %fill([x;flipud(x)],[yr-dyr;flipud(yr+dyr)],[.9 .9 .9],'linestyle','none');
    % plot
    errorbar(x,yr,dyr,...
        'o-','color',cmap(q,:))
  
    
        
%     errorbar(x,yt,dyt,...
%          marker_tfp,'LineWidth', linewidth_tfp, ...
%          'Color', color_tfp_in, ...
%          'MarkerEdgeColor', color_tfp_out, ...
%          'MarkerFaceColor',color_tfp_in, ...
%          'MarkerSize', markersize_tfp)

    % fill([x;flipud(x)],[yt-dyt;flipud(yt+dyt)],[.9 .9 .9],'linestyle','none');

  
   % axis([0 10.02 0 18]);
  %  title(expt_conditions_string{q});
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

