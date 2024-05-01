% This program is free to use for acedemic purposes only
% Users must refer to AiNU manual for disclaim and copyright notice.
% Author: Hualin Zhan, Australian National University

clear;
datafilename = 'ainu_result_p0h_trpl_22_1960_240148_09_asEnmmBC_';

scathist_pdf(['../data/output/',datafilename,'.mat'],3,2,'off',['../tif/',datafilename,'main_mm.pdf']);

close all;

%% the following function is explicitly included here for the convenience in customization
function scathist_pdf(mat_file,x_ind,y_ind,visible_flag,save_path)
    
    load(mat_file);

    % the following parameters can be adjusted for customized visualization, e.g., different color scale
    cost_prec = 2; % decimal points shown in the clorbar
    x_fit_range = [10 20]; % x-range (power of 10) for linear fit to t&s
    curve_shape_x = 0.05; % change for bandwidth value for the yellow curves, smaller value = curve becomes more discrete
    curve_shape_y = 0.05;
    r_fit_lin = 0.3; % fot Nt and sigma fitting only: z-range (percentage between z_low and z_high), recommend between 0.05 and 0.2
    top_dev_fixed = 10;
    scale_ratio_linear = 1.05;
    scale_bound = -1;
    scale_ratio_log2 = 0;
    
    cost_recalc = 0; % cost recalculation flag
    c_matrix = struct('ratio_rms',[1, 0, 0, 1]...
        ,'norm_flag',true...               % flag to normalize the caculated data using the first data point, set true for trpl
        ,'intp_scale','linear'...
        ,'intp_method','linear'...
        );

    n_tick_x = 3;
    n_tick_y = 3;
    
    all_data = sortrows(all_data,"cost_bayes","ascend");
    x_2v_name = all_data.Properties.VariableNames{end-x_ind};
    y_2v_name = all_data.Properties.VariableNames{end-y_ind};
    x_2v = all_data.(x_2v_name);
    y_2v = all_data.(y_2v_name);

    % re-scale Et_r to the range of [0, 0.5]
    n_tests = size(all_data,1);
    new_Et = zeros(n_tests,1);
    if strcmp(x_2v_name,'Et_r')
        new_Et(x_2v<=0.5) = x_2v(x_2v<=0.5);
        new_Et(x_2v>0.5) = 1 - x_2v(x_2v>0.5);
        x_2v = new_Et;
    elseif strcmp(y_2v_name,'Et_r')
        new_Et(y_2v<=0.5) = y_2v(y_2v<=0.5);
        new_Et(y_2v>0.5) = 1 - y_2v(y_2v>0.5);
        y_2v = new_Et;
    end

    %% convert categorical to numerical

    if isa(x_2v,'categorical')
        [x_dict,~,x_2v]=unique(x_2v); 
    end
    if isa(y_2v,'categorical')
        [y_dict,~,y_2v]=unique(y_2v); 
    end
    
    %% cost recalculation
    
    if cost_recalc
        
        % 1, set the universal time axis and the cost matrix
        t_uni_file = '../data/output/t_universal.mat'; %#ok<UNRCH> 
        load(t_uni_file);
    
        c_results = cell(1,4);
        c_results{1,1} = all_data;
        c_results{1,2} = exp_data;
    
        n_input = size(c_results{1,2}{1}.in_target,2); % the number of fluences in one fitting file
        c_results{1,3} = struct('y_target_intp',cell(n_input,1)...
            ,'y_target_intp_lg',cell(n_input,1)...
            ,'dydx_target_intp',cell(n_input,1)...
            ,'dydx_target_intp_lg',cell(n_input,1)...
            );
        run_size = size(c_results{1,1}.step_bayes,1);
        c_results{1,4} = nan(n_input,run_size);
    
        for i_input = 1:1:n_input
    
            % 2, recalculate the target corresponding to the new time - t_uni
            c_results{1,3}(i_input).y_target_intp = interp1(c_results{1,2}{1}.x_target{1,i_input},c_results{1,2}{1}.y_target{1,i_input},t_uni,'spline','extrap');
            if c_matrix.norm_flag
                c_results{1,3}(i_input).y_target_intp = c_results{1,3}(i_input).y_target_intp./c_results{1,3}(i_input).y_target_intp(1);
            end
    
            c_results{1,3}(i_input).dydx_target_intp = d_fun(c_results{1,3}(i_input).y_target_intp,1);
            c_results{1,3}(i_input).y_target_intp_lg = log10(c_results{1,3}(i_input).y_target_intp);
            c_results{1,3}(i_input).dydx_target_intp_lg = d_fun(c_results{1,3}(i_input).y_target_intp_lg,1);
    
            % 3, recalculate the cost for each run
            for i_run = 1:1:run_size
                if ~isnan(c_results{1,1}.cost_bayes(i_run))
                    c_results{1,4}(i_input,i_run) = cost_fun(c_results{1,1}.x_sim{i_run,1}{1,1}{i_input,1}...
                        ,c_results{1,1}.y_sim{i_run,1}{1,1}{i_input,1}...
                        ,t_uni...
                        ,c_results{1,3}(i_input).y_target_intp...
                        ,c_results{1,3}(i_input).dydx_target_intp...
                        ,c_results{1,3}(i_input).y_target_intp_lg...
                        ,c_results{1,3}(i_input).dydx_target_intp_lg...
                        ,c_matrix.ratio_rms...
                        ,c_matrix.norm_flag...
                        ,c_matrix.intp_scale...
                        ,c_matrix.intp_method...
                        );
                end
            end
        end
    
        z_2v = mean(c_results{1,4},1)'; % average cost for all fluences
    else
        z_2v = all_data.cost_bayes; %#ok<UNRCH> 
    end
    %% set the lower and upper bound for the plot
    
    z_low = ceil(min(z_2v)*10^cost_prec)/10^cost_prec;
    if z_low >0
        z_high = z_low *scale_ratio_linear;
    elseif z_low <= -1.5
        z_high = scale_bound;
    else
        z_high = scale_ratio_log2;
    end
    %% calculate for colored scatter and kernel distribution
    
    [xyz_table,cm,color_index,colorlabel]=scatterc3(x_2v,y_2v,z_2v,[z_low z_high],true,cost_prec);
    
    index_dist = find(z_2v>(z_low-1)&z_2v<(z_low+0.1*(z_high-z_low)));
    %% Rescale the data and set the plot scale parameter
    
    x_axis_label_linear = false;
    y_axis_label_linear = false;
    for i_var_opt = 1:1:size(all_var,2)
        if strcmp(all_var(i_var_opt).Name,x_2v_name)
            xl_s=all_var(i_var_opt).Range;
            if isa(xl_s,'cell')
                n_tick_x = size(x_dict,1);
                xl_s=1:1:n_tick_x;
                x_axis_label_linear = true;
            end
            if strcmp(all_var(i_var_opt).Transform,'none')
                x_axis_label_linear = true;
            end
        end
        if strcmp(all_var(i_var_opt).Name,y_2v_name)
            yl_s=all_var(i_var_opt).Range;
            if isa(yl_s,'cell')
                n_tick_y = size(y_dict,1);
                yl_s=1:1:n_tick_y;
                y_axis_label_linear = true;
            end
            if strcmp(all_var(i_var_opt).Transform,'none')
                y_axis_label_linear = true;
            end
        end
    end
    
    if strcmp(x_2v_name,'Etrap') || strcmp(x_2v_name,'Et_r') || x_axis_label_linear || any(xl_s==0)
        x_scale = 'linear';
        x_dist = x_2v(index_dist);
        xl_s = [0 0.5];
        xl_h = xl_s;
        tick_x = linspace(xl_h(1),xl_h(end),n_tick_x);
        format_x = '%g';
        curve_shape_x = 0.005;
    else
        x_scale = 'log';
        x_dist = log10(x_2v(index_dist));
        xl_h = log10(xl_s);
        tick_x = logspace(xl_h(1),xl_h(end),(xl_h(end)-xl_h(1)+1)); %((xl_h(end)-xl_h(1))/2+1))
        format_x = '%.1e';
    end
    
    if strcmp(y_2v_name,'Etrap') || strcmp(y_2v_name,'Et_r') || y_axis_label_linear || any(yl_s==0)
        y_scale = 'linear';
        y_dist = y_2v(index_dist);
        yl_s = [0 0.5];
        yl_h = yl_s;
        tick_y = linspace(yl_h(1),yl_h(end),n_tick_y);
        format_y = '%g';
        curve_shape_y = 0.005;
    else
        y_scale = 'log';
        y_dist = log10(y_2v(index_dist));
        yl_h = log10(yl_s);
        tick_y = logspace(yl_h(1),yl_h(end),(yl_h(end)-yl_h(1)+1)); %((yl_h(end)-yl_h(1))/2+1))
        format_y = '%.1e';
    end
    %% plot
    
    figure('visible',visible_flag);
    
    sh = tiledlayout(5,5,'TileSpacing','tight');
    sh.Units = "centimeters";
    sh.Position = [3,3,3.65,2];
    
    nexttile([1 4]); %
    [f_x,x_k] = ksdensity(x_dist,'Bandwidth',curve_shape_x); % use histfit or fitdist for more options
    area(x_k,f_x,'EdgeColor','#EDB120','FaceColor','#FFFF00','LineWidth',0.5);
    xlim(xl_h);
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    axis off;
    
    x_max_i = x_k(f_x==max(f_x));
    
    if all_data.cost_bayes(10)-all_data.cost_bayes(1)<0.001 && all_data.cost_bayes(1)>=0
        top_dev_index=find(all_data.cost_bayes>all_data.cost_bayes(1)+0.001);
        top_dev = top_dev_index(1)-1;
    elseif all_data.cost_bayes(1)<-2.5
        top_dev_index=find(all_data.cost_bayes>(all_data.cost_bayes(1)+0.04));
        top_dev = top_dev_index(1)-1;
        if top_dev <= 10
            top_dev = top_dev_fixed;
        end
    else
        top_dev = top_dev_fixed;
    end
    x_top = mean(x_2v(1:top_dev));
    x_top_s = std(x_2v(1:top_dev));

    if x_top_s >= 0.9*x_top
        x_top_max = max(x_2v(1:top_dev));
        x_top_min = min(x_2v(1:top_dev));
        if x_axis_label_linear
            x_t_max = x_top_max;
            x_t_min = x_top_min;
        else
            x_t_max = log10(x_top_max);
            x_t_min = log10(x_top_min);
        end

        if xl_h(end) ~= 0
            dx_max = abs((x_t_max-xl_h(end))/xl_h(end));
        else
            dx_max = abs((x_t_max-xl_h(end))/x_t_max);
        end
        if xl_h(1) ~= 0
            dx_min = abs((x_t_min-xl_h(1))/xl_h(1));
        else
            dx_min = abs((x_t_min-xl_h(1))/x_t_min);
        end

        if dx_max<=0.1 && dx_min>0.1
            x_txt=['\geq',sprintf('%.2g',x_top_min)];
            txt_x_position = 'right';
        elseif dx_max>0.1 && dx_min<=0.1
            x_txt=['\leq',sprintf('%.2g',x_top_max)];
            txt_x_position = 'left';
        elseif dx_max>0.1 && dx_min>0.1
            x_txt=sprintf('[%.2g, %.2g]',x_top_min,x_top_max);
            txt_x_position = 'left';
        else
            x_txt=[sprintf('%.3g',x_top),'\pm',sprintf('%.3g',x_top_s)];
            txt_x_position = 'center';
        end
    else
        x_txt=[sprintf('%.3g',x_top),'\pm',sprintf('%.3g',x_top_s)];
        txt_x_position = 'center';
    end

    text(x_max_i(1),max(f_x)/3,x_txt,'HorizontalAlignment',txt_x_position,'fontsize',7)
    
    
    nexttile([1 1]);
    cbs = colorbar('Ticks',[0,0.5,1],'TickLabels',colorlabel,'Location','west','fontsize',7);%,'Location','manual','Position',[0.81 0.8 0.04 0.1]
    cbs.Layout.Tile = 5;
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    axis off;
    
    nexttile([4 4]);
    
    scatter3(xyz_table.x,xyz_table.y,xyz_table.z,20,cm(color_index,:),'filled','MarkerEdgeColor','k','LineWidth',0.4);
    % linear fitting between t and s
    if (strcmp(y_2v_name,'Nt')||strcmp(y_2v_name,'sigma_n'))&&(strcmp(x_2v_name,'Nt')||strcmp(x_2v_name,'sigma_n'))
        index_for_fit_lin = find( z_2v>=z_low & z_2v<=(z_low+r_fit_lin*(z_high-z_low)) );
        x_by_z = log10(x_2v(index_for_fit_lin));
        y_by_z = log10(y_2v(index_for_fit_lin));
        index_x_fit = find(x_by_z>x_fit_range(1) & x_by_z<x_fit_range(2));
        x_fit_lin = x_by_z(index_x_fit);
        y_raw = y_by_z(index_x_fit);
        fit_coefficents = polyfit(x_fit_lin,y_raw,1);
        y_fit_lin = polyval(fit_coefficents,x_fit_lin);
        fit_eqn2 = string(" y="+10^fit_coefficents(2) + "x^{" + fit_coefficents(1)) + "}";

        hold on
        plot3(10.^x_fit_lin,10.^y_fit_lin,(z_low-0.5)*ones(size(x_fit_lin)),'r','LineWidth',2);
        text(min(xyz_table.x),min(xyz_table.y),z_low,fit_eqn2,'Color','r','FontWeight','bold');
        hold off
    end
    view(0,-90);
    if strcmp(x_2v_name,'ts_a')
        x_axis_name = '\it a \rm (cm^{-1})';
    elseif strcmp(x_2v_name,'Nt')
        x_axis_name = '\it N_t \rm (cm^{-3})';
    elseif strcmp(x_2v_name,'sigma_n')
        x_axis_name = '\it \sigma \rm (cm^{2})';
    elseif strcmp(x_2v_name,'Etrap') || strcmp(x_2v_name,'Et_r')
        x_axis_name = '\it E_t';
    elseif strcmp(x_2v_name,'NA_ABS')
        x_axis_name = '\it N_d \rm (cm^{-3})';
    elseif strcmp(x_2v_name,'mun_ABS')
        x_axis_name = '\it \mu_n \rm [cm^{2}/(V\cdots)]';
    elseif strcmp(x_2v_name,'mup_ABS')
        x_axis_name = '\it \mu_p \rm [cm^{2}/(V\cdots)]';
    elseif strcmp(x_2v_name,'B_rad')
        x_axis_name = '\it B_r \rm (cm^{3}/s)';
    elseif strcmp(x_2v_name,'C_Aug')
        x_axis_name = '\it C_A \rm (cm^{6}/s)';
    end
    if strcmp(y_2v_name,'ts_a')
        y_axis_name = '\it a \rm (cm^{-1})';
    elseif strcmp(y_2v_name,'Nt')
        y_axis_name = '\it N_t \rm (cm^{-3})';
    elseif strcmp(y_2v_name,'sigma_n')
        y_axis_name = '\it \sigma \rm (cm^{2})';
    elseif strcmp(y_2v_name,'Etrap') || strcmp(y_2v_name,'Et_r')
        y_axis_name = '\it E_t';
    elseif strcmp(y_2v_name,'NA_ABS')
        y_axis_name = '\it N_d \rm (cm^{-3})';
    elseif strcmp(y_2v_name,'mun_ABS')
        y_axis_name = '\it \mu_n \rm [cm^{2}/(V\cdots)]';
    elseif strcmp(y_2v_name,'mup_ABS')
        y_axis_name = '\it \mu_p \rm [cm^{2}/(V\cdots)]';
    elseif strcmp(y_2v_name,'B_rad')
        y_axis_name = '\it B_r \rm (cm^{3}/s)';
    elseif strcmp(y_2v_name,'C_Aug')
        y_axis_name = '\it C_A \rm (cm^{6}/s)';
    end
    xlabel(x_axis_name); 
    ylabel(y_axis_name); 
    set(gca,'Xscale',x_scale,'Yscale',y_scale);
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',7); 
    box on;
    xlim(xl_s);
    ylim(yl_s);
    tick_y = [1e-2 1e-1]; % change this
    tick_x = [1e0 1e1]; % change this
    xticks(tick_x);
    if isa(xl_s,'cell') %x_axis_label_linear
        xticklabels(cellstr(x_dict'));
    else
        xticklabels(tick_x);
        xtickformat(format_x);
    end
    yticks(tick_y);
    if isa(yl_s,'cell') %y_axis_label_linear
        yticklabels(cellstr(y_dict'));
    else
        yticklabels(tick_y);
        ytickformat(format_y);
    end
    grid off;
    
    nexttile([4 1]); %
    [f_y,y_k] = ksdensity(y_dist,'Bandwidth',curve_shape_y);
    area(y_k,f_y,'EdgeColor','#EDB120','FaceColor','#FFFF00','LineWidth',0.5);
    xlim(yl_h);
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    axis off;
    view(90,90);
    
    y_max_i = y_k(f_y==max(f_y));

    y_top = mean(y_2v(1:top_dev));
    y_top_s = std(y_2v(1:top_dev));

    if y_top_s >= 0.9*y_top
        y_top_max = max(y_2v(1:top_dev));
        y_top_min = min(y_2v(1:top_dev));
        if y_axis_label_linear
            y_t_max = y_top_max;
            y_t_min = y_top_min;
        else
            y_t_max = log10(y_top_max);
            y_t_min = log10(y_top_min);
        end

        if yl_h(end) ~= 0
            dy_max = abs((y_t_max-yl_h(end))/yl_h(end));
        else
            dy_max = abs((y_t_max-yl_h(end))/y_t_max);
        end
        if yl_h(1) ~= 0
            dy_min = abs((y_t_min-yl_h(1))/yl_h(1));
        else
            dy_min = abs((y_t_min-yl_h(1))/y_t_min);
        end

        if dy_max<=0.1 && dy_min>0.1
            y_txt=['\geq',sprintf('%.2g',y_top_min)];
            txt_y_position = 'right';
        elseif dy_max>0.1 && dy_min<=0.1
            y_txt=['\leq',sprintf('%.2g',y_top_max)];
            txt_y_position = 'left';
        elseif dy_max>0.1 && dy_min>0.1
            y_txt=sprintf('[%.2g, %.2g]',y_top_min,y_top_max);
            txt_y_position = 'center';
        else
            y_txt=[sprintf('%.3g',y_top),'\pm',sprintf('%.3g',y_top_s)];
            txt_y_position = 'left';
        end
    else
        y_txt=[sprintf('%.3g',y_top),'\pm',sprintf('%.3g',y_top_s)];
        txt_y_position = 'left';
    end

    txt_y_hdl = text(y_max_i(1)-0.3,max(f_y)/3,y_txt,'HorizontalAlignment',txt_y_position,'fontsize',7);
    set(txt_y_hdl,'Rotation',270);
    
    exportgraphics(sh,save_path,'ContentType','vector')

end