% This program is free to use for acedemic purposes only
% Users must refer to AiNU manual for disclaim and copyright notice.
% Author: Hualin Zhan, Australian National University

clear;

img_name = 'kernel_cmpr.pdf';
g_set = struct('project_dir','../data/output/'...
              ,'sample_name','sim1_'...
              ,'methods',{'trpl'}... 
              ,'m_units',{''}... 
              ,'m_values',{{'22','1960'}}... % {'22'},the size of this cell determines the number of figures
              );

mat_head = 'ainu_result_';
mat_end = {'exp','sqexp','rq','mt32','mt52'};
m_dates = {'240430','240430','240430','240430','230724'};
kn_names = {'ARD Exponential','ARD Squared Exponential','ARD Rational Quadratic','ARD Matern 3/2','ARD Matern 5/2'};

n_mruns = 22; % number of runs for each same settings/file
n_steps = 20; % number of convergence steps in each run
mat_para = 't';
c_fun = 5;
f_type = 6;
data_for_excel_f = true;
color_shade = {[1.0000    0.8314    0.7608], [0.96,0.99,1], [0.9608    0.9020    1.0000], [1.0000    0.9529    0.8392], [0.9255    1.0000    0.8314]};
color_edge = {[1.0000    0.4118    0.1608], [0.16,0.67,0.92], [0.8235    0.5412    1.0000], [0.9686    0.7686    0.3020], [0.5725    0.8118    0.2588]};
color_mean = {[0.8510    0.3255    0.0980], [0,0.4471,0.7412], [0.7176    0.2745    1.0000], [0.8784    0.6235    0.0275], [0.2588    0.4314    0.0392]};

m_methods_all = strjoin({g_set(1).methods},'_');
n_fluences = size(g_set,2);
n_kns = size(mat_end,2);
mat_names = cell(n_fluences*n_kns,n_mruns);
cost_values_all = cell(n_fluences*n_kns,1); 
cost_values_mean = cost_values_all; 
cost_values_std = cost_values_all;
cost_values_max = cost_values_all;
cost_values_min = cost_values_all;
cost_values_pos = cost_values_all;
cost_values_neg = cost_values_all;
data_for_excel = cell(1,n_fluences);

%% put all mat file names in mat_names

for i_kn = 1:1:n_kns
    m_values_cell = cell(1,n_fluences);
    for i_fluence = 1:1:n_fluences
        m_values_cell{i_fluence} = strjoin(g_set(i_fluence).m_values,'_');
        for i_run_c = 1:1:n_mruns
            if strcmp(mat_end{i_kn},'mt52')
                if strcmp(m_values_cell{i_fluence}, '22')
                    i_file = 3*c_fun + 1;
                elseif strcmp(m_values_cell{i_fluence}, '1960')
                    i_file = 3*c_fun + 2;
                elseif strcmp(m_values_cell{i_fluence}, '22_1960')
                    i_file = 3*c_fun + 3;
                end
                n_file = num2str(i_file,'%02d');
                
                mat_names{(i_kn-1)*n_fluences+i_fluence,i_run_c} = [g_set(1).project_dir,mat_head,g_set(1).sample_name,m_methods_all,'_',m_values_cell{i_fluence},'_',g_set(1).m_units,m_dates{i_kn},'_',n_file,'_',mat_para,'_',num2str(i_run_c,'%02d'),'.mat'];
            else
                mat_names{(i_kn-1)*n_fluences+i_fluence,i_run_c} = [g_set(1).project_dir,mat_head,g_set(1).sample_name,m_methods_all,'_',m_values_cell{i_fluence},'_',g_set(1).m_units,m_dates{i_kn},'_',num2str(i_run_c,'%02d'),'_c',num2str(c_fun),'_kernel_',mat_end{i_kn},'.mat'];
            end
        end
    end
end

%% start plots

load(mat_names{1,1});
t_uni = exp_data{1}.t_intp_trim{1};

for i_fluence = 1:1:n_fluences % i_fluence is also the index of the plot
    if data_for_excel_f
        data_for_excel{i_fluence} = (1:1:n_steps)';
    end
    kn_fig = figure('visible',false);
    kn_fig.Units = "centimeters";
    kn_fig.Position = [3,3,8.3,9.3/1.618];
    for i_kn = 1:1:n_kns
        i_row_in_mat = (i_kn-1)*n_fluences+i_fluence;
        cost_values_all{i_row_in_mat,1} = nan(n_mruns,n_steps); 
        for i_run_c = 1:1:n_mruns
            if isfile(mat_names{i_row_in_mat,i_run_c})
                [steps, cost_r] = calc_convergence_cost(mat_names{i_row_in_mat,i_run_c},t_uni,false,f_type);
                if isrow(cost_r)
                    cost_values_all{i_row_in_mat,1}(i_run_c,:) = cost_r;
                else
                    cost_values_all{i_row_in_mat,1}(i_run_c,:) = cost_r';
                end
            else
                disp([mat_names{i_row_in_mat,i_run_c},' does not exist!']);
                continue
            end
        end
        cost_values_mean{i_row_in_mat,1} = mean(cost_values_all{i_row_in_mat,1}, 'omitnan');
        cost_values_std{i_row_in_mat,1} = std(cost_values_all{i_row_in_mat,1}, 'omitnan');
        cost_values_max{i_row_in_mat,1} = max(cost_values_all{i_row_in_mat,1});
        cost_values_min{i_row_in_mat,1} = min(cost_values_all{i_row_in_mat,1});
        cost_values_pos{i_row_in_mat,1} = cost_values_mean{i_row_in_mat,1}+ 0.5*(cost_values_max{i_row_in_mat,1}-cost_values_mean{i_row_in_mat,1});
        cost_values_neg{i_row_in_mat,1} = cost_values_mean{i_row_in_mat,1}- 0.5*(cost_values_mean{i_row_in_mat,1}-cost_values_min{i_row_in_mat,1});
        steps = steps';

        if data_for_excel_f
            data_for_excel{i_fluence} = [data_for_excel{i_fluence} cost_values_mean{i_row_in_mat,1}' cost_values_std{i_row_in_mat,1}' cost_values_max{i_row_in_mat,1}' cost_values_min{i_row_in_mat,1}'];
        end

        if ~all(isnan(cost_values_all{i_row_in_mat,1}(:)))
            semilogy(steps,cost_values_neg{i_row_in_mat,1},'-','Color',color_edge{i_kn},'HandleVisibility','off','LineWidth',0.75);
            hold on
            semilogy(steps,cost_values_pos{i_row_in_mat,1},'-','Color',color_edge{i_kn},'HandleVisibility','off','LineWidth',0.75);
            yr_area = patch([steps fliplr(steps)], [cost_values_neg{i_row_in_mat,1} fliplr(cost_values_pos{i_row_in_mat,1})],color_shade{i_kn},'HandleVisibility','off');
            yr_area.EdgeColor = 'none';
            semilogy(steps,cost_values_mean{i_row_in_mat,1},'-','Color',color_mean{i_kn},'DisplayName',kn_names{i_kn},'LineWidth',1);
            lgd = legend;
        end
    end
    xlim([1 n_steps]);
    xlabel('Steps');
    ylabel('RMSE');
    set(gca,'FontSize',9);
    lgd.Location = 'northeast';
    lgd.FontSize = 5;
    hold off
    exportgraphics(kn_fig,['../tif/',img_name],'ContentType','vector');
end