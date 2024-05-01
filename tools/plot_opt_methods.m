% This program is free to use for acedemic purposes only
% Users must refer to AiNU manual for disclaim and copyright notice.
% Author: Hualin Zhan, Australian National University

clear;
target_value = [5e11 1e14 0.01 1e-11];
n_mruns = 22; % number of runs for each optimization method
img_name = 'opt_cmpr.pdf';

mat_head = {'ainu_result_','comsol_gbm_','comsol_sa_'};
m_date = {'240420','240423','240422'};
g_set = struct('project_dir','../data/output/'...
              ,'sample_name','sim1_'...
              ,'methods',{'trpl'}... 
              ,'m_units',{''}... 
              ,'m_values',{{'22','1960'}}... 
              );

mat_para = 'tnmB';
m_methods_all = strjoin({g_set(1).methods},'_');
n_fluences = size(g_set,2);
n_opts = size(mat_head,2);
mat_names = cell(n_fluences*n_opts,n_mruns);
converged_values_all = cell(n_fluences*n_opts,1); % all 20 (n_mruns) repeated runs for all cost functions and all injection levels
v_c = cell(n_opts,10); % 1st column: f; 2nd: best; 3rd: max; 4th: min; 5th: mean; 6th: std)
n_para = size(mat_para,2);

%% put all mat file names in mat_names

for i_opt = 1:1:n_opts
    m_values_cell = cell(1,n_fluences);
    for i_fluence = 1:1:n_fluences
        m_values_cell{i_fluence} = strjoin(g_set(i_fluence).m_values,'_');
        for i_run_c = 1:1:n_mruns
                mat_names{(i_opt-1)*n_fluences+i_fluence,i_run_c} = [g_set(1).project_dir,mat_head{i_opt},g_set(1).sample_name,m_methods_all,'_',m_values_cell{i_fluence},'_',g_set(1).m_units,m_date{i_opt},'_01_',mat_para,'__',num2str(i_run_c,'%02d'),'.mat'];
        end
    end
end

%% read all data

for i_fluence = 1:1:n_fluences % i_fluence is also the index of the plot
    for i_opt = 1:1:n_opts
        i_row_in_mat = (i_opt-1)*n_fluences+i_fluence;
        converged_values_all{i_row_in_mat,1} = nan(n_mruns,n_para); % row: runs; column: steps; value: cost value
        for i_run_c = 1:1:n_mruns
            if isfile(mat_names{i_row_in_mat,i_run_c})
                    clearvars all_data comsol_gbo comsol_sa prm_best best_prm;
                    load(mat_names{i_row_in_mat,i_run_c});
                if ~isempty(strfind(mat_head{i_opt},'ainu'))
                    all_data = sortrows(all_data,"cost_bayes","ascend");
                    for i_prm = 1:1:n_para
                        converged_values_all{i_row_in_mat,1}(i_run_c,i_prm) = all_data.(all_data.Properties.VariableNames{end-n_para+i_prm})(1)/target_value(i_prm); %parametric_values_all{i_row_in_mat,1}(i_run_c,1);
                     end
                else 
                    best_prm = 10.^prm_best;
                    converged_values_all{i_row_in_mat,1}(i_run_c,:) = best_prm./target_value;
                end
            else
                disp([mat_names{i_row_in_mat,i_run_c},' does not exist!']);
                continue
            end
        end
        conv1 = converged_values_all{i_row_in_mat,1}-1;
        v_c{i_opt,1} = i_opt;
        v_c{i_opt,2} = nan(1,n_para);
        for ib_prm = 1:1:n_para
            ib_prm_run = abs(conv1(:,ib_prm))==min(abs(conv1(:,ib_prm)));
            v_c{i_opt,2}(1,ib_prm) = converged_values_all{i_row_in_mat,1}(ib_prm_run,ib_prm);
        end
        v_c{i_opt,3} = max(converged_values_all{i_row_in_mat,1});
        v_c{i_opt,4} = min(converged_values_all{i_row_in_mat,1});
        v_c{i_opt,5} = mean(converged_values_all{i_row_in_mat,1});
        v_c{i_opt,6} = std(converged_values_all{i_row_in_mat,1});

    end

end

%% box plot

data_index_factor = 0;
data_box = [];
data_index = [];
for ip_prm=1:1:n_para
    for ip_opt = 1:1:n_opts
        data_index_factor = data_index_factor + 1;
        data_box = [data_box; converged_values_all{ip_opt,1}(:,ip_prm)];
        data_index = [data_index; data_index_factor*ones(size(converged_values_all{ip_opt,1}(:,ip_prm)))];
    end
end
opt_fig = figure('visible',false);
opt_fig.Units = "centimeters";
opt_fig.Position = [3,3,4.5,3];
boxplot(data_box,data_index,'Colors','rgb','Symbol','+m','OutlierSize',2,'Whisker',10) %
h=findobj('LineStyle','--');
set(h, 'LineStyle',':');
ylabel('Extracted');
set(gca, 'YScale', 'log')
ax = gca;
ax.XRuler.TickLabelInterpreter = 'tex';
xticklabels({'\it N_t \rm (AiNU)','\it N_t \rm (GB)','\it N_t \rm (SA)' ...
            ,'\it N_d \rm (AiNU)','\it N_d \rm (GB)','\it N_d \rm (SA)' ...
            ,'\it \mu_p \rm (AiNU)','\it \mu_p \rm (GB)','\it \mu_p \rm (SA)' ...
            ,'\it B_r \rm (AiNU)','\it B_r \rm (GB)','\it B_r \rm (SA)' ...
            });
xtickangle(90);
ylim([1e-6, 1e6]);
yticks([1e-6, 1e-3, 1e0, 1e3, 1e6]);
set(gca,'FontSize',6);
exportgraphics(opt_fig,['../tif/',img_name],'ContentType','vector');