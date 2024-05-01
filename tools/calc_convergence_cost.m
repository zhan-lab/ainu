% input: file name, t_uni_file
% output: recalculated convergence cost (c6) vs steps
% visible_flag set to true for debug
%% set the keywords for the file path and name 

function [steps, cost_r] = calc_convergence_cost(mat_file,t_uni,visible_flag,f_type)
%% set the universal time axis and the cost matrix for cost recalculation/comparison

% t_uni_file = mat_files{1};
% load(t_uni_file);
% t_uni = exp_data{1}.t_intp_trim{1};

c_matrix = struct('ratio_rms',[1, 0, 0, f_type]...
                 ,'norm_flag',true...               % flag to normalize the caculated data using the first data point, set true for trpl
                 ,'intp_scale','linear'...             
                 ,'intp_method','linear'...
                 );

% t_bounds = [1, 1.995e3];
% n_pts_intp = 1001;
% t_uni = logspace(log10(t_bounds(1)),log10(t_bounds(2)),n_pts_intp)*1e-9;

% mat_size = size(mat_file,2);
c_results = cell(1,5);
%% read all relevant files

if visible_flag
    fig_conv = figure;
    % fig_conv = figure('visible',visible_flag);
end

% for 1 = 1:1:mat_size
    % if ~isfile(mat_file)
    %     continue
    % end
    load([mat_file]);

    c_results{1,1} = all_data;
    c_results{1,2} = exp_data;
%% recalculate the cost

    n_input = size(c_results{1,2}{1}.in_target,2); % the number of fluences in one fitting file
    c_results{1,3} = struct('y_target_intp',cell(n_input,1)...
        ,'y_target_intp_lg',cell(n_input,1)...
        ,'dydx_target_intp',cell(n_input,1)...
        ,'dydx_target_intp_lg',cell(n_input,1)...
        );
    run_size = size(c_results{1,1}.step_bayes,1);
    c_results{1,4} = nan(n_input,run_size);
    c_results{1,5} = cell(1,2);

    for i_input = 1:1:n_input

        % 1, recalculate the target corresponding to the new time - t_uni
        if strcmp(c_matrix.intp_scale,'linear')
            c_results{1,3}(i_input).y_target_intp = interp1(c_results{1,2}{1}.x_target{1,i_input},c_results{1,2}{1}.y_target{1,i_input},t_uni,c_matrix.intp_method,'extrap');
        elseif strcmp(c_matrix.intp_scale,'log')
            y_target_intp_lg = interp1(log10(c_results{1,2}{1}.x_target{1,i_input}),log10(c_results{1,2}{1}.y_target{1,i_input}),log10(t_uni),c_matrix.intp_method,'extrap');
            c_results{1,3}(i_input).y_target_intp = 10.^y_target_intp_lg;
        else
            disp('Error: intp_scale unrecognizable!');
        end
        if c_matrix.norm_flag
            c_results{1,3}(i_input).y_target_intp = c_results{1,3}(i_input).y_target_intp./c_results{1,3}(i_input).y_target_intp(1);
        end

        c_results{1,3}(i_input).dydx_target_intp = d_fun(c_results{1,3}(i_input).y_target_intp,1);
        c_results{1,3}(i_input).y_target_intp_lg = log10(c_results{1,3}(i_input).y_target_intp);
        c_results{1,3}(i_input).dydx_target_intp_lg = d_fun(c_results{1,3}(i_input).y_target_intp_lg,1);

        % 2, recalculate the cost for each run
        for i_run = 1:1:run_size
            if ~isnan(c_results{1,1}.cost_bayes(i_run))
                % [c_results{1,4}(i_input,i_run),v_intp,dvdt_intp,rms_mat]
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

    c_results{1,5}{1,1} = mean(c_results{1,4},1); % average cost for all fluences
    c_results{1,5}{1,2} = cummin(c_results{1,5}{1,1}); % cumulative minimum for convergence plot
%% plot the convergence

    lgd_name = erase(mat_file,["data","baseline","/","comsol_conv_","sim_",".mat"]);
    lgd_name  = strrep(lgd_name,'_',' ');
    steps = c_results{1,1}.step_bayes;
    cost_r = c_results{1,5}{1,2};

    if visible_flag
        semilogy(steps,cost_r,'DisplayName',lgd_name);
        lgd = legend;
        hold on;

        % end

        hold off;
        xlabel('Step');
        ylabel('Cost - normalized norm');
        lgd.Location = 'southwest';
    end

% ylim([0 1000]);
% ylim([0 10]);
% disp(num2str(1));
% saveas(fig_conv,save_path)
end