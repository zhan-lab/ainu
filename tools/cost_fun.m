function [cost_val,v_intp,dvdt_intp,rms_mat] = cost_fun(t_raw,v_raw,t_target_intp,v_target_intp,dvdt_target,v_target_intp_lg,dvdt_target_lg,ratio_rms,norm_flag,intp_scale,intp_method)
    if strcmp(intp_scale,'linear')
        v_intp = interp1(t_raw,v_raw,t_target_intp,intp_method,'extrap');
    elseif strcmp(intp_scale,'log')
        v_intp_lg10 = interp1(log10(t_raw),log10(v_raw),log10(t_target_intp),intp_method,'extrap');
        v_intp = 10.^v_intp_lg10;
    else
        disp('Error: intp_scale unrecognizable!');
    end
    if norm_flag
        v_intp = v_intp/v_intp(1); % only enable this for normalized data, e.g., trpl
    end
    n_intp_pts = size(v_intp,2);
    n_sqrt = sqrt(n_intp_pts);
    v_intp_lg = log10(v_intp);
    delta_v = v_intp - v_target_intp;
    delta_v_norm = norm(delta_v);
    
    rms_v = delta_v_norm/abs(v_target_intp(end)); 
    rms_v_lg = norm(v_intp_lg - v_target_intp_lg)/abs(v_target_intp_lg(end));
    rms_lg_v = norm(log10(abs(delta_v)))/n_sqrt;
    rms_v4 = norm(delta_v./v_target_intp);
    rmse = delta_v_norm/n_sqrt; 
    
    i_max_intp = find(v_intp==max(v_intp), 1);
    i_max_target = find(v_target_intp==max(v_target_intp), 1);
    rms_t_peak = sqrt((i_max_intp-i_max_target)^2)/i_max_target; 
    
    dvdt_intp = d_fun(v_intp,1); 
    dvdt_intp_lg = d_fun(v_intp_lg,1); 
    rms_dvdt = norm(dvdt_intp - dvdt_target)/abs(max(dvdt_target));
    rms_dvdt_lg = norm(dvdt_intp_lg - dvdt_target_lg)/abs(max(dvdt_target_lg)); 
    
    rms_vector = [rms_v, rms_t_peak, rms_dvdt];
    rmse_vector = [rmse, rms_t_peak, rms_dvdt];
    rms_vector4 = [rms_v4, rms_t_peak, rms_dvdt];
    rms_vector_lg = [rms_v_lg, rms_t_peak, rms_dvdt_lg];
    rms_lg_vector = [rms_lg_v, rms_t_peak, rms_dvdt_lg];
    
    rms_mat = [rms_v, rms_v4, rms_dvdt, rms_vector_lg, rms_lg_v]; 
    
    if ratio_rms(end)==0 
        cost_val = ratio_rms(1:3)*rms_vector';
    elseif ratio_rms(end)==1 
        cost_val = log10( ratio_rms(1:3)*rms_vector' );
    elseif ratio_rms(end)==2 
        cost_val = ratio_rms(1:3)*rms_vector_lg';
    elseif ratio_rms(end)==3 
        cost_val = log10( ratio_rms(1:3)*rms_vector_lg' );
    elseif ratio_rms(end)==4 
        cost_val = ratio_rms(1:3)*rms_vector4';
    elseif ratio_rms(end)==5 
        cost_val = log10( ratio_rms(1:3)*rms_vector4' );
    elseif ratio_rms(end)==6 
        cost_val = ratio_rms(1:3)*rmse_vector';
    elseif ratio_rms(end)==7 
        cost_val = log10( ratio_rms(1:3)*rmse_vector' );
    end
end