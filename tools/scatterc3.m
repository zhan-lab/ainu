% This program is free to use for acedemic purposes only
% Users must refer to AiNU manual for disclaim and copyright notice.
% Author: Hualin Zhan, Australian National University

% color of each data point is defined by the maximum and the minimum value of c_range
% color tick in the color bar is defined by colorlabel
% makesure c_range is the same as the color label

function [xyz_table,cm,color_index,colorlabel]=scatterc3(x, y, z, plot_range, red_on_top, cost_prec)
% x, y, z are n*1 double array
% plot_range: 1*2 double array, leave empty for default values
% red_on_top: set true to put red color (minimum value) on top of the colorbar

f_color = 1000;
xyz_table = table(x, y, z);
xyz_table = rmmissing(xyz_table);
xyz_table = sortrows(xyz_table,'z','ascend');

z_raw = xyz_table.z;
if ~isempty(plot_range)
    z_plot_index = find(z_raw>plot_range(1) & z_raw<plot_range(2));
    if isempty(z_plot_index)
        disp('no cost value is selected, possibly because z range is too small. attempt to adjust the z range ...');
        z_min = min(z_raw);
        if z_min >0
            z_max = 3 * z_min;
        else
            z_max = z_min + 100*abs(z_min);
        end
        
        z_plot_index = find(z_raw>=z_min & z_raw<=z_max);
        if isempty(z_plot_index)
            disp('z range is still too small. try manually adjusting the z_high value in scathist_tif.m. start this process by loading mat_path.');
            disp(['current values: z_low = ',num2str(plot_range(1)),' and z_high = ',num2str(plot_range(2))]);
        else
            disp(['adjusted values: z_low = ',num2str(z_min),' and z_high = ',num2str(z_max)]);
            % disp(['index of selected points: ',num2str(z_plot_index')]);
        end
    else
        z_min = plot_range(1);
        z_max = plot_range(2);
    end
    % disp(num2str(z_plot_index));
    index_max = size(z_raw,1);
    if z_plot_index(1) ~= 1
        z_raw(1:z_plot_index(1)-1) = z_min * ones((z_plot_index(1)-1),1);
    end
    if z_plot_index(end) ~= index_max
        z_raw((z_plot_index(end)+1):end) = z_max * ones((index_max-z_plot_index(end)),1);
    end
else
    z_min = min(z_raw);
    z_max = max(z_raw);    
end
color_total = ceil((z_max - z_min)*f_color) + 1;
cm = colormap(parula(color_total)); % color scheme: turbo, parula, hot, spring -> winter, jet

colorlabel_raw = [z_min, 0.5*(z_max+z_min), z_max];
colorlabel_num = round(colorlabel_raw*10^cost_prec)/10^cost_prec;

z_color = ceil((z_raw - z_min)*f_color) + 1;
if red_on_top
    color_index = color_total - z_color + 1;
    colorlabel_num = flip(colorlabel_num,2);
else
    color_index = z_color;
end
colorlabel0 = arrayfun(@num2str, colorlabel_num, 'UniformOutput', 0);
colorlabel = [['\geq',colorlabel0{1}] colorlabel0(2:end)];

end