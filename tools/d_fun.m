% numerical derivative: df/dx
function [dfdx] = d_fun(f_x, x_in)
    if isrow(f_x)
        df = [f_x(2:end),f_x(end)] - [f_x(1),f_x(1:end-1)];
    else
        df = [f_x(2:end);f_x(end)] - [f_x(1);f_x(1:end-1)];
    end
    if max(size(x_in))<=1
        dx = x_in;
    else
        if isrow(x_in)
            dx = [x_in(2:end),x_in(end)] - [x_in(1),x_in(1:end-1)];
        else
            dx = [x_in(2:end);x_in(end)] - [x_in(1);x_in(1:end-1)];
        end
    end
    dfdx = df ./ dx;
end