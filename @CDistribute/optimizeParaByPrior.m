function obj = optimizeParaByPrior(obj)
%根据先验知识，优化分布参数。
[d_size,component] = size(obj.c_log_lamda);
t_log_lamda = obj.c_log_lamda;
b_x_prior = 1e-3;
pi_x_prior = 1;
a_x_prior = 1e-3;
distType = obj.DistriType;

MANUAL = 1;

switch distType
    case 'Gaussian'
        obj.opt_b_x = b_x_prior + reshape(sum(exp(t_log_lamda),1),1,component)/2;
        obj.opt_pi_x = pi_x_prior + reshape(sum(exp(t_log_lamda),1),1,component);
        obj.opt_ba_x = zeros(1,component);
        for alpha = 1:component
            obj.opt_ba_x(1,alpha) = obj.opt_b_x(1,alpha) ./ (a_x_prior + sum(exp(t_log_lamda(:,alpha)).*obj.mx2,1)/2);
        end
    case {'Exponential','Laplacian'}
        obj.opt_b_x = b_x_prior + reshape(sum(exp(t_log_lamda),1),1,component);
        obj.opt_pi_x = pi_x_prior + reshape(sum(exp(t_log_lamda),1),1,component);
        obj.opt_ba_x = zeros(1,component);
        for alpha = 1:component
            obj.opt_ba_x(1,alpha) = obj.opt_b_x(1,alpha) ./ (a_x_prior + sum(exp(t_log_lamda(:,alpha)).*obj.mx,1));
        end
    case 'Discrete'
        disp('To be settled');
    otherwise
        error('Unknown Distribution');
end


if(obj.IsManual)
    switch obj.ObjType
        case 'Image'
            obj.opt_b_x = ones(1,component)*1e-3;
            obj.opt_pi_x = [365.4685 444.4363 3.8208 186.2744];%[0.3655 0.4444 0.0038 0.1863]*1e3;
            obj.opt_ba_x = [0.6655 21.1467 3.1931e-4 0.0232];
        case 'Blur'
            obj.opt_b_x = [787.8988 201.7349 236.1948 143.1756];
            obj.opt_ba_x = [5.1143e3 5.0064e3 173.8885 50.6538];
        otherwise
    end
end
 