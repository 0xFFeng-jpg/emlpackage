function obj = calcOptCx(obj)

%计算与卷积的接口；
obj.opt_cx1 = zeros(size(obj.x1));
obj.opt_cx2 = zeros(size(obj.x2));

components = size(obj.pi_x,2);

distType = obj.DistriType;

switch distType
    case {'Gaussian'}
        for alpha = 1:components
            obj.opt_cx2 = obj.opt_cx2 + obj.ba_x(1,alpha)*exp(obj.c_log_lamda(:,alpha));
        end
    case {'Laplacian','Exponential'}
        for alpha = 1:components
            obj.opt_cx1 = obj.opt_cx1  - obj.ba_x(1,alpha)*exp(obj.c_log_lamda(:,alpha));
        end
    otherwise
        error('No such distribute');
end
