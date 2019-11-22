function obj = calcSigmaDx(obj,isUpdate)
b_sigma_prior = 1e-3;
a_sigma_prior = 1e-3;

obj.b_sigma = b_sigma_prior + obj.data_points/2;
obj.opt_ba_sigma = obj.b_sigma / (a_sigma_prior + obj.error/2);

if(isUpdate)
    obj.ba_sigma = obj.opt_ba_sigma;
end

obj.D_x = gammaln(b_sigma_prior) - gammaln(obj.b_sigma) - b_sigma_prior*log(a_sigma_prior) + ...
    obj.b_sigma*log(obj.b_sigma/obj.ba_sigma) + ...
    (obj.b_sigma/obj.opt_ba_sigma - obj.b_sigma/obj.ba_sigma)*obj.ba_sigma + obj.data_points*log(2*pi)/2;

