function obj = calcKLDistance(obj)
%计算混合分布的中参数的KL距离。感慨:发现ensemble learning这东西真灵活。
b_x = obj.b_x;
ba_x = obj.ba_x;
pi_x = obj.pi_x;

opt_b_x = obj.opt_b_x;
opt_ba_x = obj.opt_ba_x;
opt_pi_x = obj.opt_pi_x;

b_x_prior = 1e-3;
a_x_prior = 1e-3;
pi_x_prior = 1;

I1 = sum(obj.Hx,1);

I2 = sum(gammaln(b_x_prior) - b_x_prior*log(a_x_prior) -...
    gammaln(b_x) + b_x.*log(b_x./ba_x) + (b_x - opt_b_x).*(log(ba_x)-0.5./b_x) + ...
    (opt_b_x./opt_ba_x - b_x./ba_x).*ba_x,2);

I3 = sum(gammaln(pi_x_prior) - gammaln(pi_x) + (pi_x - opt_pi_x).*(log(pi_x)-0.5./pi_x),2) - ...
    gammaln(size(pi_x,2)*pi_x_prior) + gammaln(sum(pi_x,2)) + sum(pi_x - opt_pi_x,2).*(-log(sum(pi_x,2))+0.5/sum(pi_x,2));

I4 = sum(sum(obj.c_log_lamda.*exp(obj.c_log_lamda),1),2);

obj.KLDistance = I1 + I2 + I3 + I4;