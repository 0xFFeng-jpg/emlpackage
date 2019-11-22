function obj = reinitialPrior(obj)

pi_x_prior = 1;
b_x_prior = 1e-3;
a_x_prior = 1e-3;

[m,n] = size(obj.x1);
component = size(obj.pi_x,2);

tmp = m*n/component;

c_pi_x = obj.pi_x;
c_b_x = obj.b_x;
c_ba_x = obj.ba_x;

mean_scale = sum(obj.b_x(1,:)./obj.ba_x(1,:))/sum(obj.b_x(1,:));
c_pi_x(1,:) = pi_x_prior + tmp;
c_b_x(1,:) = b_x_prior + tmp;

c_ba_x(1,:) = c_b_x(1,:)./(a_x_prior + 0.5*(1:component)*mean_scale*tmp);

obj.pi_x = c_pi_x;
obj.b_x = c_b_x;
obj.ba_x = c_ba_x;