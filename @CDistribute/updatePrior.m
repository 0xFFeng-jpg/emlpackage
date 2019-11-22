function obj = updatePrior(obj)
[v_s,component] = size(obj.pi_x);
ct = obj.CompType;
if(obj.CompType > component)
    obj.pi_x = obj.opt_pi_x;
    obj.b_x  = obj.opt_b_x;
    obj.ba_x = obj.opt_ba_x;
else
    obj.pi_x(1,ct) = obj.opt_pi_x(1,ct);
    obj.b_x(1,ct) = obj.opt_b_x(1,ct);
    obj.ba_x(1,ct) = obj.opt_ba_x(1,ct);
end