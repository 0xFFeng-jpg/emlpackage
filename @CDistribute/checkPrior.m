function re = checkPrior(obj)
%检查先验是否符合条件；
component = size(obj.pi_x,2);
if(component > 1)
    
    sorted_scales = sort(obj.ba_x(1,:));
    
    restart_priors = 0;
    
    pi_x_prior = 1;
    
    
    
    if(restart_priors | any(obj.pi_x(1,:) < pi_x_prior + 1/size(obj.pi_x,2)) | any(sorted_scales(2:end) < 1.5 * sorted_scales(1:end-1)))
        re = 1;
    else
        re = 0;
    end
    
else
    re = 0;
end