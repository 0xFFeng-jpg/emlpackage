clear all;
%%
%合成数据
x = ((-3:3)'*ones(1,7));
y = (ones(7,1)*(-3:3));
me = exp(-(x.^2 + x.*y + y.^2)/2/1^2);
mx = zeros(64);
mx(ceil(rand(1,20)*64*64)) = -log(rand(1,20));

D = real(ifft2(fft2(mx).*psf2otf(me,size(mx))));
D = D/sqrt(mean(mean(D.^2)));

%%
%初始化
im_x1  = 1e4*D(:);
im_x2 = 1e4*ones(size(im_x1));

init_k = 1e-4*ones(size(me));

init_k(4,4) = 1;

kernel_x1 = 1e4*init_k(:);
kernel_x2 = 1e4*ones(size(kernel_x1));

Im = CDistribute(im_x1,im_x2,1*ones(1,2),1*ones(1,2),1*ones(1,2)); %初始化
Ker = CDistribute(kernel_x1,kernel_x2,1*ones(1,2),1*ones(1,2),1*ones(1,2)); %初始化

ImDistType = 'Gaussian';
KerDistType = 'Laplacian';

alpha = 1;
last_time_D_val = NaN;
converge_criteria = 1e-4;
state = 1;
Niter = 5000;
last_change_iter = 0;

D_log = NaN * ones(1,Niter);
gamma_log = NaN * ones(1,Niter);


%%
%这里是一个大循环：
%第一层循环主要是每个状态的的变化；每个状态变化主要包括
%1：初始模型中的参数的判定和学习，如果参数过于接近则要进行分离淘汰；
for iter = 1:Niter
    if(iter == last_change_iter + 1) %状态转换
        if(state < 3) %state 1: 更新数据 2：更新整体（包括数据分布和噪声分布）3：在新的数据和参数基础上，全面更新；
            
            Im = Im.calcExpection(ImDistType);
            Im = Im.calcLamda(ImDistType);
            Im = Im.optimizeParaByPrior(ImDistType);
            
            Im.pi_x = Im.opt_pi_x;
            Im.b_x = Im.opt_b_x;
            Im.ba_x = Im.opt_ba_x;
            
            if(Im.checkPrior())
                Im = Im.reinitialPrior();
            end
            
            
            Ker = Ker.calcExpection(KerDistType);
            Ker = Ker.calcLamda(KerDistType);
            Ker = Ker.optimizeParaByPrior(KerDistType);
            
            Ker.pi_x = Ker.opt_pi_x;
            Ker.b_x = Ker.opt_b_x;
            Ker.ba_x = Ker.opt_ba_x;
            
            if(Ker.checkPrior())
                Ker = Ker.reinitialPrior();
            end
            
        end
        Im = Im.calcLamda(ImDistType);
        Im = Im.optimizeParaByPrior(ImDistType);
        Im = Im.calcOptCx(ImDistType);
        Im = Im.calcKLDistance();
        
        if(state >= 2)
            grad_im_pi = Im.opt_pi_x - Im.pi_x;
            grad_im_b = Im.opt_b_x - Im.b_x;
            grad_im_ba = Im.opt_ba_x - Im.ba_x;
        else
            grad_im_pi = 0;
            grad_im_b = 0;
            grad_im_ba = 0;
        end
        
        Ker = Ker.calcLamda(KerDistType);
        Ker = Ker.optimizeParaByPrior(KerDistType);
        Ker = Ker.calcOptCx(KerDistType);
        Ker = Ker.calcKLDistance();
        
        if(state >= 2)
            grad_ker_pi = Ker.opt_pi_x - Ker.pi_x;
            grad_ker_b = Ker.opt_b_x - Ker.b_x;
            grad_ker_ba = Ker.opt_ba_x -Ker.ba_x;
        else
            grad_ker_pi = 0;
            grad_ker_b = 0;
            grad_ker_ba = 0;
        end
        
        oDeconv = CDeconv(Im,Ker,D,ones(size(me)));
        oDeconv = oDeconv.blindDeconv;
        oDeconv = oDeconv.calcSigmaDx((state == 3));
        oDeconv = oDeconv.updateExpection;
        
        grad_im_x1 = oDeconv.im_x1 - Im.x1;
        grad_im_x2 = oDeconv.im_x2 - Im.x2;
        
        grad_kernel_x1 = oDeconv.kernel_x1 - Ker.x1;
        grad_kernel_x2 = oDeconv.kernel_x2 - Ker.x2;
        
        direction_im_x1 = alpha * grad_im_x1;
        direction_im_x2 = alpha * grad_im_x2;
        direction_im_b = alpha * grad_im_b;
        direction_im_ba = alpha * grad_im_ba;
        direction_im_pi = alpha * grad_im_pi;
        
        direction_kernel_x1 = alpha * grad_kernel_x1;
        direction_kernel_x2 = alpha * grad_kernel_x2;
        
        direction_ker_b = alpha * grad_ker_b;
        direction_ker_ba = alpha * grad_ker_ba;
        direction_ker_pi = alpha * grad_ker_pi;
        
    end
    
    %在同一状态下，根据梯度形式进行最优化。
    step_len = 1.0;
    
    D_val = (Im.KLDistance + Ker.KLDistance + oDeconv.D_x)/oDeconv.data_points;
    
    new_Im = Im;
    new_Ker = Ker;
 
    %梯度更新;
    new_Im.x1 = Im.x1 + step_len * direction_im_x1;
    new_Im.x2 = abs(Im.x2 + step_len * direction_im_x2);
    
    new_Im.b_x = abs(Im.b_x + step_len * direction_im_b);
    new_Im.ba_x = abs(Im.ba_x + step_len * direction_im_ba);
    new_Im.pi_x = abs(Im.pi_x + step_len * direction_im_pi);
    
    new_Im = new_Im.calcExpection(ImDistType);
    new_Im = new_Im.calcLamda(ImDistType);
    new_Im = new_Im.optimizeParaByPrior(ImDistType);
    new_Im = new_Im.calcKLDistance();
    new_Im = new_Im.calcOptCx(ImDistType);
    
   
    
    new_Ker.x1 = Ker.x1 + step_len * direction_kernel_x1;
    new_Ker.x2 = abs(Ker.x2 + step_len * direction_kernel_x2);
    
    new_Ker.b_x = abs(Ker.b_x + step_len * direction_ker_b);
    new_Ker.ba_x = abs(Ker.ba_x + step_len * direction_ker_ba);
    new_Ker.pi_x = abs(Ker.pi_x + step_len * direction_ker_pi);
    
    new_Ker = new_Ker.calcExpection(KerDistType);
    new_Ker = new_Ker.calcLamda(KerDistType);
    new_Ker = new_Ker.optimizeParaByPrior(KerDistType);
    new_Ker = new_Ker.calcKLDistance();
    new_Ker = new_Ker.calcOptCx(KerDistType);
    
    new_Deconv = oDeconv;
    
    new_Deconv.ImD = new_Im;
    new_Deconv.KernelD = new_Ker;
    
    new_Deconv = new_Deconv.blindDeconv;
    new_Deconv = new_Deconv.calcSigmaDx(state == 3);
    new_Deconv = new_Deconv.updateExpection;
    
    %计算新梯度；
    if(state >= 2)
        new_grad_im_pi = new_Im.opt_pi_x - new_Im.pi_x;
        new_grad_im_ba = new_Im.opt_ba_x - new_Im.ba_x;
        new_grad_im_b = new_Im.opt_b_x - new_Im.b_x;
        
        new_grad_ker_pi = new_Ker.opt_pi_x - new_Ker.pi_x;
        new_grad_ker_ba = new_Ker.opt_ba_x - new_Ker.ba_x;
        new_grad_ker_b  = new_Ker.opt_b_x - new_Ker.b_x;
    else
        new_grad_im_pi = 0;
        new_grad_im_ba = 0;
        new_grad_im_b = 0;
        
        new_grad_ker_pi = 0;
        new_grad_ker_ba = 0;
        new_grad_ker_b = 0;
    end
    
    new_grad_im_x1 = new_Deconv.im_x1 - new_Im.x1;
    new_grad_im_x2 = new_Deconv.im_x2 - new_Im.x2;
    
    new_grad_kernel_x1 = new_Deconv.kernel_x1 - Ker.x1;
    new_grad_kernel_x2 = new_Deconv.kernel_x2 - Ker.x2;
    
    new_D_val = (new_Im.KLDistance + new_Ker.KLDistance + new_Deconv.D_x)./new_Deconv.data_points;
    
    if(new_D_val > D_val)
        direction_im_x1 = alpha*grad_im_x1;
        direction_im_x2 = alpha*grad_im_x2;
        
        direction_im_pi = alpha*grad_im_pi;
        direction_im_ba = alpha * grad_im_ba;
        direction_im_b  = alpha * grad_im_b;
        
        direction_kernel_x1 = alpha*grad_kernel_x1;
        direction_kernel_x2 = alpha * grad_kernel_x2;
        
        direction_ker_pi = alpha * grad_ker_pi;
        direction_ker_ba = alpha * grad_ker_ba;
        direction_ker_b  = alpha * grad_ker_b;
        
        step_len = 2 * step_len;
        while(new_D_val > D_val + converge_criteria/1e4) %方向不变，收敛速度变化
            step_len = 0.5 * step_len;
            new_Im = Im;
            new_Ker = Ker;
            new_Im.x1 = Im.x1 + step_len * direction_im_x1;
            new_Im.x2 = abs(Im.x2 + step_len * direction_im_x2);
            
            new_Im.b_x = abs(Im.b_x + step_len * direction_im_b);
            new_Im.ba_x = abs(Im.ba_x + step_len * direction_im_ba);
            new_Im.pi_x = abs(Im.pi_x + step_len * direction_im_pi);
            
            new_Im = new_Im.calcExpection(ImDistType);
            new_Im = new_Im.calcLamda(ImDistType);
            new_Im = new_Im.optimizeParaByPrior(ImDistType);
            new_Im = new_Im.calcKLDistance();
            new_Im = new_Im.calcOptCx(ImDistType);
            
            
            
            new_Ker.x1 = Ker.x1 + step_len * direction_kernel_x1;
            new_Ker.x2 = abs(Ker.x2 + step_len * direction_kernel_x2);
            
            new_Ker.b_x = abs(Ker.b_x + step_len * direction_ker_b);
            new_Ker.ba_x = abs(Ker.ba_x + step_len * direction_ker_ba);
            new_Ker.pi_x = abs(Ker.pi_x + step_len * direction_ker_pi);
            
            new_Ker = new_Ker.calcExpection(KerDistType);
            new_Ker = new_Ker.calcLamda(KerDistType);
            new_Ker = new_Ker.optimizeParaByPrior(KerDistType);
            new_Ker = new_Ker.calcKLDistance();
            new_Ker = new_Ker.calcOptCx(KerDistType);
            
            new_Deconv = oDeconv;
            
            new_Deconv.ImD = new_Im;
            new_Deconv.KernelD = new_Ker;
            
            new_Deconv = new_Deconv.blindDeconv;
            new_Deconv = new_Deconv.calcSigmaDx(state == 3);
            new_Deconv = new_Deconv.updateExpection;
            
            %计算新梯度；
            if(state >= 2)
                new_grad_im_pi = new_Im.opt_pi_x - new_Im.pi_x;
                new_grad_im_ba = new_Im.opt_ba_x - new_Im.ba_x;
                new_grad_im_b = new_Im.opt_b_x - new_Im.b_x;
                
                new_grad_ker_pi = new_Ker.opt_pi_x - new_Ker.pi_x;
                new_grad_ker_ba = new_Ker.opt_ba_x - new_Ker.ba_x;
                new_grad_ker_b  = new_Ker.opt_b_x - new_Ker.b_x;
            else
                new_grad_im_pi = 0;
                new_grad_im_ba = 0;
                new_grad_im_b = 0;
                
                new_grad_ker_pi = 0;
                new_grad_ker_ba = 0;
                new_grad_ker_b = 0;
            end
            
            new_grad_im_x1 = new_Deconv.im_x1 - new_Im.x1;
            new_grad_im_x2 = new_Deconv.im_x2 - new_Im.x2;
            
            new_grad_kernel_x1 = new_Deconv.kernel_x1 - Ker.x1;
            new_grad_kernel_x2 = new_Deconv.kernel_x2 - Ker.x2;
            
            new_D_val = (new_Im.KLDistance + new_Ker.KLDistance + new_Deconv.D_x)./new_Deconv.data_points;
        end
    end
    
    %更新方向；
    beta = 0.9;
    direction_im_x1 = alpha * new_grad_im_x1 + beta * direction_im_x1;
    direction_im_x2 = alpha * new_grad_im_x2 + beta * direction_im_x2;
    
    direction_im_b = alpha * new_grad_im_b + beta * direction_im_b;
    direction_im_ba = alpha * new_grad_im_ba + beta * direction_im_ba;
    direction_im_pi = alpha * new_grad_im_pi + beta * direction_im_pi;
    
    direction_kernel_x1 = alpha * new_grad_kernel_x1 + beta * direction_kernel_x1;
    direction_kernel_x2 = alpha * new_grad_kernel_x2 + beta * direction_kernel_x2;
    
    direction_ker_b = alpha * new_grad_ker_b + beta * direction_ker_b;
    direction_ker_ba = alpha * new_grad_ker_ba + beta * direction_ker_ba;
    direction_ker_pi = alpha * new_grad_ker_pi + beta * direction_ker_pi;
    
    step_len = min(1,1.1*step_len);
    
    dD_val = new_D_val - last_time_D_val;%last 为上次循环更新后的总KL距离；
    
    %重新初始化；
    Im = new_Im;
    Ker = new_Ker;
    oDeconv = new_Deconv; %
    
    grad_im_x1 = new_grad_im_x1;
    grad_im_x2 = new_grad_im_x2;
    grad_im_pi = new_grad_im_pi;
    grad_im_b = new_grad_im_b;
    grad_im_ba = new_grad_im_ba;
    
    grad_kernel_x1 = new_grad_kernel_x1;
    grad_kernel_x2 = new_grad_kernel_x2;
    grad_ker_pi = new_grad_ker_pi;
    grad_ker_b = new_grad_ker_b;
    grad_ker_ba = new_grad_ker_ba;
    
    
    D_log(1,iter) = D_val;
    gamma_log(1,iter) = oDeconv.ba_sigma;
    
    last_time_D_val = D_val;
    
    fprintf(1,['Iteration %4d State %1d Evidence %11.6f bits per pixel,' ...
        'Delta = %11.6e, Noise = %11.6e'],iter,state,-D_log(1,iter),-dD_val,gamma_log(1,iter)^-0.5);
    fprintf(1,'\n');
   
    
    %判定收敛条件；
    converged = 0;
    if(iter > 10)
        last_dD_log = D_log(1,iter-9:iter) - D_log(1,iter-10:iter-1);
        if(all(last_dD_log < converge_criteria/1e4 & last_dD_log > -converge_criteria))
            converged = 1;
        end
    end
    
    if(converged)
        if(state == 3)
            break;
        elseif(state == 2 && oDeconv.opt_ba_sigma < 1.1*oDeconv.ba_sigma)
            state = 3;
        else
            state = 3 - state;
        end
        last_change_iter = iter;
        if(state == 1)
            oDeconv.ba_sigma = oDeconv.opt_ba_sigma;
        end
        
    end
end
    
D_log = D_log(1:iter);
gamma_log = gamma_log(1,1:iter);
        
    
    
    
    
    
    
    
    
            
        

