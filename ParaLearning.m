function [mx_est,me_est] = ParaLearning(mx_init,me_init,mx_org,IsStructure)

%FFT模式；

[K_I,K_J] = size(me_init);
[I_I,I_J] = size(mx_init);
D_I = 2*I_I;
D_J = 2*I_J;

Dp = zeros(D_I,D_J);

half_k = floor(K_I/2);
half_l = floor(K_J/2);

% Dp(half_k+1:I_I+half_k,half_l+1:I_J/2+half_l) = 1;
% Dp(half_k+1:I_I+half_k,K_J+1+I_J/2:I_J+K_J) = 1;

Dp(K_I:I_I,K_J:I_J/2) = 1;
Dp(K_I:I_I,K_J+I_J/2:I_J) = 1;

D = padarray(mx_org,size(mx_org),0,'post');
D = D./sqrt(mean(mean(D.^2)));

mx = mx_init;
mx = mx_init./sqrt(mean(mean(mx_init.^2)));
me = me_init;

ensemble_x1 = 0;%min(mx(:))*1e4;
ensemble_x2 = 1e4;

im_x1  = 1e4*mx(:);% - 1e4*min(mx(:));
im_x2 = 1e4*ones(size(im_x1));

% init_k = 1e-4*ones(size(me));
% 
% init_k(4,4) = 1;
init_k = me;

kernel_x1 = 1e4*init_k(:);
kernel_x2 = 1e4*ones(size(kernel_x1));

%外界口；

EnDistType = 'Gaussian';
ImDistType = 'Gaussian';
KerDistType = 'Laplacian';

IMCOMPONENTS = 4;
KERCOMPONENTS = 4;

Ensemble = CDistribute(ensemble_x1,ensemble_x2,1,1,1,EnDistType,0,'Ensemble',1,1,0);
EnGrad = struct('x1',0,'x2',0,'pi_x',0,'b_x',0,'ba_x',0);

Im = CDistribute(im_x1,im_x2,1*ones(1,IMCOMPONENTS),1*ones(1,IMCOMPONENTS),1*ones(1,IMCOMPONENTS),ImDistType,0,'Image',I_I,I_J,IsStructure); %初始化
ImGrad = struct('x1',zeros(size(im_x1)),...
    'x2',zeros(size(im_x2)),...
    'pi_x',zeros(1,IMCOMPONENTS),...
    'b_x',zeros(1,IMCOMPONENTS),...
    'ba_x',zeros(1,IMCOMPONENTS));

Ker = CDistribute(kernel_x1,kernel_x2,1*ones(1,KERCOMPONENTS),1*ones(1,KERCOMPONENTS),1*ones(1,KERCOMPONENTS),KerDistType,0,'Blur',K_I,K_J,0); %初始化
KerGrad = struct('x1',zeros(size(kernel_x1)),...
    'x2',zeros(size(kernel_x2)),...
    'pi_x',zeros(1,KERCOMPONENTS),...
    'b_x',zeros(1,KERCOMPONENTS),...
    'ba_x',zeros(1,KERCOMPONENTS));

distriList = [Ensemble Im Ker];
gradList = [EnGrad ImGrad KerGrad];

Deconv = CDeconv(distriList,D,Dp,ones(size(me)));

state = 1;
ObjEnvidence = CEnvidence(distriList,gradList,Deconv,state,D,ones(size(me)));

alpha = 1;
beta = 0.9;
last_time_D_val = NaN;
converge_criteria = 5e-4;
state = 1;
Niter = 5000;
last_change_iter = 0;

D_log = NaN * ones(1,Niter);
gamma_log = NaN * ones(1,Niter);

inDistriList = distriList;
directionList = gradList;
inDeconv = Deconv;

for iter = 1:Niter
    if(iter == last_change_iter + 1)
        if(state < 3)
            step_len = 0.0;
            [oDistriList,oGradList,oDeconv,oEnvidence] = ObjEnvidence.calcParaEnvidence(inDistriList,directionList,inDeconv,state,step_len,[5 5 5]);
            directionList = oGradList;
            inDistriList = oDistriList;
            inDeconv = oDeconv;
       
            for n_iter = 1:size(oDistriList,2)
                DistriObj = oDistriList(1,n_iter);
                
                DistriObj = DistriObj.updatePrior();

                if(DistriObj.checkPrior())
                    DistriObj = DistriObj.reinitialPrior();
                end
                inDistriList(1,n_iter) = DistriObj;
            end
        end
        step_len = 0.0;
        
        clear oDistriList oDeconv;
        [oDistriList,oGradList,oDeconv,oEnvidence] = ObjEnvidence.calcParaEnvidence(inDistriList,directionList,inDeconv,state,step_len,[5 5 5]);
      
        inDistriList = oDistriList;
        inGradList = oGradList;
        inDeconv = oDeconv;
        
       
        e_D_val = oEnvidence;
        for n_iter = 1:size(oGradList,2)
        directionList(1,n_iter).x1 = alpha*oGradList(1,n_iter).x1;
        directionList(1,n_iter).x2 = alpha*oGradList(1,n_iter).x2;
        directionList(1,n_iter).b_x = alpha*oGradList(1,n_iter).b_x;
        directionList(1,n_iter).ba_x = alpha*oGradList(1,n_iter).ba_x;
        directionList(1,n_iter).pi_x = alpha*oGradList(1,n_iter).pi_x;
        end
        step_len = 1.0;
    end
    [oDistriList,oGradList,oDeconv,oEnvidence] = ObjEnvidence.calcParaEnvidence(inDistriList,directionList,inDeconv,state,step_len,[5 5 5]);
    t_D_val = oEnvidence;
    if (t_D_val > e_D_val)
        for n_iter = 1:size(oGradList,2)
            directionList(1,n_iter).x1 = alpha * inGradList(1,n_iter).x1;
            directionList(1,n_iter).x2 = alpha * inGradList(1,n_iter).x2;
            directionList(1,n_iter).b_x = alpha * inGradList(1,n_iter).b_x;
            directionList(1,n_iter).ba_x = alpha * inGradList(1,n_iter).ba_x;
            directionList(1,n_iter).pi_x = alpha * inGradList(1,n_iter).pi_x;
        end
        step_len = 2.0 * step_len;
        while(t_D_val > e_D_val + converge_criteria/1e4)
            step_len = 0.5 * step_len;
            [oDistriList,oGradList,oDeconv,oEnvidence] = ObjEnvidence.calcParaEnvidence(inDistriList,directionList,inDeconv,state,step_len,[5 5 5]);
            t_D_val = oEnvidence;
        end
    end
    
    for n_iter = 1:size(oGradList,2)
        directionList(1,n_iter).x1 = alpha * oGradList(1,n_iter).x1 + beta * directionList(1,n_iter).x1;
        directionList(1,n_iter).x2 = alpha * oGradList(1,n_iter).x2 + beta * directionList(1,n_iter).x2;
        directionList(1,n_iter).b_x = alpha * oGradList(1,n_iter).b_x + beta * directionList(1,n_iter).b_x;
        directionList(1,n_iter).ba_x = alpha * oGradList(1,n_iter).ba_x + beta * directionList(1,n_iter).ba_x;
        directionList(1,n_iter).pi_x = alpha * oGradList(1,n_iter).pi_x + beta * directionList(1,n_iter).pi_x;
    end
    step_len = min(1,1.1*step_len);
    
    dD_val = t_D_val - last_time_D_val;
    inDistriList = oDistriList;
    inGradList = oGradList;
    inDeconv = oDeconv;
    e_D_val = t_D_val;
    
    D_log(1,iter) = e_D_val;
    gamma_log(1,iter) = inDeconv.ba_sigma;
    
    last_time_D_val = e_D_val;
    
    fprintf(1,['Iteration %4d State %1d Evidence %11.6f bits per pixel,' ...
        'Delta = %11.6e, Noise = %11.6e'],iter,state,-D_log(1,iter),-dD_val,gamma_log(1,iter)^-0.5);
    fprintf(1,'\n');
    
    
    %判定收敛条件；
    converged = 0;
    if(iter > 3)
        last_dD_log = D_log(1,iter-2:iter) - D_log(1,iter-3:iter-1);
        if(all(last_dD_log < converge_criteria/1e4 & last_dD_log > -converge_criteria))
            converged = 1;
        end
    end
    
    if(converged)
        if(state == 3)
            break;
        elseif(state == 2 && inDeconv.opt_ba_sigma < 1.1*inDeconv.ba_sigma)
            state = 3;
        else
            state = 3 - state;
        end
        last_change_iter = iter;
        if(state == 1)
            inDeconv.ba_sigma = inDeconv.opt_ba_sigma;
        end
        
    end
end

D_log = D_log(1:iter);
gamma_log = gamma_log(1,1:iter);

[m,n] = size(mx_init);
[k,l] = size(me_init);

out_mmu = inDistriList(1,1).mx
out_mx = reshape(inDistriList(1,2).mx,m,n);
out_me = reshape(inDistriList(1,3).mx,k,l);

%mD = real(ifft2(fft2(out_mx).*psf2otf(out_me,[m,n]))) + out_mmu;

mx_est = out_mx;
me_est = out_me;

% figure;
% subplot(2,3,1);
% imagesc(out_mx)
% axis square
% colorbar
% title('Recovered Source');
% subplot(2,3,2)
% imagesc(out_me)
% axis square
% colorbar
% title('Recovered Filter');
% subplot(2,3,3)
% imagesc(mD)
% axis square
% colorbar
% title('Reconstructed Data');
% subplot(2,3,4);
% imagesc(D)
% axis square
% colorbar
% title('Data');
% subplot(2,3,5)
% imagesc(me_init)
% axis square
% colorbar
% title('Ori Me');
% subplot(2,3,6)
% imagesc(mx_init)
% axis square
% colorbar
% title ('Ori Mx');