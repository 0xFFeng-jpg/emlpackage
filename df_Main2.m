%clear all;

%%
%合成数据
x = ((-3:3)'*ones(1,7));
y = (ones(7,1)*(-3:3));
me = exp(-(x.^2 + x.*y + y.^2)/2/1^2);
mx_org = zeros(64);
mx_org(ceil(rand(1,20)*64*64)) = -log(rand(1,20));

%mx(1:100:end) = 1;



%%


[K_I,K_J] = size(me);
[I_I,I_J] = size(mx_org);
D_I = 2*I_I;
D_J = 2*I_J;

Dp = zeros(D_I,D_J);

half_k = floor(K_I/2);
half_l = floor(K_J/2);


%理论上的完全对应关系；
Dp(half_k+1:half_k+I_I,half_l+1:half_l+I_J) = 1;


D = real(ifft2(fft2(mx_org,D_I,D_J).*fft2(me,D_I,D_J)));
D = D/sqrt(mean(mean(D.^2)));
mx = D(half_k+1:half_k+I_I,half_l+1:half_l+I_J);

%实际情况下的模拟；
%  D = real(ifft2(fft2(mx_org,I_I,I_J).*fft2(me,I_I,I_J)));
%  D = D/sqrt(mean(mean(D.^2)));
%  mx = D;
% 
%  D = padarray(D,size(D),0,'post');
%  Dp(K_I+1:K_I+I_I,K_J+1:K_J+I_J) = 1;



%me = me_init;

ensemble_x1 = min(mx(:))*1e4;
ensemble_x2 = 1e4;

im_x1  = 1e4*mx(:) - 1e4*min(mx(:));
im_x2 = 1e4*ones(size(im_x1));

init_k = 1e-4*ones(size(me));

init_k(4,4) = 1;
%init_k = me;

kernel_x1 = 1e4*init_k(:);
kernel_x2 = 1e4*ones(size(kernel_x1));

%外界口；

EnDistType = 'Gaussian';
ImDistType = 'Gaussian';
KerDistType = 'Laplacian';

IMCOMPONENTS = 2;
KERCOMPONENTS = 2;

Ensemble = CDistribute(ensemble_x1,ensemble_x2,1,1,1,EnDistType,0,'Ensemble');
EnGrad = struct('x1',0,'x2',0,'pi_x',0,'b_x',0,'ba_x',0);

Im = CDistribute(im_x1,im_x2,1*ones(1,IMCOMPONENTS),1*ones(1,IMCOMPONENTS),1*ones(1,IMCOMPONENTS),ImDistType,0,'Image'); %初始化
ImGrad = struct('x1',zeros(size(im_x1)),...
    'x2',zeros(size(im_x2)),...
    'pi_x',zeros(1,IMCOMPONENTS),...
    'b_x',zeros(1,IMCOMPONENTS),...
    'ba_x',zeros(1,IMCOMPONENTS));

Ker = CDistribute(kernel_x1,kernel_x2,1*ones(1,KERCOMPONENTS),1*ones(1,KERCOMPONENTS),1*ones(1,KERCOMPONENTS),KerDistType,0,'Blur'); %初始化
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
converge_criteria = 1e-4;
state = 1;
Niter = 5000;
last_change_iter = 0;

D_log = NaN * ones(1,Niter);
gamma_log = NaN * ones(1,Niter);

inDistriList = distriList;
directionList = gradList;
inDeconv = Deconv;

comList = [1 IMCOMPONENTS KERCOMPONENTS];



%%
%查看数据变化情况；
time_v = clock;
OutDirName = fullfile(['y' int2str(time_v(1)) 'mm' int2str(time_v(2)) 'd' int2str(time_v(3))]);
s = 0;%0 all update,1 component wise
SingleName = fullfile(['S_' int2str(s) 'h' int2str(time_v(4)) 'm' int2str(time_v(5)) 's' int2str(round(time_v(6)))]);
mkdir(['ParaCompare\' OutDirName]);
ParaFileName = ['ParaCompare\' OutDirName '\' SingleName '.mat'];
% fid = fopen(ParaFileName,'w+');

ParaList = [];
OptParaList = [];
IsParaCompare = 1;
%%




for iter = 1:Niter
    %compNum = (iter-1) - floor((iter-1)./comList).*comList + 1;
    compNum = [5 5 5];
    if(iter == last_change_iter + 1)
        if(state < 3)
            step_len = 0.0;
            %%
            
            for compIter = 1:IMCOMPONENTS
                compNum(1,2) = compIter;
                [oDistriList,oGradList,oDeconv,oEnvidence] = ObjEnvidence.calcEnvidence(inDistriList,directionList,inDeconv,state,step_len,compNum);
                
                if(IsParaCompare)
                    tmp = inDistriList(1,:);
                    tmpInList = [tmp.pi_x tmp.b_x tmp.ba_x];
                    
                    tmp = oDistriList(1,:);
                    tmpOutList = [tmp.pi_x tmp.b_x tmp.ba_x];
                    tmpOptInList = [tmp.opt_pi_x tmp.opt_b_x tmp.opt_ba_x];
                    ParaList = [ParaList;tmpInList tmpOutList];
                    OptParaList = [OptParaList; tmpOptInList];
                end
                
                directionList = oGradList;
                inDistriList = oDistriList;
                inDeconv = oDeconv;
                
                for n_iter = 1:size(oDistriList,2)
                    DistriObj = oDistriList(1,n_iter);
                    
                    DistriObj = DistriObj.updatePrior();
                    if(n_iter ~= 10)
                        if(DistriObj.checkPrior())
                            DistriObj = DistriObj.reinitialPrior();
                        end
                    end
                    inDistriList(1,n_iter) = DistriObj;
                end
            end
            %%
        end
        step_len = 0.0;
        
        compNum = [5,5,5];
        %更新后对参数进行优化；
        clear oDistriList oDeconv;
        [oDistriList,oGradList,oDeconv,oEnvidence] = ObjEnvidence.calcEnvidence(inDistriList,directionList,inDeconv,state,step_len,compNum);
        if(IsParaCompare)
            tmp = inDistriList(1,:);
                tmpInList = [tmp.pi_x tmp.b_x tmp.ba_x];
                
                tmp = oDistriList(1,:);
                tmpOutList = [tmp.pi_x tmp.b_x tmp.ba_x];
                tmpOptInList = [tmp.opt_pi_x tmp.opt_b_x tmp.opt_ba_x];
                ParaList = [ParaList;tmpInList tmpOutList];
                OptParaList = [OptParaList; tmpOptInList];
        end
      
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
    [oDistriList,oGradList,oDeconv,oEnvidence] = ObjEnvidence.calcEnvidence(inDistriList,directionList,inDeconv,state,step_len,[5 5 5]);
    if(IsParaCompare)
        tmp = inDistriList(1,:);
                tmpInList = [tmp.pi_x tmp.b_x tmp.ba_x];
                
                tmp = oDistriList(1,:);
                tmpOutList = [tmp.pi_x tmp.b_x tmp.ba_x];
                tmpOptInList = [tmp.opt_pi_x tmp.opt_b_x tmp.opt_ba_x];
                ParaList = [ParaList;tmpInList tmpOutList];
                OptParaList = [OptParaList; tmpOptInList];
    end
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
            [oDistriList,oGradList,oDeconv,oEnvidence] = ObjEnvidence.calcEnvidence(inDistriList,directionList,inDeconv,state,step_len,[5 5 5]);
            if(IsParaCompare)
               tmp = inDistriList(1,:);
                tmpInList = [tmp.pi_x tmp.b_x tmp.ba_x];
                
                tmp = oDistriList(1,:);
                tmpOutList = [tmp.pi_x tmp.b_x tmp.ba_x];
                tmpOptInList = [tmp.opt_pi_x tmp.opt_b_x tmp.opt_ba_x];
                ParaList = [ParaList;tmpInList tmpOutList];
                OptParaList = [OptParaList; tmpOptInList];
            end
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

save(ParaFileName,'ParaList','OptParaList');
D_log = D_log(1:iter);
gamma_log = gamma_log(1,1:iter);

out_mmu = inDistriList(1,1).mx
out_mx = reshape(inDistriList(1,2).mx,64,64);
out_me = reshape(inDistriList(1,3).mx,7,7);

mD = real(ifft2(fft2(out_mx).*fft2(out_me,64,64))) + out_mmu;

figure;
subplot(2,3,1);
imagesc(out_mx)
axis square
colorbar
title('Recovered Source');
subplot(2,3,2)
imagesc(out_me)
axis square
colorbar
title('Recovered Filter');
subplot(2,3,3)
imagesc(mD)
axis square
colorbar
title('Reconstructed Data');
subplot(2,3,4);
imagesc(D)
axis square
colorbar
title('Data');
subplot(2,3,5)
imagesc(me)
axis square
colorbar
title('Ori Me');
subplot(2,3,6)
imagesc(mx_org)
axis square
colorbar
title ('Ori Mx');


