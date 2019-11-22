function [outGrad,outPSF] = InitialPara(org_grad,new_grad,new_psf,IsStructure)

%¿‡À∆DFLEARNING£¨deconv÷––ﬁ∏ƒrenew£ª
D = org_grad;
mx = new_grad;
me = new_psf./sum(new_psf(:));

[mx_est,me_est] = ParaLearning(mx,me,D,IsStructure);

outPSF = me_est;

outGrad = mx_est;