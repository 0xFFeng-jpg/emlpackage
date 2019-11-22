function mixturePlot(para,en,im,ker)
%
%en 是ensemble的组件数；
%im 是图像
%ker 是核；
paraUpdate = para.ParaList;
paraOpt = para.OptParaList;


compNum = [en im ker];
s_c = sum(compNum);

[m,n] = size(paraUpdate);

iterNum = m;

pi_old = [paraUpdate(m,1:s_c)];
pi_new = [paraUpdate(m,3*s_c+1:4*s_c)];
b_old = [paraUpdate(m,1*s_c+1:2*s_c)];
b_new = [paraUpdate(m,4*s_c+1:5*s_c)];
ba_old = [paraUpdate(m,2*s_c+1:3*s_c)];
ba_new = [paraUpdate(m,5*s_c+1:6*s_c)];

opt_pi = [paraOpt(m,1:s_c)];
opt_b = [paraOpt(m,1*s_c+1:2*s_c)];
opt_ba = [paraOpt(m,2*s_c+1:3*s_c)];

distriPlot(pi_old(en+1:en+im),zeros(1,im),ba_old(en+1:en+im));
%distriPlot(pi_old(en+im+1:en+im+ker),zeros(1,ker),ba_old(en+1+im:en+im+ker));
distriPlot(pi_new(en+1:en+im),zeros(1,im),ba_new(en+1:en+im));
distriPlot(opt_pi(en+1:en+im),zeros(1,im),opt_ba(en+1:en+im));