function a = paraPlot(para,en,im,ker)
%
%en 是ensemble的组件数；
%im 是图像
%ker 是核；
paraUpdate = para.ParaList;
paraOpt = para.OptParaList;


compNum = [en im ker];
s_c = sum(compNum);
%pi 0; b 1;ba 2;
pi = [paraUpdate(:,1:s_c) paraUpdate(:,3*s_c+1:4*s_c)];
b = [paraUpdate(:,1*s_c+1:2*s_c) paraUpdate(:,4*s_c+1:5*s_c)];
ba = [paraUpdate(:,2*s_c+1:3*s_c) paraUpdate(:,5*s_c+1:6*s_c)];

opt_pi = [paraOpt(:,1:s_c)];
opt_b = [paraOpt(:,1*s_c+1:2*s_c)];
opt_ba = [paraOpt(:,2*s_c+1:3*s_c)];

figure;
subplot(3,3,1);
plot(pi(:,1:en));
title('EnPiOld');
subplot(3,3,4)
plot(pi(:,s_c+1:s_c+en));
title('EnPiNew');
subplot(3,3,2)
plot(pi(:,en+1:en+im));
title('ImPiOld');
subplot(3,3,5)
plot(pi(:,en+1+s_c:en+im+s_c));
title('ImPiNew');
subplot(3,3,3);
plot(pi(:,en+im+1:en+im+ker));
title('KerPiNew');
subplot(3,3,6);
plot(pi(:,en+im+1+s_c:en+im+ker+s_c));
title('KerPiOld');
subplot(3,3,7);
plot(opt_pi(:,1:en));
title('EnOptPiOld');
subplot(3,3,8)
plot(opt_pi(:,en+1:en+im));
title('ImOptPiOld');
subplot(3,3,9);
plot(opt_pi(:,en+im+1:en+im+ker));
title('KerOptPiNew');

figure;
subplot(3,3,1);
plot(b(:,1:en));
title('EnBOld');
subplot(3,3,4)
plot(b(:,s_c+1:s_c+en));
title('EnBNew');
subplot(3,3,2)
plot(b(:,en+1:en+im));
title('ImBOld');
subplot(3,3,5)
plot(b(:,en+1+s_c:en+im+s_c));
title('ImBNew');
subplot(3,3,3);
plot(b(:,en+im+1:en+im+ker));
title('KerBNew');
subplot(3,3,6);
plot(b(:,en+im+1+s_c:en+im+ker+s_c));
title('KerBOld');
subplot(3,3,7);
plot(opt_b(:,1:en));
title('EnOptBOld');
subplot(3,3,8)
plot(opt_b(:,en+1:en+im));
title('ImOptBOld');
subplot(3,3,9);
plot(opt_b(:,en+im+1:en+im+ker));
title('KerOptBNew');

figure;
subplot(3,3,1);
plot(ba(:,1:en));
title('EnBaOld');
subplot(3,3,4)
plot(ba(:,s_c+1:s_c+en));
title('EnBaNew');
subplot(3,3,2)
plot(ba(:,en+1:en+im));
title('ImBaOld');
subplot(3,3,5)
plot(ba(:,en+1+s_c:en+im+s_c));
title('ImBaNew');
subplot(3,3,3);
plot(ba(:,en+im+1:en+im+ker));
title('KerBaNew');
subplot(3,3,6);
plot(ba(:,en+im+1+s_c:en+im+ker+s_c));
title('KerBaOld');
subplot(3,3,7);
plot(opt_ba(:,1:en));
title('EnOptBaOld');
subplot(3,3,8)
plot(opt_ba(:,en+1:en+im));
title('ImOptBaOld');
subplot(3,3,9);
plot(opt_ba(:,en+im+1:en+im+ker));
title('KerOptBaNew');


