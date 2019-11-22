function distriPlot(pi_x,miu,invsigma)

%plot mixture distribution;
x_min = -20;
x_max = 20;
x_step = 0.1;

x_v = [x_min:x_step:x_max];
tmp = [];

[dim,component] = size(pi_x);
for iter = 1:component
    ctmp = pi_x(1,iter).*Posibality(x_v,'Gaussian',miu(1,iter),invsigma(1,iter));
    %ctmp = 1.*Posibality(x_v,'Gaussian',miu(1,iter),invsigma(1,iter));
    tmp = [tmp;ctmp];
end

out = sum(tmp,1);
figure;plot(x_v,[tmp;out]);


