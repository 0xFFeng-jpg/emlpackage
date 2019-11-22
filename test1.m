clear all;
distType = 'Gaussian';%'Laplacian';
% a = CDistribute([1:100],[1:100],[1:4],[1:4],[1:4]);
% %a = CDistribute([1:100],[1:100],1,1,1);
% a = calcExpection(a,distType);
% a = calcLamda(a,distType);
% a = optimizeParaByPrior(a,distType);
% a = calcKLDistance(a);
% a = calcOptCx(a,distType);
x = ((-3:3)'*ones(1,7));
y = (ones(7,1)*(-3:3));
me = exp(-(x.^2 + x.*y + y.^2)/2/1^2);
mx = zeros(64);
mx(ceil(rand(1,20)*64*64)) = -log(rand(1,20));

D = real(ifft2(fft2(mx).*psf2otf(me,size(mx))));
%D = D/sqrt(mean(mean(D.^2))) + 1 + randn(64,64)*1e-2;

im_x1  = 1e4*D(:);
im_x2 = 1e4*ones(size(im_x1));

kernel_x1 = 1e4*me(:);
kernel_x2 = 1e4*ones(size(kernel_x1));

Im = CDistribute(im_x1,im_x2,1*ones(1,2),1*ones(1,2),1*ones(1,2));
Ker = CDistribute(kernel_x1,kernel_x2,1*ones(1,2),1*ones(1,2),1*ones(1,2));

ImdistType = 'Gaussian';
Im = Im.calcExpection(ImdistType);
Im = Im.calcLamda(ImdistType);
Im = Im.optimizeParaByPrior(ImdistType);
Im = Im.calcKLDistance();
Im = Im.calcOptCx(ImdistType);

KdistType = 'Laplacian';
Ker = Ker.calcExpection(KdistType);
Ker = Ker.calcLamda(KdistType);
Ker = Ker.optimizeParaByPrior(KdistType);
Ker = Ker.calcKLDistance();
Ker = Ker.calcOptCx(KdistType);

bDeconv = CDeconv(Im,Ker,D,ones(size(me)));
bDeconv = bDeconv.blindDeconv;
bDeconv = bDeconv.calcSigmaDx(0);
bDeconv = bDeconv.updateExpection;


