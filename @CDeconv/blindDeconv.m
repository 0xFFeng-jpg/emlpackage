function obj = blindDeconv(obj)
%FFT 模式
[m_k,n_k] = size(obj.OriK);
[m_D,n_D] = size(obj.OriIm);
m_im = m_D/2;
n_im = n_D/2;

WEIGHT = obj.ValidRegion;



mmu = obj.Ensemble.mx;
me = reshape(obj.KernelD.mx,m_k,n_k);
mx = reshape(obj.ImD.mx,m_im,n_im);

mmu2 = obj.Ensemble.mx2;
me2 = reshape(obj.KernelD.mx2,m_k,n_k);
mx2 = reshape(obj.ImD.mx2,m_im,n_im);

%在频域解决卷积；
fft_Dp = fft2(WEIGHT);

fft_mx = fft2(mx,m_D,n_D);
fft_mx2 = fft2(mx2,m_D,n_D);
fft_mx3 = fft2(mx.^2,m_D,n_D);


fft_me = fft2(me,m_D,n_D);
fft_me2 = fft2(me2,m_D,n_D);
fft_me3 = fft2(me.^2,m_D,n_D);

deconvIm = real(ifft2(fft_mx.*fft_me));
fft_deconvErr = fft2(WEIGHT.*(obj.OriIm - deconvIm - mmu),m_D,n_D);

data_points = sum(WEIGHT(:));

error = sum(sum(WEIGHT.*((obj.OriIm - deconvIm - mmu).^2 + real(ifft2(fft_me2.*fft_mx2 - fft_me3.*fft_mx3))))) + data_points*(mmu2 - mmu^2);

e1 = real(ifft2(fft_deconvErr.*conj(fft_mx)));
corr = real(ifft2(fft_Dp.*conj(fft_mx3)));
e1(1:m_k,1:n_k) = e1(1:m_k,1:n_k) + me.*corr(1:m_k,1:n_k);

e2 = real(ifft2(fft_Dp.*conj(fft_mx2)));
e2 = e2(1:m_k,1:n_k);


x1 = real(ifft2(fft_deconvErr.*conj(fft_me)));
corr = real(ifft2(fft_Dp.*conj(fft_me3)));

x1(1:m_im,1:n_im) = x1(1:m_im,1:n_im) + mx.*corr(1:m_im,1:n_im);

x2 = real(ifft2(fft_Dp.*conj(fft_me2)));
x2 = x2(1:m_im,1:n_im);

obj.im_x1 = reshape(x1(1:m_im,1:n_im),m_im*n_im,1);
obj.im_x2 = reshape(x2(1:m_im,1:n_im),m_im*n_im,1);

obj.kernel_x1 = reshape(e1(1:m_k,1:n_k),m_k*n_k,1);
obj.kernel_x2 = reshape(e2(1:m_k,1:n_k),m_k*n_k,1);

obj.ensemble_x1 = sum(sum(WEIGHT.*(obj.OriIm - deconvIm)));
obj.ensemble_x2 = data_points;

obj.error = error;
obj.data_points = data_points;
