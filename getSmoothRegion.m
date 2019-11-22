function SmoothRegion = getSmoothRegion(ingray,initpsf)
%通过变化区域来获取平滑区域；
%coder: Fred 
%date : 2009.11.3
%input: ingray is image only one channel and initpsf is psfshape;
%output: SmoothRegion is which unchanged part of the image;
th = 5; %变换的阈值；

[m,n] = size(ingray);
[k,l] = size(initpsf);



SmoothRegion = ones(m,n);

half_k = ceil(k/2.0);
half_l = ceil(l/2.0);

%对图像区域取block;
blk_std = 1000;
for i = 1+half_k:m-half_k
    for j = 1+half_l:n-half_l
        blkm = ingray(i-half_k:i+half_k,j-half_l:j+half_l);
        blk = blkm(:);
       % blk_std = sqrt(sum((blk-mean(blk)).^2)/size(blk,1)/((size(blk,1)-1)));
        blk_std = std2(blkm);
        if(blk_std > th)
            SmoothRegion(i,j) = 0.0;
        end
    end
end