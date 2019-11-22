function [OutIm,OutPSF] = deconvDF(Im,PSF)

THRESHOLD = 10;
LUCY_ITS = 20;

blur_kernel = PSF/sum(PSF(:));

threshold = max(blur_kernel(:))/THRESHOLD;
z = find(blur_kernel(:) < threshold);
blur_kernel(z) = 0;

blur_kernel = blur_kernel / sum(blur_kernel(:));

obs_im_gam = edgetaper(Im,blur_kernel);

out = deconvlucy(obs_im_gam,blur_kernel,LUCY_ITS);

out = double(out);

out = out - min(out(:));
out = out / max(out(:));

out = histmatch(out,uint8(Im));

OutIm = out;
OutPSF = blur_kernel;