function MultiScaleDemo4

%% 2010.09.26;
%��߶ȸ�ԭ���⣻

PSFShape = getPSFfromImage;
InitPSF = double(PSFShape);
InitPSF = InitPSF./sum(InitPSF(:));% normlization
figure;imshow(InitPSF,[]);
[mPSF,nPSF,lPSF] = size(InitPSF);

%% Second: get the Image;
%����ͼ��rgb��ʽ��
[filename, pathname] = uigetfile({'*.bmp','BMP�ļ�(*.bmp)';'*.jpg', 'JPEG�ļ�(*.jpg)';'*.png','PNG�ļ�(*.png)';'*.tif','Tif'});
if(filename == 0), return, end
filename = [pathname filename];
InRGBget1 = imread(filename);
%����ͨ�����д����ʱ����ʱ����ʹһ��ͨ����ͼ��������ͨ��ͼ��һ�£��Ӷ�������������ƫ�
Pitch = getPitch(InRGBget1);
InRGBget = Pitch;
figure;imshow(Pitch,[]);

InRGB = rgb2gray(InRGBget);
InRGB = (double(InRGB).^2.2)/(256.^(2.2-1));


[mInRGB,nInRGB,lInRGB] = size(InRGB);
OutRGB= zeros(size(InRGB));

RESIZE_STEP = 1.4142;%�Խ��߱�����
RESIZE_MODE = 'matlab_bilinear';
RESIZE_OBJ = 'grad'; % grad �������ִ������originΪԭʼͼ��gradΪ�ݶ�����
UPSAMPLE_MODE = 'matlab';
CENTER_BLUR = 0;

for nChannel = 1:lInRGB;
    
    %% Third: Extend the Image;
    % ��չͼ�������߽�����壻
    if(mod(mPSF,2))
        mEdge = 2*mPSF;
    else
        mEdge = 2*mPSF;
    end
    
    mEdge = 0
    
    if(mod(nPSF,2))
        nEdge = 2*nPSF;
    else
        nEdge = 2*nPSF;
    end
    
    nEdge = 0;
    
    InExtend = EdgePadArray(InRGB(:,:,nChannel),mEdge,nEdge);

    [mExtend,nExtend,lExtend] = size(InExtend);
    OutExtend = zeros(mExtend,nExtend,lExtend);
    
    %% Forth: Mark the Valid Region;
    % ��ȡ��Ч����Ŀ���Ƕ���չ�߽紦�������ദ��
    ValidRegion = zeros(mExtend,nExtend);
    half_mEdge= mEdge/2;
    half_nEdge = nEdge/2;
    ValidRegion(half_mEdge+1:half_mEdge+mInRGB,half_nEdge+1:half_nEdge+nInRGB) = ones(mInRGB,nInRGB);
    
    %% Fifth: Get the Smooth Region;
    SmoothRegion = getSmoothRegion(InExtend,InitPSF);
    
    %% MultiScale: 1;
    % ��PSF����Ҫ�ֽ�ĳ߶ȣ�
    
    k_max = max(mPSF,nPSF);
    k_min = min(mPSF,nPSF);
    NumScale = floor(log(k_min)/log(RESIZE_STEP))-1;
    
     %% MultiScale: 2;
    % �ֽ�PSF��
    PSFScale{NumScale} = InitPSF;
    for s = 2:NumScale
        
        dims = size(PSFScale{NumScale}) * (1/RESIZE_STEP)^(s-1);
        dims = dims + (1-mod(dims,2)); %% make odd size     
        if (min(dims)<1)
            
            %% Blur kernel first, then resize (since it will use nearest neighbour)
            h = fspecial('gaussian',dims,1);
            PSFScale{NumScale-s+1} = imresize(conv2(InitPSF,h),dims,'nearest');
            
        else
            %% ���ý�����ʽ��������
            if strcmp(RESIZE_MODE,'matlab_nearest')
                PSFScale{NumScale-s+1} = imresize(PSFScale{NumScale-s+2},dims,'nearest');
            elseif  strcmp(RESIZE_MODE,'matlab_bilinear')
                PSFScale{NumScale-s+1} = imresize(PSFScale{NumScale-s+2},dims,'bilinear');
            elseif   strcmp(RESIZE_MODE,'matlab_bicubic')
                PSFScale{NumScale-s+1} = imresize(PSFScale{NumScale-s+2},dims,'bicubic');
            elseif   strcmp(RESIZE_MODE,'edge_Interpolation')
                PSFScale{NumScale-s+1} = imresize(PSFScale{NumScale-s+2},dims,'bicubic');
            else
                error foo
            end
            
        end
        %% ��һ����
        PSFScale{NumScale-s+1} = PSFScale{NumScale-s+1} / sum(PSFScale{NumScale-s+1}(:));
    end
    
    %% MultiScale: 3;
    % ͼ���ݶȷֽ⣻ֱ�����ţ�
    obsImScale{NumScale}   = InExtend;
    smthRgScale{NumScale}  = SmoothRegion;
    validRgScale{NumScale} = ValidRegion;
    
%     
%     obsImScale_dx{NumScale} = real(ifft2(fft2(InExtend).*psf2otf([1 -1],[mExtend,nExtend]),mExtend,nExtend));
%     obsImScale_dy{NumScale} = real(ifft2(fft2(InExtend).*psf2otf([1 -1]',[mExtend,nExtend]),mExtend,nExtend));
    
    obsIm_dx = conv2(InExtend,[1 -1],'valid');
    obsIm_dy = conv2(InExtend,[1 -1]','valid');
    
    yy = min(size(obsIm_dx,1),size(obsIm_dy,1));
    xx = min(size(obsIm_dx,2),size(obsIm_dy,2));
    
    obsImScale_dx{NumScale} = obsIm_dx(1:yy,1:xx);
    obsImScale_dy{NumScale} = obsIm_dy(1:yy,1:xx);
    
    obsImScale_grad{NumScale} = [obsImScale_dx{NumScale},obsImScale_dy{NumScale}];
    
    for s = 2:NumScale
        
%         dims = size(obsImScale{NumScale}) * (1/RESIZE_STEP)^(s-1);
%         dims = dims + (1-mod(dims,2)); %% make odd size
        dims = (1/RESIZE_STEP)^(s-1);
        
        
        if (min(dims)<0)
            
            %% Blur kernel first, then resize (since it will use nearest neighbour)
            error('the size is small');
            
        else
            

            if strcmp(RESIZE_MODE,'matlab_nearest')%�˴�Ӧ�򻯣�����ʱ�䡣
                obsImScale{NumScale-s+1} = imresize(obsImScale{NumScale},dims,'nearest');
                smthRgScale{NumScale-s+1} = imresize(smthRgScale{NumScale},dims,'nearest');
                validRgScale{NumScale-s+1} = imresize(validRgScale{NumScale},dims,'nearest');
                
                obsImScale_dx{NumScale-s+1} = imresize(obsImScale_dx{NumScale},dims,'nearest');
                obsImScale_dy{NumScale-s+1} = imresize(obsImScale_dy{NumScale},dims,'nearest');
                
            elseif  strcmp(RESIZE_MODE,'matlab_bilinear')
                obsImScale{NumScale-s+1} = imresize(obsImScale{NumScale},dims,'bilinear');
                smthRgScale{NumScale-s+1} = imresize(smthRgScale{NumScale},dims,'bilinear');
                validRgScale{NumScale-s+1} = imresize(validRgScale{NumScale},dims,'bilinear');
                
                obsImScale_dx{NumScale-s+1} = imresize(obsImScale_dx{NumScale},(1/RESIZE_STEP)^(s-1),'bilinear');
                obsImScale_dy{NumScale-s+1} = imresize(obsImScale_dy{NumScale},(1/RESIZE_STEP)^(s-1),'bilinear');
            elseif   strcmp(RESIZE_MODE,'matlab_bicubic')
                obsImScale{NumScale-s+1} = imresize(obsImScale{NumScale},dims,'bicubic');
                smthRgScale{NumScale-s+1} = imresize(smthRgScale{NumScale},dims,'bicubic');
                validRgScale{NumScale-s+1} = imresize(validRgScale{NumScale},dims,'bicubic');
                
                obsImScale_dx{NumScale-s+1} = imresize(obsImScale_dx{NumScale},dims,'bicubic');
                obsImScale_dy{NumScale-s+1} = imresize(obsImScale_dy{NumScale},dims,'bicubic');
            elseif   strcmp(RESIZE_MODE,'edge_Interpolation')
                obsImScale{NumScale-s+1} = AdapInterpolation(obsImScale{NumScale},dims);
                smthRgScale{NumScale-s+1} = AdapInterpolation(smthRgScale{NumScale},dims);
                validRgScale{NumScale-s+1} = AdapInterpolation(validRgScale{NumScale},dims);
                
                obsImScale_dx{NumScale-s+1} = AdapInterpolation(obsImScale_dx{NumScale},dims);
                obsImScale_dy{NumScale-s+1} = AdapInterpolation(obsImScale_dy{NumScale},dims);
            else
                error foo
            end
            obsImScale_grad{NumScale-s+1} = [obsImScale_dx{NumScale-s+1},obsImScale_dy{NumScale-s+1}];
        end
    end
    
    if(1)
        for s = 1:NumScale
            obsImScale_grad_old{s} = obsImScale_grad{s};
            db = delta_kernel(size(PSFScale{s},1));
            obsImScale_grad{s} = real(ifft2(fft2(obsImScale_grad{s}).*fft2(db,size(obsImScale_grad{s},1),size(obsImScale_grad{s},2))));
        end
    end
    
    %% MultiScale: 4;
    % ���ÿ���߶Ƚ��м��㣻���ÿ���߶ȵĴ�����������������չ
    for s = 1:NumScale
        if strcmp(RESIZE_OBJ,'grad')
            if(s == 1)
                [OutImScale_grad{s},OutPSFScale{s}] = SingleScaleProcess6(obsImScale_grad{s},...
                   obsImScale_grad{s}, PSFScale{s});     
            else
                new_kernel = imresize(OutPSFScale{s-1},size(PSFScale{s}),'nearest');
                newImScale_grad = imresize(OutImScale_grad{s-1},size(obsImScale_grad{s}),'bilinear');
                [m_k,n_k] = size(PSFScale{s});
                [m_im,n_im] = size(obsImScale_grad{s});
                
                [newImScale_grad,new_kernel] = move_level(OutImScale_grad{s-1},OutPSFScale{s-1},...
                    m_k,n_k,m_im,n_im,UPSAMPLE_MODE,RESIZE_STEP,CENTER_BLUR);
                
                [OutImScale_grad{s},OutPSFScale{s}] = SingleScaleProcess6(obsImScale_grad{s},...
                    newImScale_grad,new_kernel);
                
            end
        else
            error foo;
        end
        
        %%
        %write to file;
        time_v = clock;
        OutDirName = fullfile(['y' int2str(time_v(1)) 'mm' int2str(time_v(2)) 'd' int2str(time_v(3))]);
        SingleName = fullfile(['S_' int2str(s) 'h' int2str(time_v(4)) 'm' int2str(time_v(5)) 's' int2str(round(time_v(6)))]);
        %write PSF
        tmpPSFScale = getNormImg(OutPSFScale{s});
        mkdir(['tmpfdir\' OutDirName]);
        imwrite(uint8(tmpPSFScale*200),['tmpfdir\' OutDirName '\' SingleName '.jpg'],'jpg');
        %write Result
        tmpOutdx = getNormImg(OutImScale_grad{s});
        mkdir(['tmpgrad\' OutDirName]);
        imwrite(uint8(tmpOutdx*200),['tmpgrad\' OutDirName '\' SingleName '.jpg'],'jpg');
        
       
        
        tmpOridx = getNormImg(obsImScale_grad{s});
        mkdir(['tmpOrigrad\' OutDirName]);
        imwrite(uint8(tmpOridx*200),['tmpOrigrad\' OutDirName '\' SingleName '.jpg'],'jpg');
        
      
        
        
       
    end
    %��������µ�PSF���ݶ�ͼ���ؽ�ͼ��
    PSF = OutPSFScale{NumScale};
    [OutImage,OutPSF] = deconvDF(InRGBget1,PSF);
end
time_v = clock;
%imwrite(uint8(OutRgb),['Rdir\' 'h' int2str(time_v(4)) 'm' int2str(time_v(5)) 's' int2str(round(time_v(6))) '.bmp'],'bmp');
OutDirName = fullfile(['y' int2str(time_v(1)) 'mm' int2str(time_v(2)) 'd' int2str(time_v(3))]);
SingleName = fullfile(['h' int2str(time_v(4)) 'm' int2str(time_v(5)) 's' int2str(round(time_v(6)))]);
mkdir(['Rdir\' OutDirName]);
imwrite(uint8(OutImage),['Rdir\' OutDirName '\' SingleName '.jpg'],'jpg');

mkdir(['fdir\' OutDirName]);
imwrite(uint8(OutPSF*200),['fdir\' OutDirName '\' SingleName '.jpg'],'jpg');
    
    
    