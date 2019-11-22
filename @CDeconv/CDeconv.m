classdef CDeconv
    %因为其为一个交互过程，所以我们将其提出来进行处理
    %两个分布
    properties
        im_x1  %图像梯度均值与方差乘积
        im_x2  %方差
        
        kernel_x1  %卷积核的均值与方差的乘积
        kernel_x2  %卷积核的方差
        
        ensemble_x1
        ensemble_x2
        
        
        error      %重建误差
        
        data_points  %有效图像像素数；
    end
    
    
    properties
        ImD                 %class CDistribute
        KernelD             %class CDistribute
        Ensemble
        
        OriIm             %原始图像；
        OriK              %原始卷积核
        ValidRegion
    end
    
    properties
        b_sigma
        ba_sigma = 1;
        opt_ba_sigma
    end
    
    properties
        D_x
    end
    
    methods
        function obj = CDeconv(DistriList,OriginalIm,ValidRegion,OriginalK)
            obj.ImD= DistriList(1,2);
            obj.KernelD = DistriList(1,3);
            obj.Ensemble = DistriList(1,1);
            obj.OriIm = OriginalIm;
            obj.OriK  = OriginalK;
            obj.ValidRegion = ValidRegion;
        end
    end
    
    %卷积；
    methods 
        obj = blindDeconv(obj);
        obj = calcSigmaDx(obj,isSigmaUpdate); %计算噪声方差和相应KL距离；
        obj = updateExpection(obj); %更新信号与卷积核的均值与方差等；
    end
    
    methods 
        [outDeconv,outGradList] = renew(obj,inDistriList,inGradList,state);
        [outDeconv,outGradList] = renewPara(obj,inDistriList,inGradList,state);
    end
end