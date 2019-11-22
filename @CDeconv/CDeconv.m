classdef CDeconv
    %��Ϊ��Ϊһ���������̣��������ǽ�����������д���
    %�����ֲ�
    properties
        im_x1  %ͼ���ݶȾ�ֵ�뷽��˻�
        im_x2  %����
        
        kernel_x1  %����˵ľ�ֵ�뷽��ĳ˻�
        kernel_x2  %����˵ķ���
        
        ensemble_x1
        ensemble_x2
        
        
        error      %�ؽ����
        
        data_points  %��Чͼ����������
    end
    
    
    properties
        ImD                 %class CDistribute
        KernelD             %class CDistribute
        Ensemble
        
        OriIm             %ԭʼͼ��
        OriK              %ԭʼ�����
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
    
    %�����
    methods 
        obj = blindDeconv(obj);
        obj = calcSigmaDx(obj,isSigmaUpdate); %���������������ӦKL���룻
        obj = updateExpection(obj); %�����ź������˵ľ�ֵ�뷽��ȣ�
    end
    
    methods 
        [outDeconv,outGradList] = renew(obj,inDistriList,inGradList,state);
        [outDeconv,outGradList] = renewPara(obj,inDistriList,inGradList,state);
    end
end