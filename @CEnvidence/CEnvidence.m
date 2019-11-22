classdef CEnvidence
    %主要用于计算整体学习过程；
    properties
        DitribList   %几种分布对象列表；
        GradList     %几种梯度分布列表；
        DeconvObj    %卷积及噪声对象；
    end
    
    properties
        OriIm
        OriKer
    end
    
    properties
        state   %状态
    end
    
    properties
        envidenceList
    end
    
    methods
        function obj = CEnvidence(DitriList,GradList,DeconvObj,state,OriIm,OriKer)
            obj.DitribList = DitriList;
            obj.GradList = GradList;
            obj.DeconvObj = DeconvObj;
            obj.state = state;
            
            obj.OriIm = OriIm;
            obj.OriKer = OriKer;
        end
    end
    
    methods %公共接口
        [oDitriList,oGradList,oDeconvObj,oEnvidence] = calcEnvidence(obj,DitriList,GradList,DeconvObj,state,steplen,CompType);
        [oDitriList,oGradList,oDeconvObj,oEnvidence] = calcParaEnvidence(obj,DitriList,GradList,DeconvObj,state,steplen,CompType);
    end
    
    methods %私有接口
       
    end
    
end
        
        