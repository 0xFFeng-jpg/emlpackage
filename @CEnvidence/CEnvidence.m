classdef CEnvidence
    %��Ҫ���ڼ�������ѧϰ���̣�
    properties
        DitribList   %���ֲַ������б�
        GradList     %�����ݶȷֲ��б�
        DeconvObj    %�������������
    end
    
    properties
        OriIm
        OriKer
    end
    
    properties
        state   %״̬
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
    
    methods %�����ӿ�
        [oDitriList,oGradList,oDeconvObj,oEnvidence] = calcEnvidence(obj,DitriList,GradList,DeconvObj,state,steplen,CompType);
        [oDitriList,oGradList,oDeconvObj,oEnvidence] = calcParaEnvidence(obj,DitriList,GradList,DeconvObj,state,steplen,CompType);
    end
    
    methods %˽�нӿ�
       
    end
    
end
        
        