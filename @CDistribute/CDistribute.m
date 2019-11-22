classdef CDistribute
    %�ֲ��ࣻ
    %�����������ݵķֲ���
    %�Ե�ģ�ͣ��������������б�ʾ��
    properties
        x1 %�����������ݵľ�ֵ��EM�㷨��������������ֵ��һ�ף�
        x2 %�����������ݵķ����������ĳ˻������Ƶ���ʽ��أ�
        opt_cx1
        opt_cx2
        row
        column
    end
    
    properties
        pi_x %��Ϸֲ��ķֲ����ӣ�
        b_x  %ϡ��ֲ��У���������ֲ��Ĳ���,gamma�ֲ��е�b��
        ba_x %ϡ��ֲ��У���������ֲ��Ĳ�����gamma�ֲ��е�b/a;
        opt_pi_x
        opt_b_x
        opt_ba_x
    end
    
    properties
        DistriType
        ObjType
        IsManual
        IsStructure
        CompType %���������ĸ���������С�ڵ���������������ÿ�����������������������ͬʱ����ȫ����
    end
    
    properties
        c_log_lamda = []; % log(��Ϸֲ���������
        CataNum;
        CataDegree;
    end
    
    properties
        Hx
        mx
        mx2
    end
    
    properties
        KLDistance
    end
    
    methods
        function obj = CDistribute(init_x1,init_x2,init_pi,init_b,init_ba,ditriType,IsManual,objType,row,column,IsStructure)
            obj.x1 = init_x1;
            obj.x2 = init_x2;
            obj.pi_x = init_pi;
            obj.b_x = init_b;
            obj.ba_x = init_ba;
            obj.DistriType = ditriType;
            obj.IsManual = IsManual;
            obj.ObjType = objType;
            obj.row = row;
            obj.column = column;
            obj.IsStructure = IsStructure;
        end     
    end
    

    
    methods
        obj = calcKLDistance(obj); %�����ʼ�ֲ����Ż����KL���롣
        obj = opitmizeParabyPrior(obj); %�Ż���������ö�����Ż���
        obj = calcLamda(obj);   %�����м����lamda;
        obj = calcExpection(obj);     %����x1,x2����༶������
        obj = calcOptCx(obj);
        
        re = checkPrior(obj);
        [t,sm] = calcNeigbour(obj,mx2);
        obj = reinitialPrior(obj);
        obj = updatePrior(obj);
    end
    
    methods
        [outDistri,outGrad] = renew(obj,inGrad,state,steplen,CompType);
    end
end
    