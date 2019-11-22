classdef CDistribute
    %分布类；
    %用于描述数据的分布；
    %对单模型，用两个参数进行表示；
    properties
        x1 %用以描述数据的均值，EM算法就是算期望，均值是一阶；
        x2 %用以描述数据的方差与期望的乘积，与推导公式相关；
        opt_cx1
        opt_cx2
        row
        column
    end
    
    properties
        pi_x %混合分布的分布因子：
        b_x  %稀疏分布中，描述方差分布的参数,gamma分布中的b；
        ba_x %稀疏分布中，描述方差分布的参数，gamma分布中的b/a;
        opt_pi_x
        opt_b_x
        opt_ba_x
    end
    
    properties
        DistriType
        ObjType
        IsManual
        IsStructure
        CompType %描述更新哪个组件，如果小于等于组件数，则更新每个组件，如果大于组件数，则同时更新全部；
    end
    
    properties
        c_log_lamda = []; % log(混合分布）期望；
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
        obj = calcKLDistance(obj); %计算初始分布与优化后的KL距离。
        obj = opitmizeParabyPrior(obj); %优化参数，获得对外的优化项
        obj = calcLamda(obj);   %计算中间参数lamda;
        obj = calcExpection(obj);     %根据x1,x2计算多级期望；
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
    