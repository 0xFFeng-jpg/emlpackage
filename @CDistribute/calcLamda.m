function obj = calcLamda(obj)
%������ϵ��log(lamda);
[dim_vector,component] = size(obj.pi_x);%dim_vector ��Զ�ά�������ã���ͼ�������ݲ����ǣ�
num_varients = size(obj.x1,1);
t_pi_x = obj.pi_x;
t_b_x = obj.b_x;
t_ba_x = obj.ba_x;

t_mx = obj.mx;
t_mx2 = obj.mx2;
%��ͼ������ƶ��������ܹ����д����ĸ�����
if(strcmp(obj.ObjType,'Image') & obj.IsStructure)
    [st_mx2,sm_mx2] = obj.calcNeigbour(t_mx2);
end


c_log_lamda = zeros(num_varients,component);%��ʼ��0���������Ҳ������

distributeType = obj.DistriType;
if(component <= 1)
    %disp('No mixture,lamda is zeros');
else
    if (strcmp(obj.ObjType,'Image') & obj.IsStructure)
        switch distributeType
            case 'Gaussian'
                for alpha = 1:component
                    c_log_lamda(:,alpha) = log(t_pi_x(1,alpha)) - 0.5/t_pi_x(1,alpha)+...
                        0.5*log(t_ba_x(1,alpha))-0.25/t_b_x(1,alpha) - 0.5*t_mx2(:)*t_ba_x(1,alpha) - st_mx2(:)./t_ba_x(1,alpha);
                end
            case 'Laplacian' %����������һ�£�
                for alpha = 1:component
                    c_log_lamda(:,alpha) = log(t_pi_x(1,alpha)) - 0.5/t_pi_x(1,alpha)+...
                        log(t_ba_x(1,alpha)) - 0.5/t_b_x(1,alpha) - t_mx(:)*t_ba_x(1,alpha);
                end
            case 'Exponential'
                for alpha = 1:component
                    c_log_lamda(:,alpha) = log(t_pi_x(1,alpha)) - 0.5/t_pi_x(1,alpha) + ...
                        0.5*log(t_ba_x(1,alpha)) - 0.25/t_ba_x(1,alpha) - 0.5*abs(t_mx(:))*t_ba_x(1,alpha);
                end
            otherwise
                error('No such distribution');
        end
    else
        switch distributeType
            case 'Gaussian'
                for alpha = 1:component
                    c_log_lamda(:,alpha) = log(t_pi_x(1,alpha)) - 0.5/t_pi_x(1,alpha)+...
                        0.5*log(t_ba_x(1,alpha))-0.25/t_b_x(1,alpha) - 0.5*t_mx2(:)*t_ba_x(1,alpha);
                end
            case 'Laplacian' %����������һ�£�
                for alpha = 1:component
                    c_log_lamda(:,alpha) = log(t_pi_x(1,alpha)) - 0.5/t_pi_x(1,alpha)+...
                        log(t_ba_x(1,alpha)) - 0.5/t_b_x(1,alpha) - t_mx(:)*t_ba_x(1,alpha);
                end
            case 'Exponential'
                for alpha = 1:component
                    c_log_lamda(:,alpha) = log(t_pi_x(1,alpha)) - 0.5/t_pi_x(1,alpha) + ...
                        0.5*log(t_ba_x(1,alpha)) - 0.25/t_ba_x(1,alpha) - 0.5*abs(t_mx(:))*t_ba_x(1,alpha);
                end
            otherwise
                error('No such distribution');
        end
    end
    %��һ��Լ��lamda;�˴���������˼����
    max_c_log_lamda = max(c_log_lamda,[],2);
    for alpha = 1:component
        c_log_lamda(:,alpha) = c_log_lamda(:,alpha) - max_c_log_lamda;
    end
    log_sum_lamda = log(sum(exp(c_log_lamda),2));
    for alpha = 1:component
        c_log_lamda(:,alpha) = c_log_lamda(:,alpha) - log_sum_lamda;
    end    
end

obj.c_log_lamda = c_log_lamda;

obj.CataNum = sum(exp(c_log_lamda),1);

I1 = (obj.CataNum - size(c_log_lamda,1).*(obj.pi_x./sum(obj.pi_x))).^2;
I2 = size(c_log_lamda,1).*(obj.pi_x./sum(obj.pi_x));

obj.CataDegree = sum(I1./I2);



