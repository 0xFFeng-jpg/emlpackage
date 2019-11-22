function shapediag = rectdiag(shape,direction)
%
%生成对角线，如果direction为1，则为主，否则为辅
%

[m,n] = size(shape);
if(n == 0 || m == 0)
    disp('Can not calc the diag');
    return;
end

shapediag = zeros(m,n);

k = double(m)./double(n);

if(direction)
    if(k < 1)%说明为扁形；
        for j = 1:n
            i = floor(double(j) * k)+1;
            if(i > 0 && i <= m)
                shapediag(i,j) = 1;
            end
        end
    else
        for i = 1:m
            j = ceil(double(i) / k);
            if( j > 0 && j <= n)
                shapediag(i,j) = 1;
            end
        end
    end
else
    if(k < 1)
        for j = 1:n
            i = floor(double(j) * (-k)) + m+1;
            if( i > 0 && i <= m)
                shapediag(i,j) = 1;
            end
        end
    else
        for i = 1:m
            j = floor(double(i)/(-k)) + n + 1;
            if( j > 0 && j <= n)
                shapediag(i,j) = 1;
            end
        end
    end
end

