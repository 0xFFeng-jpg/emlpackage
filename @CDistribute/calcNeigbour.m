function [st_mx2,sm_mx2] = calcNeigbour(obj,mx2)
%%
%first to get the 2-order neigbourhood 
m = obj.row;
n = obj.column;
matrix_m = reshape(mx2,m,n);

down_m = zeros(m,n);
down_m(1:m-1,:) = matrix_m(2:m,:);

up_m = zeros(m,n);
up_m(2:m,:) = matrix_m(1:m-1,:);

left_m = zeros(m,n);
left_m(:,2:n) = matrix_m(:,1:n-1);

right_m = zeros(m,n);
right_m(:,1:n-1) = matrix_m(:,2:n);

epslon = 1e-2;

%为了速度，我们选用4邻域；
v1 = abs(matrix_m - down_m) < epslon;
v2 = abs(matrix_m - up_m) < epslon;
v3 = abs(matrix_m - left_m) < epslon;
v4 = abs(matrix_m - right_m) < epslon;

st_mx2 = v1(:) + v2(:) + v3(:) + v4(:);
st_mx2 = st_mx2.*(st_mx2 < 4);

sm_mx2 = [];