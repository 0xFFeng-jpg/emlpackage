function InExtend = EdgePadArray(In,edgeSizeM,edgeSizeN)
%通过对图像边界宽度的定义，来扩展图像；
[mim,nim,lin] = size(In);

mextend = mim + edgeSizeM;
nextend = nim + edgeSizeN;

halfedgesizem = edgeSizeM/2;
halfedgesizen = edgeSizeN/2;

InExtend = zeros(mextend,nextend,lin);

%InExtend = padarray(In,[halfedgesizem,halfedgesizen],'circular','both');

%中心位置
%中心位置为实际图像；
InExtend(halfedgesizem+1:halfedgesizem+mim,halfedgesizen+1:halfedgesizen+nim,:) = In;

%采用PadArray对前后一行数据进行复制；
if(edgeSizeM ~= 0 || edgeSizeN ~= 0)
InExtend(1:halfedgesizem,halfedgesizen+1:halfedgesizen+nim,:) = padarray(In(1,:,:),[halfedgesizem-1;0],'circular','pre');
InExtend(mextend:-1:mextend-halfedgesizem+1,halfedgesizen+1:halfedgesizen+nim,:) = padarray(In(mim,:,:),[halfedgesizem-1;0],'circular','post');

%采用padarray对前后一列数据进行复制；
InExtend(:,1:halfedgesizen,:) = padarray(InExtend(:,halfedgesizen+1,:),[0,halfedgesizen-1],'circular','pre');
InExtend(:,nextend:-1:nextend-halfedgesizen+1,:) = padarray(InExtend(:,nextend-halfedgesizen-1,:),[0,halfedgesizen-1],'circular','post');
end

