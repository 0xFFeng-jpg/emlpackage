function InExtend = EdgePadArray(In,edgeSizeM,edgeSizeN)
%ͨ����ͼ��߽��ȵĶ��壬����չͼ��
[mim,nim,lin] = size(In);

mextend = mim + edgeSizeM;
nextend = nim + edgeSizeN;

halfedgesizem = edgeSizeM/2;
halfedgesizen = edgeSizeN/2;

InExtend = zeros(mextend,nextend,lin);

%InExtend = padarray(In,[halfedgesizem,halfedgesizen],'circular','both');

%����λ��
%����λ��Ϊʵ��ͼ��
InExtend(halfedgesizem+1:halfedgesizem+mim,halfedgesizen+1:halfedgesizen+nim,:) = In;

%����PadArray��ǰ��һ�����ݽ��и��ƣ�
if(edgeSizeM ~= 0 || edgeSizeN ~= 0)
InExtend(1:halfedgesizem,halfedgesizen+1:halfedgesizen+nim,:) = padarray(In(1,:,:),[halfedgesizem-1;0],'circular','pre');
InExtend(mextend:-1:mextend-halfedgesizem+1,halfedgesizen+1:halfedgesizen+nim,:) = padarray(In(mim,:,:),[halfedgesizem-1;0],'circular','post');

%����padarray��ǰ��һ�����ݽ��и��ƣ�
InExtend(:,1:halfedgesizen,:) = padarray(InExtend(:,halfedgesizen+1,:),[0,halfedgesizen-1],'circular','pre');
InExtend(:,nextend:-1:nextend-halfedgesizen+1,:) = padarray(InExtend(:,nextend-halfedgesizen-1,:),[0,halfedgesizen-1],'circular','post');
end

