function psfshape = getPSFfromImage
%从图像中获取PSF形状；
IsGet = 0;
if(IsGet)
    [filename, pathname] = uigetfile({'*.bmp','BMP文件(*.bmp)';'*.jpg', 'JPEG文件(*.jpg)';'*.png','PNG文件(*.png)'});
    if(filename == 0), return, end
    filename = [pathname filename];
    psfImage = imread(filename);
    figure;imshow(psfImage);
    [x y] = ginput(2);
    direction = ((y(1)-y(2))/(x(1) - x(2))>0);
    shape = ones(abs(y(1)-y(2)),abs(x(1) - x(2)));
    shapediag = rectdiag(shape,direction);
    [m,n] = size(shapediag);
    outshape = zeros(m+10,n+10);
    m_start = 5;
    n_start = 5;
    outshape(m_start+1:m_start+m,n_start+1:n_start+n) = shapediag;
else
    outshape = zeros(35,35);
    outshape(18,18) = 1;
end
figure;imshow(outshape,[]);
psfshape = outshape;