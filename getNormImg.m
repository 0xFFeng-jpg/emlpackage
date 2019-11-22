function tmp = getNormImg(in)
tmp = (in - min(in(:)))./(max(in(:)) - min(in(:)));