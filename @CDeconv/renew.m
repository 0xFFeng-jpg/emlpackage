function [outDeconv,outGradList] = renew(obj,inDistriList,inGradList,state)

%
outGradList = inGradList;


outDeconv = obj;
outDeconv.Ensemble = inDistriList(1,1);
outDeconv.ImD = inDistriList(1,2);
outDeconv.KernelD = inDistriList(1,3);

outDeconv = outDeconv.blindDeconv();
outDeconv = outDeconv.calcSigmaDx((state == 3));
outDeconv = outDeconv.updateExpection();

outGradList(1,1).x1 = outDeconv.ensemble_x1 - inDistriList(1,1).x1;
outGradList(1,1).x2 = outDeconv.ensemble_x2 - inDistriList(1,1).x2;


outGradList(1,2).x1 = outDeconv.im_x1 - inDistriList(1,2).x1;
outGradList(1,2).x2 = outDeconv.im_x2 - inDistriList(1,2).x2;

outGradList(1,3).x1 = outDeconv.kernel_x1 - inDistriList(1,3).x1;
outGradList(1,3).x2 = outDeconv.kernel_x2 - inDistriList(1,3).x2;

