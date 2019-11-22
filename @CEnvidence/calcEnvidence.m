function [oDitriList,oGradList,oDeconvObj,oEnvidence] = calcEnvidence(obj,DitriList,GradList,DeconvObj,state,steplen,CompType)

%
oEnvidence = 0;



DitriNum = size(DitriList,2);
tmpDitriList = DitriList;
tmpGradList = GradList;
for iter = 1:DitriNum
    inDitri = DitriList(1,iter);
    inGrad  = GradList(1,iter);
    [outDitri,outGrad] = inDitri.renew(inGrad,state,steplen,CompType(iter));
    tmpDitriList(1,iter) = outDitri;
    tmpGradList(1,iter) = outGrad;
end


[oDeconvObj,oGradList] = DeconvObj.renew(tmpDitriList,tmpGradList,state);
%º∆À„envidence;
for iter = 1:DitriNum
    oEnvidence = oEnvidence + tmpDitriList(1,iter).KLDistance;
end

oDitriList = tmpDitriList;

oEnvidence = oEnvidence + oDeconvObj.D_x;
oEnvidence = oEnvidence./oDeconvObj.data_points*log2(exp(1));



