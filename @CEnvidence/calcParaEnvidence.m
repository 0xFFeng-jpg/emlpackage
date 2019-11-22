function [oDitriList,oGradList,oDeconvObj,oEnvidence] = calcParaEnvidence(obj,DitriList,GradList,DeconvObj,state,steplen,ComType)

oEnvidence = 0;

DitriNum = size(DitriList,2);
tmpDitriList = DitriList;
tmpGradList = GradList;
for iter = 1:DitriNum
    inDitri = DitriList(1,iter);
    inGrad  = GradList(1,iter);
    [outDitri,outGrad] = inDitri.renew(inGrad,state,steplen,ComType(1,iter));
    tmpDitriList(1,iter) = outDitri;
    tmpGradList(1,iter) = outGrad;
end


[oDeconvObj,oGradList] = DeconvObj.renewPara(tmpDitriList,tmpGradList,state);
%º∆À„envidence;
for iter = 1:DitriNum
    oEnvidence = oEnvidence + tmpDitriList(1,iter).KLDistance;
end

oDitriList = tmpDitriList;

oEnvidence = oEnvidence + oDeconvObj.D_x;
oEnvidence = oEnvidence./oDeconvObj.data_points*log2(exp(1));