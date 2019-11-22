function [OutImDx,OutImDy,OutPSF] = SingleScaleProcess5(OriginImDx,OriginImDy,NewImDx,NewImDy,InPSF)

D = [OriginImDx OriginImDy];
mx = [NewImDx NewImDy];
me = InPSF./sum(InPSF(:));

[mx_est,me_est] = DFLearning(mx,me,D);

OutPSF = me_est;

OutImDx = mx_est(:,1:size(OriginImDx,2));
OutImDy = mx_est(:,size(OriginImDx,2)+1:end);