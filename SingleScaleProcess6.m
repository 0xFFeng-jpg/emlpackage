function [OutImGrad,OutPSF] = SingleScaleProcess6(OriginImGrad,NewImGrad,InPSF,IsStructure)

D = OriginImGrad;
mx = NewImGrad;
me = InPSF./sum(InPSF(:));

[mx_est,me_est] = DFLearning(mx,me,D,IsStructure);

OutPSF = me_est;

OutImGrad = mx_est;