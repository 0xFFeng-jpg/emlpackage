function [outDistri,outGrad] = renew(obj,inGrad,state,steplen,CompType)

%根据状况；
%inGrad 是优化的梯度方向，state是状态，steplen是步长；
outDistri = obj; %%??
outDistri.CompType = CompType;
outGrad = inGrad;

outDistri.x1 = obj.x1 + steplen * inGrad.x1;
outDistri.x2 = abs(obj.x2 + steplen * inGrad.x2);

outDistri.b_x = abs(obj.b_x + steplen * inGrad.b_x);
outDistri.ba_x = abs(obj.ba_x + steplen * inGrad.ba_x);
outDistri.pi_x = abs(obj.pi_x + steplen * inGrad.pi_x);


outDistri = outDistri.calcExpection();%此处与components无关

outDistri = outDistri.calcLamda();%此处也不应与componentswise的变化有关？？

outDistri = outDistri.optimizeParaByPrior();%此处只是获取优化值；是否改变还要根据其state来判定；



outDistri = outDistri.calcKLDistance();
outDistri = outDistri.calcOptCx();

if(state >= 2)
    outGrad.pi_x = outDistri.opt_pi_x - outDistri.pi_x;
    outGrad.ba_x = outDistri.opt_ba_x - outDistri.ba_x;
    outGrad.b_x = outDistri.opt_b_x - outDistri.b_x;
else
    outGrad.pi_x = 0;
    outGrad.ba_x = 0;
    outGrad.b_x = 0;
end
