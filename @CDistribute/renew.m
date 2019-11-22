function [outDistri,outGrad] = renew(obj,inGrad,state,steplen,CompType)

%����״����
%inGrad ���Ż����ݶȷ���state��״̬��steplen�ǲ�����
outDistri = obj; %%??
outDistri.CompType = CompType;
outGrad = inGrad;

outDistri.x1 = obj.x1 + steplen * inGrad.x1;
outDistri.x2 = abs(obj.x2 + steplen * inGrad.x2);

outDistri.b_x = abs(obj.b_x + steplen * inGrad.b_x);
outDistri.ba_x = abs(obj.ba_x + steplen * inGrad.ba_x);
outDistri.pi_x = abs(obj.pi_x + steplen * inGrad.pi_x);


outDistri = outDistri.calcExpection();%�˴���components�޹�

outDistri = outDistri.calcLamda();%�˴�Ҳ��Ӧ��componentswise�ı仯�йأ���

outDistri = outDistri.optimizeParaByPrior();%�˴�ֻ�ǻ�ȡ�Ż�ֵ���Ƿ�ı仹Ҫ������state���ж���



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
