function obj = updateExpection(obj)

%做出来是信号的期望与卷积的期望。
obj.im_x1 = obj.ImD.opt_cx1 + obj.ba_sigma * obj.im_x1;
obj.im_x2 = obj.ImD.opt_cx2 + obj.ba_sigma * obj.im_x2;

obj.kernel_x1 = obj.KernelD.opt_cx1 + obj.ba_sigma *obj.kernel_x1;
obj.kernel_x2 = obj.KernelD.opt_cx2 + obj.ba_sigma * obj.kernel_x2;

obj.ensemble_x1 = obj.Ensemble.opt_cx1 + obj.ba_sigma * obj.ensemble_x1;
obj.ensemble_x2 = obj.Ensemble.opt_cx2 + obj.ba_sigma * obj.ensemble_x2;


