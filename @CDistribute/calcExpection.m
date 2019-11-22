function obj = calcExpection(obj)
%DistributionType Îª×Ö·û´®ÀàÐÍ
x1 = obj.x1;
x2 = obj.x2;
DistributionType = obj.DistriType;
switch DistributionType
    case 'Gaussian'
        mx = x1./x2;
        mx2 = x1.^2.*x2.^-2+1./x2;
        Hx = -0.5 + 0.5*log(x2);
    case 'Laplacian'
        t=-x1./sqrt(2*x2);
        erf_table=erfcx(t);
        
        mx=(x1./x2+sqrt(2/pi./x2)./erf_table).*(t<=25)+(-1./x1+2*x2.*x1.^-3-10*x2.^2.*x1.^-5).*(t>25);
        mx2=(x1.^2.*x2.^-2+1./x2+2*x1./x2./sqrt(2*pi*x2)./erf_table).*(t<=25)+(2.*x1.^-2-10*x2.*x1.^-4+74*x2.^2.*x1.^-6).*(t>25);
        Hx=(-log(erfc(min(t,25)))+0.5*log(2*x2/pi)-0.5+x1./sqrt(2*pi*x2)./erf_table).*(t<25)+(log(-x1)-1+2*x2.*x1.^-2-15*x2.^2.*x1.^-4/2+148*x2.^3.*x1.^-6/3).*(t>=25);
        
    case 'Rectified Gaussian'
        %Rectified Gaussian
        t=-x1./sqrt(2*x2);
        erf_table=erfcx(t);
        
        mx=(x1./x2+sqrt(2/pi./x2)./erf_table).*(t<=25)+(-1./x1+2*x2.*x1.^-3-10*x2.^2.*x1.^-5).*(t>25);
        mx2=(x1.^2.*x2.^-2+1./x2+2*x1./x2./sqrt(2*pi*x2)./erf_table).*(t<=25)+(2.*x1.^-2-10*x2.*x1.^-4+74*x2.^2.*x1.^-6).*(t>25);
        Hx=(-log(erfc(min(t,25)))+0.5*log(2*x2/pi)-0.5+x1./sqrt(2*pi*x2)./erf_table).*(t<25)+(log(-x1)-1+2*x2.*x1.^-2-15*x2.^2.*x1.^-4/2+148*x2.^3.*x1.^-6/3).*(t>=25)+0.5*log(pi/2);
        
    case 'Discrete'
        %Discrete (1,-1)
        mx=tanh(x1);
        mx2=ones(size(x1));
        Hx=x1.*mx-abs(x1)-log(1+exp(-2*abs(x1)))+log(2);
    otherwise 
        error('No the distribution');
end
%set(obj,'Hx',Hx);
% set.mx(mx);
% set.mx2(mx2);
obj.Hx = Hx;
obj.mx = mx;
obj.mx2 = mx2;