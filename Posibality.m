function a = Posibality(x,type,par1,par2)
%
switch type
    case 'Gaussian'
        a = sqrt(par2./(2*pi)).*exp(-(x-par1).^2*0.5*par2);
    otherwise
        disp('No this distribution');
end

        