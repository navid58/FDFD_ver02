function fs = source(f,f0,sname)

% SOURCE function generates a wavelet pulse in the frequency domain.
% this function generates both Ricker and first derivative of Gaussian
% wavelet.
%
% INPUTS
% ======
% f      : frequencies vector
% f0     : dominant frequency of the source
% sname  : name of wavelet
%
% OUTPUT
% ======
% fs    : source wavelet pressure in frequency domain.
%
% By: Navid Amini
% email: amini_navid@yahoo.com

if strcmp(sname,'gaussian')
    
    t0 = -2/(pi*f0);
    fs = 1i*(2*pi*f).*sqrt(pi/f0).*exp(-.5*(f/f0).^2 + 1i*(2*pi*f)*t0);
    
elseif strcmp(sname,'ricker')
    
    fs = sqrt(4/(pi*f0^2))*((f/f0).^2).*exp(-(f/f0).^2).*exp(-1i*2*pi*f/f0);
    
end