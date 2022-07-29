function [pt,t] = four2time(pf_p,tmax,twrap,freq_zpad)

% FOUR2TIME converts pressure wavefield from frequency domain to time domain
%  
% INPUTS
% ======
% pf_p : pressure wavefield in frequency domain (f > 0, positive frequencies)
%        DIM1 size: number of receivers
%        DIM2 size: number of frequency components (f > 0)
%        DIM3 size: number of sources
%
% tmax       : maximum time of modeling
% twrap      : damping to suppress time aliasing
% freq_zpad  : number of zeros to pad in the frequency domain to have a smooth waveform in time domain
%
% OUTPUT
% ======
% pt   : pressure wavefield in time domain
%        DIM1 size: number of receivers
%        DIM2 size: number of times samples 
%        DIM3 size: number of sources
% t    : time vector
%
% By: Navid Amini
% email: amini_navid@yahoo.com

[nr,nf,ns] = size(pf_p);

% forming negative, zero and positive frequencies together
pf = cat(2 , conj(pf_p(:,end:-1:1,:)) , zeros(nr,1,ns) , pf_p(:,1:end,:));

% pad freq_zpad zeros to the positive and negative frequencies
pf_pad = cat(2 , zeros(nr,freq_zpad,ns) ,pf ,  zeros(nr,freq_zpad,ns));

% inverse Fourier transfrom
pf = ifftshift(pf_pad,2);
pt = real(ifft(pf,[],2)); % pressure wavefield in the time domain

% time vector
t = linspace(0,tmax, 2*(nf+freq_zpad)+1 );

% undamp wavefieled in time domain
undamp = repmat(exp(twrap*t/tmax),[nr,1,ns]);
pt = pt.*undamp;

fprintf (['Frequency domian to time domain is done! \n \n'])



