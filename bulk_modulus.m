function blkm = bulk_modulus(vp,rho,q,w,wref,atten_opt)

% BULK_MODULUS computes complex bulk modulus according to attenuation mechanism
%
% INPUTS
% ======
% vp  : velocity model
% rho : density model 
% q   : quality factor model 
% w   : angular frequency 
% wref : reference frequency for Kolsky-Futterman attenuation mechanism
% atten_opt : option for attenuation mechanism
% atten_opt = 'KF' : Kolsky-Futterman frequency-dependent attenuation mechanism
% atten_opt = 'no_atten' :  no attenuation, bulk modulus is real in this case
%
% OUTPUT
% ======
% blkm : complex bulk modulus 
%
% By: Navid Amini
% email: amini_navid@yahoo.com


if strcmp(atten_opt,'KF')
    
    % frequency dependent (Kolsky-Futterman)
    vp = 1./( (1./vp)+(1./(pi*vp.*q))*log(wref./w) + 1i*0.5./(vp.*q) );
    blkm = rho.*vp.*vp;
    
elseif strcmp(atten_opt,'no_atten')
    
    blkm = rho.*vp.*vp;
end

