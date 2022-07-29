function [pf_p,w] = fdfd(vp,rho,q,atten_opt,wref,dx,tmax,twrap,...
                         sname,f0,fmax,L,alpha,top_bc,Rx,Rz,Sx,Sz,use_parfor)

% FDFD solves Ap = s, using direct solver ("\" operator) for each frequency
% with for loop. A is the impedance matrix and s is the source term and p is
% pressure wavefield. Frequency components are independent and "parfor" can
% be used to parallelize over frequencies.
%
% INPUTS
% ======
% vp : Vp velocity model (nz*nx)
% rho : density model (nz*nx)
% q : quality factor model (nz*nx)
% atten_opt: option for attenuation 
% wref : reference angular frequency in Kolsky-Futterman attenuation model
% dx : grid-point interval
% tmax : maximum time of modeling
% twrap: damping to suppress time aliasing
% sname : type of source wavelet
% f0 : dominant frequency of the source wavelet
% fmax : maximum frequency
% L : width of PML layer (number of grid-points)
% alpha : amplitude of PML damping cosine function
% top_bc: boundary condition at top boundary
% Rx : horizontal index of receivers
% Rz : vertical index of receivers
% Sx : horizontal index of sources
% Sz : vertical index of sources
% use_parfor : if 0 sequential "for" loop over frequencies and 
% if 1 parallelizes over frequencies using "parfor" loop. Use “parpool” 
% command, to start MATLAB workers. For more information see MATLAB 
% parallel processing toolbox documents.
%
% OUTPUT
% ======
% pf_r : pressure wavefield in frequency domain (f > 0, positive frequencies)
%        DIM1 size: number of receivers
%        DIM2 size: number of frequencies
%        DIM3 size: number of sources
% w : angular frequency vector
%
% By: Navid Amini
% email: amini_navid@yahoo.com

df = 1/tmax;
f = df:df:fmax; % list of frequencies 
nf = length(f);
fs = source(f,f0,sname); % source wavelet

[nz,nx] = size(vp);
rec_ind = sub2ind([nz,nx],Rz,Rx);
ns = length(Sx);    % number of sources
nr = length(Rx(:)); % number of receivers

% Build sparse RHS term
s = sparse(nz*nx,ns);
for i = 1:ns
    s(sub2ind([nz,nx],Sz(i),Sx(i)),i) = 1;
end

w = 2*pi*f + 1i*twrap/tmax;

pf_p = zeros(nr,nf,ns);
% loop over frequencies
if use_parfor == 0
    
    fprintf ('Starting sequential loop over frequencies ... \n \n')
    
    for k = (1:length(w))
        blkm = bulk_modulus(vp,rho,q,real(w(k)),wref,atten_opt); % bulk modulus
        A = imp_nine(w(k),dx,L,alpha,vp,rho,blkm,top_bc); % impedance matrix
        RHS = s*fs(k);
        temp = full(A\RHS); % direct solver
        pf_p(:,k,:) = temp(rec_ind(:),:);      % extract wavefield at receiver locations
        fprintf ('f%1d  = %0.2f (Hz) solved. \n', k, f(k))
    end
    
elseif use_parfor == 1
   
    fprintf ('Starting parallel loop over frequencies ... \n \n')
    
    parfor k = (1:length(w))
        blkm = bulk_modulus(vp,rho,q,real(w(k)),wref,atten_opt); % bulk modulus
        A = imp_nine(w(k),dx,L,alpha,vp,rho,blkm,top_bc); % impedance matrix
        RHS = s*fs(k);
        temp = full(A\RHS); % direct solver
        pf_p(:,k,:) = temp(rec_ind(:),:);      % extract wavefield at receiver locations
        fprintf ('f%1d  = %0.2f (Hz) solved. \n', k, f(k))
    end
    
end

fprintf ('\nFDFD for all frequencies is done! \n \n')

