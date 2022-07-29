function bu = rho_stg(rho,m,n)

% RHO_STG interpolates buoyancy (inverse of density) at intermediate 
% positions of 9-point stencil
%
% INPUTS
% ======
% rho : density model (nz*nx)
% m   : vertical index of center node of the stencil
% n   : horizontal index of center node of the stencil
%
% OUTPUT
% ======
% bu : interpolated buoyancy values 
%
% By: Navid Amini
% email: amini_navid@yahoo.com

b = 1./rho;
[nz,nx] = size(rho);

if (m-1)>=1  &&  (n-1)>=1,  buNW = b(m-1,n-1); else buNW = b(m,n); end
if               (n-1)>=1,   buW = b(m,n-1);   else buW  = b(m,n); end
if (m+1)<=nz &&  (n-1)>=1,  buSW = b(m+1,n-1); else buSW = b(m,n); end

if (m-1)>=1,               buN  = b(m-1,n);   else buN = b(m,n); end
if (m+1)<=nz,              buS  = b(m+1,n);   else buS = b(m,n); end

if (m-1)>=1  && (n+1)<=nx,  buNE = b(m-1,n+1); else buNE = b(m,n); end
if              (n+1)<=nx,   buE = b(m,n+1);   else buE  = b(m,n); end
if (m+1)<=nz && (n+1)<=nx,  buSE = b(m+1,n+1); else buSE = b(m,n); end

bu(5)  = b(m,n);

bu(1) = 0.5*(buNW+bu(5));
bu(2) = 0.5*(buW+bu(5));
bu(3) = 0.5*(buSW+bu(5));

bu(4) = 0.5*(buN+bu(5));
bu(6) = 0.5*(buS+bu(5));

bu(7) = 0.5*(buNE+bu(5));
bu(8) = 0.5*(buE+bu(5));
bu(9) = 0.5*(buSE+bu(5));