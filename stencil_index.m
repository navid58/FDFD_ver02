function r = stencil_index(nz,nx,m,n)

% STENCIL_INDEX gives single index of 9-point stencil nodes, 
% given the subscript of the center node (m,n) of the stencil.
%
% INPUTS
% ======
% nz     : number of grids in z direction
% nx     : number of grids in x direction
% m     : vertical index of center node of the stencil
% n     : horizontal index of center node of the stencil
%
% OUTPUT
% ======
% r : single index value of nodes inside stencil

% By: Navid Amini
% email: amini_navid@yahoo.com

r(1) = nz*(n-2)+(m-1); 
r(2) = nz*(n-2)+(m);   
r(3) = nz*(n-2)+(m+1);

r(4) = nz*(n-1)+(m-1);
r(5) = nz*(n-1)+(m);
r(6) = nz*(n-1)+(m+1); 

r(7) = nz*(n)+(m-1); 
r(8) = nz*(n)+(m); 
r(9) = nz*(n)+(m+1); 

r(r<=0) = 0;
r(r>nz*nx)= 0;

