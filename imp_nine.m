function A = imp_nine(w,dx,L,alpha,vp,rho,blkm,top_bc)

% IMP_NINE constructs impedance matrix based on mixed-grid approach in the
% interior domain and hybrid PML+ABC as the absorbing boundary condition.
%
% For more information about the theory of this work please see the following paper:
%
% Amini, N. and Javaherian, A., “A MATLAB-based frequency-domain 
% finite-difference finite-difference package for solving 2D visco-acoustic 
% wave equation”,Waves in Random and Complex Media, vol. 21, no. 1, pp. 161–183, 2011.
% doi:10.1080/17455030.2010.537708.
%
% Please cite the above paper when reporting, reproducing or extending the results.
%
% INPUTS
% ======
% w : angular frequency
% dx : grid-point interval
% L : width of PML layer (number of grid-points)
% alpha : amplitude of PML damping cosine function
% vp :  Vp velocity
% rho : density
% blkm : bolk mudulus
% top_bc: boundary condition at top of the model
% if:
% top_bc = 'PML' : PML absorbing layer at top
% top_bc = 'Dirichlet' : Dirichlet boundary condition at top (fixed boundary)
% top_bc = 'Neumann' : Neumann boundary condition at top (free boundary)
%
% OUTPUT
% ======
% A : impedance matrix
%
% By: Navid Amini
% email: amini_navid@yahoo.com

[nz,nx]=size(vp);

% pre-allocate the triplets I, J, and V such that A(I(k),J(k)) = V(k)
I = zeros(9*nz*nx,1);
J = zeros(9*nz*nx,1);
V = zeros(9*nz*nx,1);

% PML damping templates
ind = 1:L;
damp  = alpha*(1-cos((L-ind+1)*pi/2/(L)));
damp_z = [damp,zeros(1,nx-2*L),damp(end:-1:1)]; % vertical damping template

if strcmp(top_bc,'PML')
    damp_x = [damp,zeros(1,nz-2*L),damp(end:-1:1)];       % horizontal damping template
else
    damp_x = [zeros(1,L),zeros(1,nz-2*L),damp(end:-1:1)]; % horizontal damping template
end

k = 1; % non-zero entries counter

%%%% Top left corner (1,1)
m = 1; n = 1;
r = stencil_index(nz,nx,m,n);
tmp = -w/vp(m,n);

I(k) = r(5); J(k) = r(5); V(k) = - 1/dx - 1i*tmp*sqrt(2)/4; k = k+1;
I(k) = r(5); J(k) = r(9); V(k) = + 1/dx - 1i*tmp*sqrt(2)/4; k = k+1;
I(k) = r(5); J(k) = r(6); V(k) =  - 1i*tmp*sqrt(2)/4; k = k+1;
I(k) = r(5); J(k) = r(8); V(k) =  - 1i*tmp*sqrt(2)/4; k = k+1;

%%%% Left boundary (n=1)
n = 1;
for m = 2:nz-1
    r = stencil_index(nz,nx,m,n);
    tmp = -w/vp(m,n);
    
    I(k) = r(5); J(k) = r(5); V(k) = - 2*1i*tmp/dx +1i/(tmp*dx^3) + tmp*tmp -1.5/dx^2; k = k+1;
    I(k) = r(5); J(k) = r(8); V(k) = + 2*1i*tmp/dx -1i/(tmp*dx^3) + tmp*tmp -1.5/dx^2; k = k+1;
    I(k) = r(5); J(k) = r(4); V(k) = - 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
    I(k) = r(5); J(k) = r(6); V(k) = - 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
    I(k) = r(5); J(k) = r(9); V(k) = + 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
    I(k) = r(5); J(k) = r(7); V(k) = + 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
end

%%% PML Left boundary block (n=2:L)
for n = 2:L
    for m = 2:nz-1
        Ez = 1 + 1i*damp_z(n)/w;
        Ez_L = 0.5*(2 + 1i*(damp_z(n)+damp_z(n-1))/w);
        Ez_R = 0.5*(2 + 1i*(damp_z(n)+damp_z(n+1))/w);
        
        Ex = 1 + 1i*damp_x(m)/w;
        Ex_T = 0.5*(2 + 1i*(damp_x(m)+damp_x(m-1))/w);
        Ex_B = 0.5*(2 + 1i*(damp_x(m)+damp_x(m+1))/w);
        
        r = stencil_index(nz,nx,m,n);
        
        I(k) = r(5); J(k) = r(2); V(k) = + 1/(Ez*Ez_L*rho(m,n)*dx^2); k = k+1;
        I(k) = r(5); J(k) = r(4); V(k) = + 1/(Ex*Ex_T*rho(m,n)*dx^2); k = k+1;
        I(k) = r(5); J(k) = r(5); V(k) = + w*w/blkm(m,n) -1/(Ex*Ex_T*rho(m,n)*dx^2) ...
            -1/(Ex*Ex_B*rho(m,n)*dx^2) ...
            -1/(Ez*Ez_L*rho(m,n)*dx^2) ...
            -1/(Ez*Ez_R*rho(m,n)*dx^2);   k = k+1;
        I(k) = r(5); J(k) = r(6); V(k) = + 1/(Ex*Ex_B*rho(m,n)*dx^2); k = k+1;
        I(k) = r(5); J(k) = r(8); V(k) = + 1/(Ez*Ez_R*rho(m,n)*dx^2); k = k+1;
    end
end

%%%% Bottom left corner (n=1,m=nz)
m = nz; n = 1;
r = stencil_index(nz,nx,m,n);
tmp = -w/vp(m,n);

I(k) = r(5); J(k) = r(5); V(k) = + 1/dx + 1i*tmp*sqrt(2)/4; k = k+1;
I(k) = r(5); J(k) = r(7); V(k) = - 1/dx + 1i*tmp*sqrt(2)/4; k = k+1;
I(k) = r(5); J(k) = r(4); V(k) =  + 1i*tmp*sqrt(2)/4; k = k+1;
I(k) = r(5); J(k) = r(8); V(k) =  + 1i*tmp*sqrt(2)/4; k = k+1;


%%% Top boundary (m=1)
if strcmp(top_bc,'PML')
    
    m = 1;
    for n = 2:nx-1
        r = stencil_index(nz,nx,m,n);
        tmp = -w/vp(m,n);
        
        I(k) = r(5); J(k) = r(5); V(k) = - 2*1i*tmp/dx + 1i/(tmp*dx^3) + tmp*tmp - 1.5/dx^2; k = k+1;
        I(k) = r(5); J(k) = r(6); V(k) = + 2*1i*tmp/dx - 1i/(tmp*dx^3) + tmp*tmp - 1.5/dx^2; k = k+1;
        I(k) = r(5); J(k) = r(2); V(k) = - 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
        I(k) = r(5); J(k) = r(8); V(k) = - 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
        I(k) = r(5); J(k) = r(3); V(k) = + 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
        I(k) = r(5); J(k) = r(9); V(k) = + 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
    end
    
    %%% PML Top boundary block (m=2:L)
    for n = 2:nx-1
        for m = 2:L
            Ez = 1 + 1i*damp_z(n)/w;
            Ez_L = 0.5*(2 + 1i*(damp_z(n)+damp_z(n-1))/w);
            Ez_R = 0.5*(2 + 1i*(damp_z(n)+damp_z(n+1))/w);
            
            Ex = 1 + 1i*damp_x(m)/w;
            Ex_T = 0.5*(2 + 1i*(damp_x(m)+damp_x(m-1))/w);
            Ex_B = 0.5*(2 + 1i*(damp_x(m)+damp_x(m+1))/w);
            
            r = stencil_index(nz,nx,m,n);
            
            I(k) = r(5); J(k) = r(2); V(k) = + 1/(Ez*Ez_L*rho(m,n)*dx^2); k = k+1;
            I(k) = r(5); J(k) = r(4); V(k) = + 1/(Ex*Ex_T*rho(m,n)*dx^2); k = k+1;
            I(k) = r(5); J(k) = r(5); V(k) = + w*w/blkm(m,n) -1/(Ex*Ex_T*rho(m,n)*dx^2) ...
                -1/(Ex*Ex_B*rho(m,n)*dx^2) ...
                -1/(Ez*Ez_L*rho(m,n)*dx^2) ...
                -1/(Ez*Ez_R*rho(m,n)*dx^2);   k = k+1;
            I(k) = r(5); J(k) = r(6); V(k) = + 1/(Ex*Ex_B*rho(m,n)*dx^2); k = k+1;
            I(k) = r(5); J(k) = r(8); V(k) = + 1/(Ez*Ez_R*rho(m,n)*dx^2); k = k+1;
        end
    end
    
elseif strcmp(top_bc,'Dirichlet')
    m = 1;
    for n = 2:nx-1
        r = stencil_index(nz,nx,m,n);
        tmp = -w/vp(m,n);
        
        I(k) = r(5); J(k) = r(5); V(k) = -2/dx^2 + tmp^2; k = k+1;
        I(k) = r(5); J(k) = r(6); V(k) = 0; k = k+1;
    end
    
elseif strcmp(top_bc,'Neumann')
    m = 1;
    for n = 2:nx-1
        r = stencil_index(nz,nx,m,n);
        tmp = -w/vp(m,n);
        
        I(k) = r(5); J(k) = r(5); V(k) = -2/dx^2 + tmp^2; k = k+1;
        I(k) = r(5); J(k) = r(6); V(k) =  2/dx^2 ; k = k+1;
    end
    
end

%%%% Intorior
% Jo el al. (1996) 9-point scheme weights
a = 0.5461;
b = 0.6248;
c = 0.25*(1-b);
d = 0;

if strcmp(top_bc,'PML')
    LL  = L;
else
    LL = 1;
end

for n = L+1:nx-L
    for m = LL+1:nz-L
        
        r = stencil_index(nz,nx,m,n);
        bu = rho_stg(rho,m,n);
        
        I(k) = r(5); J(k) = r(1); V(k) = + (1-a)*bu(1)/(2*dx^2);                     k = k+1;
        I(k) = r(5); J(k) = r(2); V(k) = + (w^2)*c/blkm(m,n-1)+(a*bu(2)/(dx^2));     k = k+1;
        I(k) = r(5); J(k) = r(3); V(k) = + (1-a)*bu(3)/(2*dx^2);                     k = k+1;
        I(k) = r(5); J(k) = r(4); V(k) = + (w^2)*c/blkm(m-1,n) + (a*bu(4)/(dx^2));   k = k+1;
        I(k) = r(5); J(k) = r(5); V(k) = +  w*w*(b/blkm(m,n)) - a*(bu(2)+bu(4)+bu(6)+bu(8))/dx^2 ...
                                         -  0.5*(1-a)*(bu(1)+bu(3)+bu(7)+bu(9))/dx^2; k = k+1;
        I(k) = r(5); J(k) = r(6); V(k) = + (w^2)*c/blkm(m+1,n) + (a*bu(6)/(dx^2));    k = k+1;
        I(k) = r(5); J(k) = r(7); V(k) = + (1-a)*bu(7)/(2*dx^2);                      k = k+1;
        I(k) = r(5); J(k) = r(8); V(k) = + (w^2)*c/blkm(m,n+1) +(a*bu(8)/(dx^2));     k = k+1;
        I(k) = r(5); J(k) = r(9); V(k) = + (1-a)*bu(9)/(2*dx^2);                      k = k+1;
    end
end

%%%% Bottom boundary (m=nz)
m = nz;
for n = 2:nx-1
    r = stencil_index(nz,nx,m,n);
    tmp = w/vp(m,n);
    
    I(k) = r(5); J(k) = r(5); V(k) = + 2*1i*tmp/dx - 1i/(tmp*dx^3) + tmp*tmp - 1.5/dx^2;   k = k+1;
    I(k) = r(5); J(k) = r(4); V(k) = - 2*1i*tmp/dx + 1i/(tmp*dx^3) + tmp*tmp - 1.5/dx^2;   k = k+1;
    I(k) = r(5); J(k) = r(2); V(k) = + 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
    I(k) = r(5); J(k) = r(8); V(k) = + 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
    I(k) = r(5); J(k) = r(1); V(k) = - 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
    I(k) = r(5); J(k) = r(7); V(k) = - 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
end

%%% PML Bottom boundary block (m=2:nz-1)
for n = 2:nx-1
    for m = nz-L+1:nz-1
        Ez = 1 + 1i*damp_z(n)/w;
        Ez_L = 0.5*(2 + 1i*(damp_z(n)+damp_z(n-1))/w);
        Ez_R = 0.5*(2 + 1i*(damp_z(n)+damp_z(n+1))/w);
        
        Ex = 1 + 1i*damp_x(m)/w;
        Ex_T = 0.5*(2 + 1i*(damp_x(m)+damp_x(m-1))/w);
        Ex_B = 0.5*(2 + 1i*(damp_x(m)+damp_x(m+1))/w);
        
        r = stencil_index(nz,nx,m,n);
        
        I(k) = r(5); J(k) = r(2); V(k) = + 1/(Ez*Ez_L*rho(m,n)*dx^2); k = k+1;
        I(k) = r(5); J(k) = r(4); V(k) = + 1/(Ex*Ex_T*rho(m,n)*dx^2); k = k+1;
        I(k) = r(5); J(k) = r(5); V(k) = + w*w/blkm(m,n) -1/(Ex*Ex_T*rho(m,n)*dx^2) ...
            -1/(Ex*Ex_B*rho(m,n)*dx^2) ...
            -1/(Ez*Ez_L*rho(m,n)*dx^2) ...
            -1/(Ez*Ez_R*rho(m,n)*dx^2);   k = k+1;
        I(k) = r(5); J(k) = r(6); V(k) = + 1/(Ex*Ex_B*rho(m,n)*dx^2); k = k+1;
        I(k) = r(5); J(k) = r(8); V(k) = + 1/(Ez*Ez_R*rho(m,n)*dx^2); k = k+1;
    end
end

%%%% Top right corner (m=1 , n =nx)
n = nx; m = 1;
r = stencil_index(nz,nx,m,n);
tmp = -w/vp(m,n);

I(k) = r(5); J(k) = r(5); V(k) = - 1/dx -1i*tmp*sqrt(2)/4; k = k+1;
I(k) = r(5); J(k) = r(3); V(k) = + 1/dx -1i*tmp*sqrt(2)/4; k = k+1;
I(k) = r(5); J(k) = r(6); V(k) = - 1i*tmp*sqrt(2)/4; k = k+1;
I(k) = r(5); J(k) = r(2); V(k) = - 1i*tmp*sqrt(2)/4; k = k+1;


%%%% Right boundary (n = nx)
n = nx;
for m = 2:nz-1
    r = stencil_index(nz,nx,m,n);
    tmp = w/vp(m,n);
    
    I(k) = r(5); J(k) = r(5); V(k) = + 2*1i*tmp/dx -1i/(tmp*dx^3) + tmp*tmp -1.5/dx^2; k = k+1;
    I(k) = r(5); J(k) = r(2); V(k) = - 2*1i*tmp/dx +1i/(tmp*dx^3) + tmp*tmp -1.5/dx^2; k = k+1;
    I(k) = r(5); J(k) = r(4); V(k) = + 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
    I(k) = r(5); J(k) = r(6); V(k) = + 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
    I(k) = r(5); J(k) = r(3); V(k) = - 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
    I(k) = r(5); J(k) = r(1); V(k) = - 1i/(2*tmp*dx^3) + 3/(4*dx^2); k = k+1;
end

%%% PML Right boundary block (n = nx-L+1:nx-1 , m = L+1:nz-L)
for n = nx-L+1:nx-1
    for m = 2:nz-1
        Ez = 1 + 1i*damp_z(n)/w;
        Ez_L = 0.5*(2 + 1i*(damp_z(n)+damp_z(n-1))/w);
        Ez_R = 0.5*(2 + 1i*(damp_z(n)+damp_z(n+1))/w);
        
        Ex = 1 + 1i*damp_x(m)/w;
        Ex_T = 0.5*(2 + 1i*(damp_x(m)+damp_x(m-1))/w);
        Ex_B = 0.5*(2 + 1i*(damp_x(m)+damp_x(m+1))/w);
        
        r = stencil_index(nz,nx,m,n);
        
        I(k) = r(5); J(k) = r(2); V(k) = + 1/(Ez*Ez_L*rho(m,n)*dx^2); k = k+1;
        I(k) = r(5); J(k) = r(4); V(k) = + 1/(Ex*Ex_T*rho(m,n)*dx^2); k = k+1;
        I(k) = r(5); J(k) = r(5); V(k) = + w*w/blkm(m,n) -1/(Ex*Ex_T*rho(m,n)*dx^2) ...
            -1/(Ex*Ex_B*rho(m,n)*dx^2) ...
            -1/(Ez*Ez_L*rho(m,n)*dx^2) ...
            -1/(Ez*Ez_R*rho(m,n)*dx^2);   k = k+1;
        I(k) = r(5); J(k) = r(6); V(k) = + 1/(Ex*Ex_B*rho(m,n)*dx^2); k = k+1;
        I(k) = r(5); J(k) = r(8); V(k) = + 1/(Ez*Ez_R*rho(m,n)*dx^2); k = k+1;
    end
end

%%%% Bottom right corner (m = nz , n = nx)
n = nx; m = nz;
r = stencil_index(nz,nx,m,n);
tmp = -w/vp(m,n);

I(k) = r(5); J(k) = r(5); V(k) = + 1/dx + 1i*tmp*sqrt(2)/4; k = k+1;
I(k) = r(5); J(k) = r(1); V(k) = - 1/dx + 1i*tmp*sqrt(2)/4; k = k+1;
I(k) = r(5); J(k) = r(4); V(k) = + 1i*tmp*sqrt(2)/4; k = k+1;
I(k) = r(5); J(k) = r(2); V(k) = + 1i*tmp*sqrt(2)/4;

ind1 = find(I==0);
I(ind1) = [];
J(ind1) = [];
V(ind1) = [];

ind2 = find(J==0);
I(ind2) = [];
J(ind2) = [];
V(ind2) = [];

A = sparse(I,J,V);

A = conj(A);

