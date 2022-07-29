function ve = ext_pml(v,L,top_bc)

% EXT_PML pads grids around the model for applying appropriate boundary condition
%
% INPUTS
% ======
% v : model (velocity or density or quality factor) 
% L : width of PML layer
% top_bc: boundary condition at top of model. if:
% top_bc = 'PML' : pads L grids to all sides of the model
% top_bc = 'Dirichlet' or  top_bc = 'Neumann' : pads L grids to the left, 
% right, and bottom boundaries and one grid to the top boundary
%
% OUTPUT
% ======
% ve : extended model
%
% By: Navid Amini
% email: amini_navid@yahoo.com

if strcmp(top_bc,'PML')
    
    ve = [v(1,1)*ones(L) , repmat(v(1,:),L,1) , v(1,end)*ones(L) ;...
        repmat(v(:,1),1,L) , v , repmat(v(:,end),1,L) ; ...
        v(end,1)*ones(L) , repmat(v(end,:),L,1) , v(end,end)*ones(L)];
    
elseif strcmp(top_bc,'Dirichlet') || strcmp(top_bc,'Neumann')
    
    ve = [repmat(v(:,1),1,L) , v , repmat(v(:,end),1,L) ; ...
        v(end,1)*ones(L) , repmat(v(end,:),L,1) , v(end,end)*ones(L)];
    ve = [ve(1,:) ; ve];
    
end