function [K_bc, M_bc, K_gl, M_gl]=amatrix(Kele,Mele, node, ndg, H,Ki,Kgb,k, Lam, I)
%amatrix initializes K and M matrices, makes the assembly of the global
%matrices and reduces the global matrices by imposing the boundary
%conditions
 NUM_NODE=length(node);   % the size of nodes
 NUM_ELEM=length(node)-1; % the size of elements
 matrix_size=ndg*NUM_ELEM+ndg;
%
% --------------------------------------------------------------------------
% Initializes K and F matrices
% --------------------------------------------------------------------------
%
 M=zeros(matrix_size,matrix_size);
 K=zeros(matrix_size,matrix_size);
%
% --------------------------------------------------------------------------
% Assembling of K and M matrices
% --------------------------------------------------------------------------
%
% (1). assembly of K matrix
 ELNO=0; % ELNO:the ith element
 for ii=1:ndg:matrix_size-ndg  
    ELNO=ELNO+1;
    [KE]=k_elemu(node(ELNO,2),node(ELNO+1,2),H,Ki,Kgb,k, Kele, ndg);
    K((ELNO*ndg-(ndg-1)):(ELNO+1)*ndg,(ELNO*ndg-(ndg-1)):(ELNO+1)*ndg)=K((ELNO*ndg-(ndg-1)):(ELNO+1)*ndg,(ELNO*ndg-(ndg-1)):(ELNO+1)*ndg)+KE;
    end % end of FOR loop -- assembly of K matrix
%
% (2). assembly of M matrix
 ELNO=0; % ELNO:the ith element
 for ii=1:ndg:matrix_size-ndg  
    ELNO=ELNO+1;
    [ME]=m_elemu(node(ELNO,2),node(ELNO+1,2),Lam, Mele, ndg, I);
     M((ELNO*ndg-(ndg-1)):(ELNO+1)*ndg,(ELNO*ndg-(ndg-1)):(ELNO+1)*ndg)=M((ELNO*ndg-(ndg-1)):(ELNO+1)*ndg,(ELNO*ndg-(ndg-1)):(ELNO+1)*ndg)+ME;
 end
 % end of FOR loop -- assembly of M matrix

% --------------------------------------------------------------------------
% Imposing the essential Boundary Conditions
% --------------------------------------------------------------------------
% In this case, one end of the beam is fixed.
% The other end is free. Therefore, FROM the K and M matrices can be deleted 
% ndg rows and ndg columns. In other words, we can delete the first,  
% second and third rows and columns of the matrices.     
%
% K_bc: apply BCs to the K_matrix
% M_bc: be reduced orders as apply BCs to the K_matrix
%
%   K_bc=zeros(matrix_size-ndg,matrix_size-ndg);
%   M_bc=zeros(matrix_size-ndg,matrix_size-ndg);
%
  K_bc=K((ndg+1):matrix_size,(ndg+1):matrix_size);
  M_bc=M((ndg+1):matrix_size,(ndg+1):matrix_size);
  
  K_gl=K;
  M_gl=M;
end
