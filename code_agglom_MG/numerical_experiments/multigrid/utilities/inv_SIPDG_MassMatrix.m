function M_inv = inv_SIPDG_MassMatrix(M, nelem, nln)


%THE MASS MATRIX IS BLOCK DIAGONAL, 
%      ___           ___
%     | M1 0  . . .  0  |
%     | 0  M2        .  |
%     | 0  0  .      .  |
% M = | .  .    .    .  |
%     | .  .      .     |
%     | .  .    .       |
%     | 0  .  .  .   Mn | 
%      ---           ---  
%
%THE INVERSE IS:
%      ___                           ___
%     | inv(M1)   0    .   .   .    0   |
%     | 0      inv(M2)              .   |
%     | 0         0    .            .   |
% M = | .         .       .         .   |
%     | .         .          .          |
%     | .         .             .       |
%     | 0         .    .  .     inv(Mn) | 
%      ---                           ---  
%
% where n is the number of elements, and
% Mi \in R^{nln x nln}, nln is the number of
% dofs inside each element

N = nelem * nln;

M_inv = sparse(N,N);

for n = 1:nelem
   
    ind_of_block = [ (n-1)*nln + 1  :   1   :   n*nln ];
    
    M_aux = M(ind_of_block, ind_of_block);
    
    M_inv(ind_of_block, ind_of_block) = inv(M_aux);

end



end
