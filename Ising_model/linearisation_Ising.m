function [B] = linearisation_Ising(sx,sz,d)
% This function saves the matrices, which are needed to represent the
% discrete Hamiltonian, in an array B


B = cell(d,d);
for ii=1:d
    B{ii,ii} = sx;
end
for ii=d+1:(2*d-1)
    B{ii,ii-d} = sz;
    B{ii,ii+1-d} = sz;
end



end
