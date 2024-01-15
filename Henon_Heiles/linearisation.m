function [B] = linearisation(D,M1,M2,M3,W,d,lambda)
% This function saves the matrices, which are needed to represent the
% discrete Hamiltonian, in an array B

B = cell(3,d);

for ii=1:d
%     [n,~] = size(D{ii});
%     Z = zeros(n,n);
    % mode k
    if ii==1
        A = D{ii} + 0.5*M1{ii} + W{ii};
    else
        A = D{ii} + 0.5*M1{ii} + W{ii} - (1/3)*lambda*M2{ii-1};
    end
    B{ii,ii} = A;
end

for ii=1:d-1
    
    B{d+ii,ii} = lambda*M1{ii};
    B{d+ii,ii+1} = M3{ii};
end
    
end
