function [D,M1,M2,M3,W] = pre_calculations_exp(a,b,K,d)

D =cell(1,d);
M1 = cell(1,d);   % q_k^2
M2 = cell(1,d-1); % q_k+1^3
M3 = cell(1,d-1); % q_k+1
W = cell(1,d);

% Initializaitions
eta = -1;
qi_l = -6;
qi_r = 6;
br = 3;
bl = 3;


n = (-(K(1)-1)/2):1:((K(1)-1)/2);

% discrete Laplacian
% D matrices
D1 = zeros(K(1),K(1));
for i=1:K(1) % -0.5 as we have -0.5\Delta
    D1(i,i) = -0.5*((2*1i*pi*n(i))/(b-a))^2;
end

% potential
M_1 = zeros(K(1),K(1));
M_2 = zeros(K(1),K(1));
M_3 = zeros(K(1),K(1));
W1 = zeros(K(1),K(1));

for i=1:K(1)
    for j=1:K(1)
        % q_k^2
        f = prod1_phi(n(i),n(j),a,b);
        Int1 = integral(f,a,b);
        
        % q_k+1
        f = prod3_phi(n(i),n(j),a,b);
        Int2 = integral(f,a,b);
        
        % q_k+1^3
        f = prod2_phi(n(i),n(j),a,b);
        Int3 = integral(f,a,b);
        
        % CAP potential
        f1 = CAP_1(n(i),n(j),a,b,qi_r,br);
        Int_C1 = integral(f1,a,b);
%         Int_C1 = integral(f1,qi_r,b);
        f2 = CAP_2(n(i),n(j),a,b,qi_l,bl);
        Int_C2 = integral(f2,a,b);
%         Int_C2 = integral(f2,a,qi_l);
%         W1(i,j) = 1i*eta*(Int_C1 + Int_C2);
        W1(i,j) = -eta*(Int_C1 + Int_C2);
        
        % create M_i matrices
        M_1(i,j) = Int1;
        M_2(i,j) = Int3;
        M_3(i,j) = Int2;
        
    end
end

for k=1:d
    D{k} = D1;
%     D{k} = zeros(size(D1));
    M1{k} = M_1;
% 	M1{k} = zeros(size(D1));
    if k<d
        M2{k} = M_2;
%         M2{k} = zeros(size(D1));
        M3{k} = M_3;
%         M3{k} = zeros(size(D1));
    end
    W{k} = W1;
%     W{k} = zeros(size(W1));
%     W{k} = -0.1*1i*ones(size(W1));
end

end


function[f] = prod1_phi(i,j,a,b)
% defines integrand of q_k^2 * phi_i*phi_j

f = @(x) x.^2 .*(1/(b-a)).*exp(-(2*1i*i*pi*(x-a))/(b-a)).*...
    exp((2*1i*j*pi*(x-a))/(b-a));
end

function[f] = prod2_phi(i,j,a,b)
% defines integrand of q_k^3 * phi_i*phi_j

f = @(x) x.^3 .*(1/(b-a)).*exp(-(2*1i*i*pi*(x-a))/(b-a)).*...
    exp((2*1i*j*pi*(x-a))/(b-a));
end

function[f] = prod3_phi(i,j,a,b)
% defines integrand of  q_k * phi_i*phi_j

f = @(x) x .*(1/(b-a)).*exp(-(2*1i*i*pi*(x-a))/(b-a)).*...
    exp((2*1i*j*pi*(x-a))/(b-a));
end

function[f] = CAP_1(i,j,a,b,qi_r,br)
% defines integrand of (q_i - qi_r)^3

f = @(x) max(0,(x-qi_r)).^br ...
    .*(1/(b-a)).*exp(-(2*1i*i*pi*(x-a))/(b-a)).*...
    exp((2*1i*j*pi*(x-a))/(b-a));
end

function[f] = CAP_2(i,j,a,b,qi_l,bl)
% defines integrand of (q_i - qi_l)^3

f = @(x) min(0,(x-qi_l)).^bl ...
    .*(1/(b-a)).*exp(-(2*1i*i*pi*(x-a))/(b-a)).*...
    exp((2*1i*j*pi*(x-a))/(b-a));
end

% function[f] = CAP_1(i,j,a,b,qi_r,br)
% % defines integrand of (q_i - qi_r)^3
% 
% f = @(x) (x-qi_r).^br ...
%     *(1/(b-a)).*exp(-(2*1i*i*pi*(x-a))/(b-a)).*...
%     exp((2*1i*j*pi*(x-a))/(b-a));
% end
% 
% function[f] = CAP_2(i,j,a,b,qi_l,bl)
% % defines integrand of (q_i - qi_l)^3
% 
% f = @(x) (x-qi_l).^bl ...
%     *(1/(b-a)).*exp(-(2*1i*i*pi*(x-a))/(b-a)).*...
%     exp((2*1i*j*pi*(x-a))/(b-a));
% end


% for k=1:d
%     n = (-(K(1)-1)/2):1:((K(1)-1)/2);
%
%     % discrete Laplacian
%     % D matrices
%     D1 = zeros(K(1),K(1));
%     for i=1:K(1) % -0.5 as we have -0.5\Delta
%         D1(i,i) = -0.5*((2*1i*pi*n(i))/(b-a))^2;
%     end
%     D{k} = D1;
%
%     % potential
%     M_1 = zeros(K(1),K(1));
%     M_2 = zeros(K(1),K(1));
%     M_3 = zeros(K(1),K(1));
%     W1 = zeros(K(1),K(1));
%
%     for i=1:K(1)
%         for j=1:K(1)
%             % q_k^2
%             f = prod1_phi(n(i),n(j),a,b);
%             Int1 = integral(f,a,b);
%
%             if k<d
%                 % q_k+1
%                 f = prod3_phi(n(i),n(j),a,b);
%                 Int2 = integral(f,a,b);
%
%                 % q_k+1^3
%                 f = prod2_phi(n(i),n(j),a,b);
%                 Int3 = integral(f,a,b);
%             end
%
%             % CAP potential
%             f1 = CAP_1(n(i),n(j),a,b,qi_r,br);
%             Int_C1 = integral(f1,qi_r,b);
%             f2 = CAP_2(n(i),n(j),a,b,qi_l,bl);
%             Int_C2 = integral(f2,a,qi_l);
%             W1(i,j) = 1i*eta*(Int_C1 + Int_C2);
%
%             % create M_i matrices
%             M_1(i,j) = Int1;
%             if k<d
%                 M_2(i,j) = Int3;
%                 M_3(i,j) = Int2;
%             end
%
%         end
%     end
%     M1{k} = M_1;
%     if k<d
%         M2{k} = M_2;
%         M3{k} = M_3;
%     end
%     W{k} = W1;
% end
