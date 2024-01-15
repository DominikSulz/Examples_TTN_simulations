function [] = pre_calculations_exp_diff_CAP(a,b,K,d)
global  D M1 M2 M3 W

D =cell(1,d);
M1 = cell(1,d);   % q_k^2
M2 = cell(1,d-1); % q_k+1^3
M3 = cell(1,d-1); % q_k+1
W = cell(1,d);

% Initializaitions
qi_0 = 8;


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
        f1 = CAP1(n(i),n(j),a,b,qi_0);
        Int_C1 = integral(f1,a,b);
        f1 = CAP2(n(i),n(j),a,b,-qi_0);
        Int_C2 = integral(f1,a,b);
        W1(i,j) = -1i*(Int_C1 + Int_C2);
        
        % create M_i matrices
        M_1(i,j) = Int1;
        M_2(i,j) = Int3;
        M_3(i,j) = Int2;
        
    end
end

for k=1:d
    D{k} = D1;
    M1{k} = M_1;
    if k<d
        M2{k} = M_2;
        M3{k} = M_3;
    end
    W{k} = W1;
%     W{k} = zeros(size(W1));
%     W{k} = -0.1*1i*ones(size(W1));
end

end


function[f] = prod1_phi(i,j,a,b)
% defines integrand of q_k^2 * phi_i*phi_j

f = @(x) x.^2 *(1/(b-a)).*exp(-(2*1i*i*pi*(x-a))/(b-a)).*...
    exp((2*1i*j*pi*(x-a))/(b-a));
end

function[f] = prod2_phi(i,j,a,b)
% defines integrand of q_k^3 * phi_i*phi_j

f = @(x) x.^3 *(1/(b-a)).*exp(-(2*1i*i*pi*(x-a))/(b-a)).*...
    exp((2*1i*j*pi*(x-a))/(b-a));
end

function[f] = prod3_phi(i,j,a,b)
% defines integrand of  q_k * phi_i*phi_j

f = @(x) x *(1/(b-a)).*exp(-(2*1i*i*pi*(x-a))/(b-a)).*...
    exp((2*1i*j*pi*(x-a))/(b-a));
end

function[f] = CAP1(i,j,a,b,qi_0)
% defines integrand of (q_i - qi_r)^3

f = @(x) max(0,(x-qi_0).*(abs(x-qi_0)).^3) ...
    *(1/(b-a)).*exp(-(2*1i*i*pi*(x-a))/(b-a)).*...
    exp((2*1i*j*pi*(x-a))/(b-a));
end

function[f] = CAP2(i,j,a,b,qi_0)
% defines integrand of (q_i - qi_r)^3

f = @(x) max(0,-(x-qi_0).*(abs(x-qi_0)).^3) ...
    *(1/(b-a)).*exp(-(2*1i*i*pi*(x-a))/(b-a)).*...
    exp((2*1i*j*pi*(x-a))/(b-a));
end

