function [Y,tau] = init_Gaussian_binary_tree_all_dim(K,r,d,a,b,q,r_tau)
% This funciton creates the initial data of a product of Gauissians, 
% d = 2^N in the form of a binary tree

% determine N and c: d = 2^N + c;
N = 0;
c = d - 2^N;
while 2^N <= c
    N = N + 1;
    c = d - 2^N;
end
    
[tau,Y] = create_binary_tree(N,c);

for k=1:d
    n = (-(K(k)-1)/2):1:((K(k)-1)/2);
    F_coef = zeros(K(k),1);
    for i=1:K(k)
        f = @(x) exp(-((x-q(k)).^2) / 2).*...
                 exp(-(2*pi*1i*n(i)*(x-a))/(b-a));
        F_coef(i) = (1/(b-a))*integral(f,a,b);
    end
    
    U = [F_coef (ones(K(k),r(k)-1)+1i*ones(K(k),r(k)-1))];
    Y = set_leaf(Y,U,k,d,r_tau,r);
end
Y{end-1} = 1;
T = zeros(r_tau,r_tau);
T(1,1) = 1;
Y{end} = tensor(T,[r_tau r_tau 1]);

Y = orth_tree(Y);

end

function [tau,Y] = create_binary_tree(k,c)

if k==1 && c==0
    tau = cell(1,2);
    Y = cell(1,4);
    tau{1} = [];
    tau{2} = [];
    Y{1} = [];
    Y{2} = [];
elseif k==1 && c==1
    tau = cell(1,2);
    Y = cell(1,4);
    tau{1} = cell(1,2);
    tau{1}{1} = [];
    tau{1}{2} = [];
    Y{1} = cell(1,4);
    Y{1}{1} = [];
    Y{1}{2} = [];
    tau{2} = [];
    Y{2} = [];
elseif k==1 && c==2
    tau = cell(1,2);
    Y = cell(1,4);
    tau{1} = cell(1,2);
    tau{1}{1} = [];
    tau{1}{2} = [];
    Y{1} = cell(1,4);
    Y{1}{1} = [];
    Y{1}{2} = [];
    tau{2} = tau{1};
    Y{2} = Y{1};
else
    tau = cell(1,2);
    Y = cell(1,4);
    if (c - 2^(k-1)) > 0
        c_right = c - 2^(k-1);
        c_left = 2^(k-1);
    else 
        c_right = 0;
        c_left = c;
    end
    [tau{1},Y{1}] = create_binary_tree(k-1,c_left);
    [tau{2},Y{2}] = create_binary_tree(k-1,c_right);
end

end

function [Z] = set_leaf(Y,S,k,d,r_tau,r)
% S is the matrix that shall be on k-th leaf

Z = Y;
if d==2
    Z{k} = S;
    if mod(k,2) == 0 % even k
        Z{end-1} = eye(r_tau,r_tau);
        M = zeros(r(k-1),r(k),r_tau);
        M(1,1,1) = 1;
        Z{end} = tensor(M);
    else % odd k
        Z{end-1} = eye(r_tau,r_tau);
        M = zeros(r(k),r(k+1),r_tau);
        M(1,1,1) = 1;
        Z{end} = tensor(M);
    end
elseif d==3 
    if k==3
        Z{2} = S;
        Z{end-1} = eye(r_tau,r_tau);
        M = zeros(r(k-1),r(k),r_tau);
        M(1,1,1) = 1;
        Z{end} = tensor(M);
    else
        Z{1} = set_leaf(Y{1},S,k,2,r_tau,r);
    end
else
    N = 0;
    c = d - 2^N;
    while 2^N <= c
        N = N + 1;
        c = d - 2^N;
    end
    
    if c<= (2^(N-1))
        dim_left = 2^(N-1) + c;
        dim_right = 2^(N-1);
    else
        dim_left = 2^N;
        dim_right = c;
    end

    if k <= dim_left
        k_new = k;
        d_new = dim_left;
        Z{1} = set_leaf(Y{1},S,k_new,d_new,r_tau,r);
    else
        k_new = k - dim_left;
        d_new = dim_right;
        Z{2} = set_leaf(Y{2},S,k_new,d_new,r_tau,r);
    end
    M = eye(r_tau,r_tau);
    Z{end-1} = M;
    M = zeros(r_tau,r_tau,r_tau);
    M(1,1,1) = 1;
    Z{end} = tensor(M);
end     
end

function [X] = orth_tree(T)
% This function orthonormalizes the core of T

X = T;
m = length(T) - 2;

for i=1:m
    if 1==iscell(T{i})
        dum2 = orth_tree(T{i});
        s = size(dum2{end});
        ss = length(s);
        Mat = double(tenmat(dum2{end},ss,1:(ss-1)));
        [Q,R] = qr(Mat.',0);
        dum2{end} = tensor(mat2tens(Q.',s,ss,1:(ss-1)));
        X{i} = dum2;
        X{end} = ttm(X{end},R,i);
    else
        [Q,R] = qr(T{i},0);
        X{i} = Q;
        X{end} = ttm(X{end},R,i);
    end
end

end