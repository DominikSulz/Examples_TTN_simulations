function [Y,tau] = init_spin_all_dim_diff_rank(r,cc,d)
% This funciton creates the initial data of a product state

% determine N and c: d = 2^N + c;
N = 0;
c = d - 2^N;
while 2^N <= c
    N = N + 1;
    c = d - 2^N;
end
    
[tau,Y] = create_binary_tree(N,c);

for k=1:d
    U = [1;0];
    for l=1:r-1
        U = [U [0;1]];
    end
    Y = set_leaf(Y,U,k,d,cc,r);
end
Y{end-1} = 1;
if 1==iscell(Y{1})
    s1 = size(Y{1}{end});
else
    [~,s1] = size(Y{1});
end
if 1==iscell(Y{2})
    s2 = size(Y{2}{end});
else
    [~,s2] = size(Y{2});
end

C = zeros(s1(end),s2(end));
C(1) = 1;
Y{end} = tensor(C,[s1(end) s2(end) 1]);

% adjust the ranks to fulfill rank condition 
s = size(Y{end});
cond = 1;
% check rank condition
for ii=1:length(s)
    tmp = s;
    tmp(ii) = [];
    if prod(tmp) < s(ii)
        cond = 0;
        ss = s(s~=1);
        C = zeros(min(ss),min(ss));
        C(1) = 1;
        Y{end} = tensor(C,[min(ss) min(ss) 1]);
        s = size(Y{end});
    end
end
if cond == 0
    Y = fulfill_rank(Y,s);
end

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

function [Z] = set_leaf(Y,S,k,d,cc,r)
% S is the matrix that shall be on k-th leaf

Z = Y;
if d==2
    Z{k} = S;
    if mod(k,2) == 0 % even k
        if cc==1
            rr = [r r r];
        elseif cc==2
            rr = [r r 2*r]; 
        elseif cc==3
            rr = [r r r^2];
        end
        Z{end-1} = eye(rr(end),rr(end));
        C = zeros(rr);
        C(1) = 1;
        Z{end} = tensor(C,rr);
    else % odd k
        if cc==1
            rr = [r r r];
        elseif cc==2
            rr = [r r 2*r]; 
        elseif cc==3
            rr = [r r r^2];
        end
        Z{end-1} = eye(rr(end),rr(end));
        C = zeros(rr);
        C(1) = 1;
        Z{end} = tensor(C,rr);
    end
elseif d==3 
    if k==3
        Z{2} = S;
        if cc==1
            rr = [r r r];
        elseif cc==2
            rr = [2*r r 2*r]; 
        elseif cc==3
            rr = [r^2 r r^2];
        end
        Z{end-1} = eye(rr(end),rr(end));
        C = zeros(rr);
        C(1) = 1;
        Z{end} = tensor(C,rr);
    else
        Z{1} = set_leaf(Y{1},S,k,2,cc,r);
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
        Z{1} = set_leaf(Y{1},S,k_new,d_new,cc,r);
    else
        k_new = k - dim_left;
        d_new = dim_right;
        Z{2} = set_leaf(Y{2},S,k_new,d_new,cc,r);
    end
    
    m = length(Y)-2;
    for ii=1:m
        if 1 == iscell(Z{ii})
            s = size(Z{ii}{end});
        else
            s = size(Z{ii});
        end
        rr(ii) = s(end);
    end
%     rr(m+1) = 2*rr(1);
    % for a test
    rr(m+1) = 4*rr(2);
    
    if length(rr) == length(find(rr)) % if no Z{ii} is empty
        C = zeros(rr);
        C(1) = 1;
        Z{end-1} = eye(rr(end),rr(end));
        Z{end} = tensor(C,rr);
    end
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

function [Y1] = fulfill_rank(Y,s)
 
    for ii=1:length(s)
        if 1==iscell(Y{ii})
            s2 = size(Y{ii}{end});
            if s2(end) ~= s(ii)
                tmp = s2;
                tmp(end) = s(ii);
                C = zeros(tmp);
                C(1) = 1;
                Y{ii}{end} = tensor(C,tmp);
                Y{ii}{end-1} = eye(s(ii),s(ii));
                s = size(Y{ii}{end});
                Y{ii} = fulfill_rank(Y{ii},s);
            end
        else
            
        end
    end
    Y1 = Y;
    
end