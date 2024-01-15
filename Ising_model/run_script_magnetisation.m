% This script implements an Ising model in a transverse field with next
% neighbor interactions. A detailed description can be found in 
% Rank-adaptive time integration of tree tensor networks, SIAM J. Numer. Anal. 61 (2023), 194-222

clear all; close all; clc
% insert path of the folders for BUG and rank-adaptive BUG, for example
addpath('...\unconventional_TTN_integrator')
addpath('...\rank_adaptive_integrator_for_TTN')

% my working directory
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\unconventional_TTN_integrator')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\rank_adaptive_integrator_for_TTN')

% number particles
d = 8;
% parameters of the model
Omega = 1;
% time step size and final time
T_end = 1;
dt = 0.025;
% rank of initial data at bottom layer
r = 2;
% bond dimension operator
r_op_max = 4;
r_op_min = 2;

% for rank-adaptive integrator
tol = 10^-8;
r_min = 2;
r_max = 25;

%%% time-step
Iter=T_end/dt;

% initialisations
sx=[0,1;1,0];  %% Pauli Matrix x
sy=[0,-1i;1i,0]; %% Pauli Matrix y
sz=[1,0;0,-1]; %% Pauli Matrix z

% % % create initial data binary tree
[X,tau] = init_spin_all_dim_diff_rank(r,2,d);

% create cell array for Hamiltonian
B = linearisation_Ising(Omega*sx,sz,d);

% make operator of Ising model in TTN representation -> change that to HSS construction
A = make_operator(X,B,tau,2*ones(d,1));
A{end} = -A{end};
A{end} = -1i*A{end};
A = rounding(A,tau);
A = truncate(A,10^-14,r_op_max,r_op_min);

% make operator magnetisation 
B2 = cell(d,d);
for ii=1:d
    B2{ii,ii} = sz;
end
Mag = make_operator(X,B2,tau,2*ones(d,1));

% time evolution
time = [];
obs_sz_ad = [];
obs_sz_BUG = [];

tmp = F_Ising(0,X,A,d);
en_ad = -1i*Mat0Mat0(X,tmp);
nn_ad = sqrt(abs(Mat0Mat0(X,X))); % times i, as we only want the Hamiltonian

en_BUG = -1i*Mat0Mat0(X,tmp);
nn_BUG = sqrt(abs(Mat0Mat0(X,X))); % times i, as we only want the Hamiltonian

r_max_ad = max_rank(X);
r_max_BUG = max_rank(X);

X_start_ad = X;
X_start_BUG = X;

tic
for it=1:Iter
    t0 = (it-1)*dt;
    t1 = it*dt;
    time(it) = t1;
    
    % rank-adaptive
    X_new_ad = TTN_integrator_complex_rank_adapt_nonglobal_spin(tau,X_start_ad,@F_Ising,t0,t1,A,d,r_min);
    X_new_ad = truncate(X_new_ad,tol,r_max,r_min);
    
    % BUG integrator
    X_new_BUG = TTN_integrator_complex_nonglobal_spin(tau,X_start_BUG,@F_Ising,t0,t1,A,d);
    
    % compute energy and norm
    en_ad(it+1) = -1i*Mat0Mat0(X_new_ad,F_Ising(0,X_new_ad,A,d)); % times i, as we only want the Hamiltonian
    nn_ad(it+1) = sqrt(abs(Mat0Mat0(X_new_ad,X_new_ad)));
    
    en_BUG(it+1) = -1i*Mat0Mat0(X_new_BUG,F_Ising(0,X_new_BUG,A,d)); % times i, as we only want the Hamiltonian
    nn_BUG(it+1) = sqrt(abs(Mat0Mat0(X_new_BUG,X_new_BUG)));
    
    % setting for next time step
    X_start_ad = X_new_ad;
    X_start_BUG = X_new_BUG;
    
    % compute magnetisation in z-direction    
    tmp = apply_operator_nonglobal(X_new_ad,Mag,d);
    obs_sz_ad(it) = (1/d)*Mat0Mat0(X_new_ad,tmp);
    
    tmp = apply_operator_nonglobal(X_new_BUG,Mag,d);
    obs_sz_BUG(it) = (1/d)*Mat0Mat0(X_new_BUG,tmp);
    
    % max rank
    r_max_ad(it+1) = max_rank(X_new_ad);
    r_max_BUG(it+1) = max_rank(X_new_BUG);
    
    it
end
total_time = toc;

if d<=8
    [obs_ref,sol] = ref_sol(d,Omega,dt,T_end);
end

figure(1)
plot(time,real(obs_ref),'Linewidth',1.5)
hold on
plot(time,real(obs_sz_ad),'-.','Linewidth',1.5)
plot(time,real(obs_sz_BUG),':','Linewidth',1.5)
xlabel('Time')
title('Magnetization in z-direction')
legend('Exact reference solution','Rank-adaptive BUG','BUG integrator')

figure(2)
plot([0 time],abs(nn_ad),'-.','Linewidth',2)
hold on
plot([0 time],abs(nn_BUG),':','Linewidth',2)
xlabel('Time')
title('Norm')
legend('Rank-adaptive BUG','BUG integrator')

figure(3)
plot([0 time],abs(en_ad),'-.','Linewidth',2)
hold on
plot([0 time],abs(en_BUG),':','Linewidth',2)
xlabel('Time')
title('Energy')
legend('Rank-adaptive BUG','BUG integrator')

figure(4) 
plot([0 time],r_max_ad,'-.','Linewidth',2)
hold on
plot([0 time],r_max_BUG,':','Linewidth',2)
xlabel('Time')
title('Max. ranks')
legend('Rank-adaptive BUG','BUG integrator')

figure(5)
plot(time,abs(obs_ref - obs_sz_BUG),'Linewidth',2)
hold on
plot(time,abs(obs_ref - obs_sz_ad),'Linewidth',2)
legend('Error BUG','Error Rank-adaptive BUG')


function [obs_exact,psi] = ref_sol(d,Omega,dt,T_end)

%%%% dimension
L=d;
h=Omega;

FinalTime=T_end;
Iter=FinalTime/dt;


sx=[0,1;1,0];
sy=[0,-i;i,0];
sz=[1,0;0,-1];

iden=[1 0;0 1];




Proj1=[1 0;0 0];
Proj2=[0 0;0 1];

tmp2=diag(sparse(ones(2^(L-1),1)));
X{1,1}=kron(sx,tmp2);
X{L,1}=kron(tmp2,sx);

Z{1,1}=kron(sz,tmp2);
Z{L,1}=kron(tmp2,sz);


for i1=2:L-1
   tmp1=diag(sparse(ones(2^(i1-1),1)));
   tmp2=diag(sparse(ones(2^(L-i1),1)));
   X{i1,1}=kron(kron(tmp1,sx),tmp2);
   Z{i1,1}=kron(kron(tmp1,sz),tmp2);
end

Mag=0*X{1,1};
H_X=0*X{1,1};

for i1=1:L
    Mag=Mag+Z{i1,1};
    H_X=H_X+X{i1,1};
end

H_NN=0*X{1,1};
for i1=1:L-1
    H_NN=H_NN+Z{i1,1}*Z{i1+1,1};
end

H_sys=-h*H_X-H_NN;

psi=zeros(2^L,1);
psi(1,1)=1;

tic
Vt=expm(-1i*dt*H_sys); % Vt=expm(-1i*dt*H_sys);
toc

for tst=1:Iter
    time(tst)=tst*dt;
    psi=Vt*psi;
    psi=psi/norm(psi);

    p1(tst)=dot(psi,Mag*psi)/L;
    energy(tst)=dot(psi,H_sys*psi)/L;
    
    mag(tst)=dot(psi,Mag*psi)/L;
end


%%%% diagonalization

Ground_State_Energy=min(real(eig(H_sys)))/L;
%%%%% diagonalization of H_eff
[v2,d2]=eig(full(H_sys));
[x2,y2]=sort(diag(d2),'ascend');


 %%% positive eigenvalue
 ind_p=y2(end);


 norm(H_sys*v2(:,ind_p)-d2(ind_p,ind_p)*v2(:,ind_p));
 Ground_state=v2(:,ind_p);
 Ground_state=Ground_state/norm(Ground_state);
 
 
 Mag_Ground_state=dot(Ground_state,Mag*Ground_state)/L;


obs_exact = mag;

end