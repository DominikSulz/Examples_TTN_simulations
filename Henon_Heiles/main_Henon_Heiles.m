% This script computes an approximation to the solution of the Schrödinger
% equation, with H = \Delta + Henon-Heiles potential. A detailed description 
%can be found in 
% Time integration of tensor trains, SIAM J. Numer. Anal. 53 (2015), 917-941
clear; clc; close all;

% insert path of the folders for BUG and rank-adaptive BUG, for example
addpath('...\unconventional_TTN_integrator')

% my working directory
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\unconventional_TTN_integrator')

% Initialisations
a = -10;
b = 10;
lamda = 0.111803;

d = 4; % DOF, i.e. number of particles
K = 31*ones(1,d); % number of basis functions in the Fourier base
r = 8*ones(1,d); % initial rank of the system
r_tau = 8;
q = 2*ones(1,d); % positions of each 1-d Gaussain

tol = 10^-8;

t0 = 0;
T_end = 10; 
h = 0.01;
r_max = 10;
r_min = 2;

% initial data
[Y0,tau] = init_Gaussian_binary_tree_all_dim(K,r,d,a,b,q,r_tau);

%% operator construction
t_steps = ceil((T_end - t0)/h);

[D,M1,M2,M3,W] = pre_calculations_exp(a,b,K,d);

% represent the Hamiltonian as a linear operator
B = linearisation(D,M1,M2,M3,W,d,lamda);
A = make_operator(Y0,B,tau,K.*ones(1,d));
A{end} = -1i*A{end};    

ac = zeros(t_steps+1,1);
Norm_con = zeros(t_steps+1,1);
Energy_con = zeros(t_steps+1,1);

ac(1) = Mat0Mat0_prod(Y0,Y0,a,b,d);
Norm_con(1) = sqrt(ac(1));
F = F_HH(t0,Y0,A,d);
Energy_con(1) = 1i*Mat0Mat0_prod(Y0,F,a,b,d); % 1i* as in F_HH I have 1/i *H


Y_start = Y0;
for i=1:t_steps    
%     BUG integrator
    Y_new = TTN_integrator_complex_nonglobal(tau,Y_start,@F_HH,(i-1)*h,i*h,A,d); 

    Norm_con(i+1) = sqrt(Mat0Mat0_prod(Y_new,Y_new,a,b,d));
    F = F_HH(t0+i*h,Y_new,A,d);
    Energy_con(i+1) = 1i*Mat0Mat0_prod(Y_new,F,a,b,d); % 1i* as in F_Ham I have 1/i *H
    
    ac(i+1) = Mat0Mat0_prod(Y_new,Y0,a,b,d);
    Y_start = Y_new;
   
    fprintf('t=%f, E=%f, norm=%f \n', i*h,  abs(Energy_con(i,1)), abs(Norm_con(i,1)))
end


%% compute vibrational spectrum
omega = 0:0.001:30; 

% compute the integral of the FT via the trapezoidal rule
FT3 = zeros(length(omega),1);
time = 0:h:T_end;
for i=1:length(omega)
    for k=1:length(ac)
        if k==1 || k==length(ac)
            FT3(i) = FT3(i) + 0.5*h*ac(k)*exp(-1i*omega(i)*time(k));
        else
            FT3(i) = FT3(i) + h*ac(k)*exp(-1i*omega(i)*time(k));
        end
    end
end
FT3 = (1/pi)*real(FT3);

% Plot the results
figure(1)
plot(omega,abs(FT3),'Linewidth',2)
legend('$S(\omega )$ BUG','Interpreter','latex','FontSize',14);
title('Vibrational spectrum using BUG integrator')

figure(2)
subplot(1,2,1); 
plot(t0:h:T_end,abs(Norm_con),'Linewidth',2)
legend('$ \vert \vert X(t) \vert \vert$','Interpreter','latex','FontSize',14)
subplot(1,2,2); 
plot(t0:h:T_end,abs(Energy_con),'Linewidth',2)
legend('$ < X(t) \vert H \vert X(t) > $','Interpreter','latex','FontSize',14)
