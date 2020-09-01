%-Non-uniform data example in Sect. 3.2------------------------------------
clear all;
close all;
addpath('./src');
addpath('./src/utils');

% Set random seed
rng(1234);

% Interval
a = 0;
b = 2*pi;

% Parameters     
N = 2^10;
M = N;
R = 8;
M_sp = 24;
tau = (1/M^2)*(pi*M_sp)/(R*(R-0.5));

% Non-uniform position vector
x = a + (b/2-a)*rand(N/2,1);
x = [x; b/2+1 + (b-b/2-1)*rand(N/2,1)];
x = sort(x);

% Data vector
amp_1 = 2.0;
amp_2 = 1.0;
f_1 = 50.0;
f_2 = 100.0;
f = 1/N * (amp_1*sin(f_1*x) + amp_2*sin(f_2*x));

% Frequencies
k = (-M/2:M/2-1)';

% NUFFT
F_nufft = nufft1d(f,x,M,R,M_sp,tau);

% Direct summation
F_ds = direct_summation(f,x);

% Relative L2 norm
r_sum_nufft = relative_error_norm(F_nufft,F_ds);
fprintf('Relative L2 norm:  %0.4d\n',r_sum_nufft);

% Single-sided spectrum
F_nufft = abs(F_nufft(M/2:end));
F_ds = abs(F_ds(M/2:end));
freq = k(M/2:end);

% Plot
figure('DefaultAxesFontSize',13);
plot(x,f,'kx-','LineWidth',1);
xlim([a,b])
xlabel('x');
ylabel('f');
grid on;

figure('DefaultAxesFontSize',13);
hold on;
plot(freq,F_ds,'b','LineWidth',1);
plot(freq,F_nufft,'rx--','LineWidth',1);
xlim([0,200])
xlabel('Frequency (Hz)');
ylabel('|F|');
legend('DS','NUFFT');
grid on;
