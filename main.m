
%%
% Clean workspace

close all
clear all
clc

%%
% Add toolbox
%addpath('../Toolboxes/Matlab_Swevers');

%%
% Initialitation of variables

mass = .03; % [kg]
spring = 50; % [N/m]
damping = 0.1; % [Ns/m]
t_max = 10; % [s]

pos_initial = 0; % [m]
vel_initial = 0; % [m/s]

num_samples = 1001;

force_constant = 0; %[N]
freq_exciting = 8; % [Hz]
amp_exciting = 100; % [N]

%%
% time steps

t_sampled = linspace(0,t_max, num_samples);
sample_period = mean(diff(t_sampled));
sample_frequency = 1/sample_period;

%%
% Resonance frequency

freq_natural_analytic = sqrt(spring/mass)/(2*pi);

%%
% Integration of the sytem (non-exited)

pos_int = zeros(size(t_sampled));
pos_int(1) = pos_initial;
vel_int = zeros(size(t_sampled));
vel_int(1) = vel_initial;
accel_int = zeros(size(t_sampled));

for n_time = (1:length(t_sampled)-1)
    [t_sim,x_sim] = ode45 (@(t,x) system_mkd(t,x,mass,spring,damping,force_constant,@(t,x) 0), [t_sampled(n_time), t_sampled(n_time+1)], [pos_int(n_time); vel_int(n_time)]);
    pos_int(n_time+1) = x_sim(end,1);
    vel_int(n_time+1) = x_sim(end,2);
end
for n_time = 1:length(t_sampled)
    state_current = system_mkd(t_sampled(n_time),[pos_int(n_time),vel_int(n_time)],mass,spring,damping,force_constant,@(t,x) 0);
    accel_int(n_time) = state_current(2);
end

%%
% Exciting system - one frequency

force_exciting = amp_exciting*sin(2*pi*freq_exciting*t_sampled);
pos_exciting = zeros(size(t_sampled));
pos_exciting(1) = pos_initial;
vel_exciting = zeros(size(t_sampled));
vel_exciting(1) = vel_initial;
accel_exciting = zeros(size(t_sampled));

for n_time = (1:length(t_sampled)-1)
    [t_sim,x_sim] = ode45 (@(t,x) system_mkd(t,x,mass,spring,damping,force_exciting(n_time),@(t,x) 0), [t_sampled(n_time), t_sampled(n_time+1)], [pos_int(n_time); vel_int(n_time)]);
    pos_exciting(n_time+1) = x_sim(end,1);
    vel_exciting(n_time+1) = x_sim(end,2);
end
for n_time = 1:length(t_sampled)
    state_current = system_mkd(t_sampled(n_time),[pos_int(n_time),vel_int(n_time)],mass,spring,damping,force_exciting(n_time),@(t,x) 0);
    accel_exciting(n_time) = state_current(2);
end

%%
% Generate identification excitation
[x_iden,X_iden,freq_iden,Xt_iden,freqt_iden] = musint2(sample_frequency,0,40,num_samples,'r','c','f',1,1);

%%
% Excite system with identification signal

force_iden = x_iden';
pos_iden = zeros(size(t_sampled));
pos_iden(1) = pos_initial;
vel_iden = zeros(size(t_sampled));
vel_iden(1) = vel_initial;
accel_iden = zeros(size(t_sampled));

for n_time = (1:length(t_sampled)-1)
    [t_sim,x_sim] = ode45 (@(t,x) system_mkd(t,x,mass,spring,damping,force_iden(n_time),@(t,x) 0), [t_sampled(n_time), t_sampled(n_time+1)], [pos_int(n_time); vel_int(n_time)]);
    pos_iden(n_time+1) = x_sim(end,1);
    vel_iden(n_time+1) = x_sim(end,2);
end
for n_time = 1:length(t_sampled)
    state_current = system_mkd(t_sampled(n_time),[pos_int(n_time),vel_int(n_time)],mass,spring,damping,force_iden(n_time),@(t,x) 0);
    accel_iden(n_time) = state_current(2);
end

%%
% Fourier transforms
NFFT = 2^nextpow2(num_samples);
pos_s = fft(pos_int,NFFT)/num_samples;
pos_exciting_s = fft(pos_exciting,NFFT)/num_samples;
pos_iden_s = fft(pos_iden,NFFT)/num_samples;
force_iden_s = fft(force_iden,NFFT)/num_samples;
transf_function = pos_iden_s ./ force_iden_s;
freq_transform = sample_frequency/2*linspace(0,1,NFFT/2+1);

%%
% Transfer function model

H_analytic_num = [1];
H_analytic_den = [mass, -damping, spring];
H_analytic = tf(H_analytic_num, H_analytic_den);

%%
% Plots

% Time domain plot
figure(1);
hold on
subplot(3,1,1);
plot(t_sampled,pos_int);
title('Integration of MBK system with ODE');
xlabel('Time [s]');
ylabel('Position [m]');
subplot(3,1,2);
plot(t_sampled,vel_int);
xlabel('Time [s]');
ylabel('Velocity [m/s]');
subplot(3,1,3);
plot(t_sampled,accel_int);
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');

% Frequency domain plot
figure(2);
plot(freq_transform,2*abs(pos_s(1:NFFT/2+1)));
title('Single-Sided Amplitude Spectrum of pos(t)');
xlabel('Frequency (Hz)');
ylabel('|pos(s)|');

% Time domain plot
figure(3);
hold on
subplot(3,1,1);
plot(t_sampled,pos_exciting);
title('Integration of MBK system with ODE');
xlabel('Time [s]');
ylabel('Position [m]');
subplot(3,1,2);
plot(t_sampled,vel_exciting);
xlabel('Time [s]');
ylabel('Velocity [m/s]');
subplot(3,1,3);
plot(t_sampled,accel_exciting);
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');

% Frequency domain plot
figure(4);
plot(freq_transform,2*abs(pos_exciting_s(1:NFFT/2+1)));
title('Single-Sided Amplitude Spectrum of pos_{exciting}(t)');
xlabel('Frequency (Hz)');
ylabel('|pos_{exciting}(s)|');

% Time domain plot
figure(5);
hold on
subplot(3,1,1);
plot(t_sampled,pos_iden);
title('Integration of MBK system with ODE');
xlabel('Time [s]');
ylabel('Position [m]');
subplot(3,1,2);
plot(t_sampled,vel_iden);
xlabel('Time [s]');
ylabel('Velocity [m/s]');
subplot(3,1,3);
plot(t_sampled,accel_iden);
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');

% Frequency domain plot
figure(6);
subplot(3,1,1);
plot(freq_transform,2*abs(pos_iden_s(1:NFFT/2+1)));
title('Single-Sided Amplitude Spectrum of pos_{identification}(t)');
xlabel('Frequency (Hz)');
ylabel('|pos_{identification}(s)|');
subplot(3,1,2);
hold on
plot(freq_transform,2*abs(force_iden_s(1:NFFT/2+1)));
plot(freq_iden,abs(X_iden),'r');
plot(freqt_iden,abs(Xt_iden),'g');
title('Single-Sided Amplitude Spectrum of force_{identification}(t)');
xlabel('Frequency (Hz)');
ylabel('|force_{identification}(s)|');
subplot(3,1,3);
plot(freq_transform,2*abs(transf_function(1:NFFT/2+1)));
title('Transfer function');
xlabel('Frequency (Hz)');
ylabel('|H(s)|');

%%
% Bode plot of the system, analytical model
figure(7);
bode(H_analytic);
title('Bode analysis of the analytic plant of MKB system');

%%
% Close unwatched figures
% close(1);
% close(2);
% close(3);
% close(4);
close(5);
close(6);
close(7);
