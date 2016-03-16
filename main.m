
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

pos_initial = 10; % [m]
vel_initial = 0; % [m/s]

num_samples = 1024;

force_constant = 0; %[N]
freq_exciting = 6.5; % [Hz]
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
% Transfer function model

H_pos_analytic_num = [1];
H_pos_analytic_den = [mass, damping, spring];
H_pos_analytic = tf(H_pos_analytic_num, H_pos_analytic_den);
H_vel_analytic_num = [1 0];
H_vel_analytic_den = [mass, damping, spring];
H_vel_analytic = tf(H_vel_analytic_num, H_vel_analytic_den);

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
[x_iden,X_iden,freq_iden] = musin_type(sample_frequency,0,40,num_samples,'r','c','f',1,1);

%%
% Excite system with identification signal

force_iden = x_iden';
pos_iden = zeros(size(t_sampled));
pos_iden(1) = pos_initial;
vel_iden = zeros(size(t_sampled));
vel_iden(1) = vel_initial;
accel_iden = zeros(size(t_sampled));

for n_time = (1:length(t_sampled)-1)
    [t_sim,x_sim] = ode23tb (@(t,x) system_mkd(t,x,mass,spring,damping,force_iden(n_time),@(t,x) 0), [t_sampled(n_time), t_sampled(n_time+1)], [pos_int(n_time); vel_int(n_time)]);
    pos_iden(n_time+1) = x_sim(end,1);
    vel_iden(n_time+1) = x_sim(end,2);
end
for n_time = 1:length(t_sampled)
    state_current = system_mkd(t_sampled(n_time),[pos_int(n_time),vel_int(n_time)],mass,spring,damping,force_iden(n_time),@(t,x) 0);
    accel_iden(n_time) = state_current(2);
end

[y_sim,t_sim] = lsim(H_pos_analytic,force_iden,t_sampled,10);
pos_iden = y_sim';

%%
% Fourier transforms
NFFT = 2^nextpow2(num_samples);
pos_s = fft(pos_int,NFFT)/num_samples;
pos_exciting_s = fft(pos_exciting,NFFT)/num_samples;
pos_iden_s = fft(pos_iden,NFFT)/num_samples;
force_iden_s = fft(force_iden,NFFT)/num_samples;
transf_function = pos_iden_s ./ force_iden_s;
freq_transform = sample_frequency/2*linspace(0,1,NFFT/2+1);
range = 1:NFFT/2+1;
range = 1:floor(length(freq_transform)*3.9/5);

%%
% Fiting of the FRF

[Bn,An,Bls,Als,Bls2,Als2] = nllsfdi(transf_function(range)',freq_transform(range)',ones(size(freq_transform(range)')),2,0,0,100,1e-6,0,'c',sample_frequency);

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
plot(freq_transform(range),2*abs(pos_iden_s(range)));
title('Single-Sided Amplitude Spectrum of pos_{identification}(t)');
xlabel('Frequency (Hz)');
ylabel('|pos_{identification}(s)|');
subplot(3,1,2);
hold on
plot(freq_transform(range),2*abs(force_iden_s(range)));
plot(freq_iden,abs(X_iden),'r');
%axis([0 40 0 1.1]);
title('Single-Sided Amplitude Spectrum of force_{identification}(t)');
xlabel('Frequency (Hz)');
ylabel('|force_{identification}(s)|');
subplot(3,1,3);
plot(freq_transform(range),2*abs(transf_function(range)));
%axis([0 40 0 1e-4]);
title('Transfer function');
xlabel('Frequency (Hz)');
ylabel('|H(s)|');

%%
% Bode plot of the system, analytical model
figure(7);
hold on
[mag_a,phase_a,wout_a] = bode(H_pos_analytic);
[mag_f,phase_f,wout_f] = bode(tf(Bn,An),2*pi*freq_transform(range));
subplot(2,1,1);
plot(wout_a',db(reshape(mag_a,1,length(mag_a))),'b', ...
    wout_f',db(reshape(mag_f,1,length(mag_f))),'g', ...
    2*pi*freq_transform(range),db(transf_function(range)),'r');
set(gca,'xscale','log');
legend('Analytic model','Fit with nllsfdi','Measured','Location','SouthWest');
subplot(2,1,2);
plot(wout_a,reshape(phase_a,1,length(phase_a)),'b', ...
    wout_f,reshape(+1*phase_f,1,length(phase_f)),'g', ...
    2*pi*freq_transform(range),(180/pi)*unwrap(angle(transf_function(range))),'r');
%legend('Analytic model','Fit with nllsfdi');
set(gca,'xscale','log');
title('Bode analysis of Position/Force analytic plant of MKB system');
title('Bode analysis of fitted solution');
figure(8);
bode(H_vel_analytic);
title('Bode analysis of Velocity/Force analytic plant of MKB system');

%%
% Close unwatched figures
close(1);
close(2);
close(3);
close(4);
close(5);
close(6);
% close(7);
close(8);
