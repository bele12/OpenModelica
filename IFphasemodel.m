%Simple demonstration for baseband phase-domain modeling
clear all
close all

%% define scenario parameters
c=3e8;     %speed of light
param.p = [0;0];     %position of radar
param.N_ant_rx = 4;

SNR = 0;  %rx snr [dB]

objects.p_0 = [0 20; 0 10]';   %object's initial (cartesian) positions
objects.v = [0 10; 0 0]';   %object's initial velocities in x and y directions in [m/s]
%objects.rcs = [10];       %object's RCS in dBm^2
objects.alpha = [1 1];   %amplitude of object, in this simple model, no radar equation is used
objects.M = length(objects.alpha);


%% signaling parameters
sig.B_sw = 1e9;    %sweep BW, [Hz]
sig.f_0 = 77e9;    %start frequency [Hz]
sig.T_sw = 24e-6;   %sweep time (of one sweep/ramp) [s]
sig.N_sw = 128;     %number of sweeps/ramps
sig.N_samp_per_ramp = 1024;  %number of fast-time samples (ADC samples per ramp)
sig.k = sig.B_sw/sig.T_sw; %sweep slope [Hz/s]
sig.lambda = c/(sig.f_0+sig.B_sw/2);

sig.f_s = 25e6;      %sampling frequency
sig.T_s = 1/sig.f_s;
sig.t_sw = 0:sig.T_s:(sig.N_samp_per_ramp-1)*sig.T_s; sig.t_sw = sig.t_sw(:);

param.ant_rx_spacing = c/sig.f_0/2;  %RX antenna spacing [m]

%% convert object parameters to channel parameters (d, v, AoA)
objects.d_0_true = zeros(objects.M,1);
objects.AoA = zeros(objects.M,1);
objects.v_r_true = zeros(objects.M, 1);

for m = 1:objects.M
    objects.d_0_true(m) = norm(param.p-objects.p_0(:,m));  %distance
    %calculate components of velocity
    objects.v_true_abs(m) = norm(objects.v(:, m)); %absolute velocity
    p_diff = objects.p_0(:,m)-param.p;  %vector from radar to object
    p_u = -p_diff/norm(p_diff);  %unit vector from object to radar
    if(objects.v_true_abs(m) > eps)  %prevent computation for static objects
        cos_phi = objects.v(:,m)'*p_u/objects.v_true_abs(m);  %cos of angle between velocity and radial velocity
        objects.v_r_true(m) = objects.v_true_abs(m)*cos_phi;  %radial velocity
    end
    objects.AoA(m) = atan2(p_diff(2), p_diff(1))-pi/2; %angle of arrival
end

%% generate IF signal
sig_IF = zeros(sig.N_samp_per_ramp, sig.N_sw, param.N_ant_rx);
%distance difference of arrival at RX antennas
delta_d = param.ant_rx_spacing .* sin(objects.AoA);

%base time matrix (same for all sweeps)
T_sw = sig.t_sw*ones(1,sig.N_sw);
%continued time matrix (continuous also along sweep dimension)
T_sw_cont = T_sw + sig.T_sw*repmat(0:(sig.N_sw-1), sig.N_samp_per_ramp, 1);

for i_ant = 1:param.N_ant_rx
    sig_temp = zeros(sig.N_samp_per_ramp, sig.N_sw);  %temporary to sum up object signals
    for m = 1:objects.M  %over objects
        d_0_i_m = objects.d_0_true(m);  %current distance of i-th antenna
        tau_0_i_m = (2*d_0_i_m + delta_d(m).*(i_ant-1))/c;   %RTDT of i-th antenna (start of first sweep)
        
        %compute phase of m-th IF signal
        Tau = tau_0_i_m + 2*objects.v_r_true(m)/c .* T_sw_cont;
        phi_m = 2*pi*sig.k.*Tau.*T_sw + 2*pi.*Tau.*sig.f_0 - pi .*sig.k.*Tau.^2;
            
        sig_temp = sig_temp + objects.alpha(m)*cos(phi_m);
        
    end  %over objects
    
    %calculate signal and noise power and generate noise
    sig_power = 10*log10(sum(abs(sig_temp(:).^2))./length(sig_temp(:)));
    noise_power = sig_power-SNR;
    sig_noise_temp = sqrt(10^(noise_power/10))*randn(size(sig_temp));
    
    %total IF signal
    sig_IF(:,:,i_ant) = sig_temp + sig_noise_temp;
    
end  %over antennae

%% resolution and basis vectors
limits.d_max = c*sig.f_s/4/abs(sig.k);
%limits.d_max = c*sig.N_samp_per_ramp/4/sig.B_sw;
limits.d_res = c/2/abs(sig.B_sw);
limits.v_max = c/(4*sig.f_0*sig.T_sw);  
limits.v_res = c/(2*sig.f_0*sig.T_sw*sig.N_sw);

%distance vector for first stage DFT
d_vec = linspace(0,1,sig.N_samp_per_ramp/2)*limits.d_max; 
%velocity vector
v_vec = [-sig.N_sw/2 : sig.N_sw/2-1].' ./sig.N_sw .* 2*limits.v_max; 

%% Range/Doppler processing
%windowing in fast-time
win_range = repmat(hann(sig.N_samp_per_ramp), 1, sig.N_sw, param.N_ant_rx);
sig_IF_win = sig_IF.*win_range;
%FFT over fast-time
S_range = 1/sqrt(sig.N_samp_per_ramp) * fft(sig_IF_win, sig.N_samp_per_ramp);
%cut out positive side of spectrum 
S_range = S_range(1:length(d_vec), :, :);

%windowing in slow-time/ramps
win_doppler = repmat(hann(sig.N_sw).', length(d_vec), 1, param.N_ant_rx);
S_range_win = S_range.*win_doppler;
%FFT over slow-time/ramps
S_RD = 1/sqrt(sig.N_sw) * fftshift(fft(S_range_win, sig.N_sw, 2),2);


%% plot
%IF signal
figure, hold on, grid on, box on
plot(sig.t_sw, sig_IF(:, 1, 1), 'b-')
xlabel('t [s]')
ylabel('r_{IF}(t)')
xlim([0 sig.T_sw])

%Range spectrum
figure, hold on, grid on, box on
%select RX antenna
S_range_plot = S_range(:, :, 1);
%scale
S_range_plot = S_range_plot./max(max(abs( S_range_plot )));
S_range_plot = 20*log10(abs(S_range_plot));
imagesc(1:sig.N_sw, d_vec, S_range_plot)
colormap(flipud(hot)), colorbar
axis([1 sig.N_sw min(d_vec) max(d_vec)])
xlabel('Ramp index')
ylabel('distance [m]')
title('After first DFT')

%Range-Doppler matrix
figure, hold on, grid on, box on
%select RX antenna
S_RD_plot = S_RD(:, :, 1);
%scale
S_RD_plot = S_RD_plot./max(max(abs( S_RD_plot )));
S_RD_plot = 20*log10(abs(S_RD_plot));
imagesc(v_vec, d_vec, S_RD_plot)
colormap(flipud(hot)), colorbar
axis([min(v_vec) max(v_vec) min(d_vec) max(d_vec)])
xlabel('velocity [m/s]')
ylabel('distance [m]')
title('Range-Doppler matrix')


