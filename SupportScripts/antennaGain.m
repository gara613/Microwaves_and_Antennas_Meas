%% Antenna_Gain
% from the S21 measurements of a reference antenna and the DUT antenna
% Friis equation: Pr=PtGtGr(lambda/4pir)^2, S21=V2-/V1+|V2+=0, S21^2 = Pr/Pt
% The S parameters come from to the CST simulation of two horn antennas with 20dBi intended gain, 
% separated at 3 and 10m respectively
clear all; close all;

%% Physical parameters 
c=3e8;                          % speed of light
r_near=3.67;                    % first case separation of antennas port references
r_far=10.67;                    % second case separation of antennas % farfield: 2D^2/lambda_min
t_tr_near=r_near/c;             % estimated traveling time from port to port plane 
t_tr_far=r_far/c;                                   
tmax=100e-9;                    % maximum plotting time for impulse responses (doesn't affect the signal support) 
tp_est_near=[tmax,t_tr_near];
tp_est_far=[tmax,t_tr_far];

%% Parameters for the definition of the truncation window in time domain
tini_tr_near=0.2*t_tr_near;     % initial time for truncation window (adjust using the presented graphics)
tfin_tr_near=2.0*t_tr_near;     % final time for truncation window 
tini_tr_far=0.4*t_tr_far;       
tfin_tr_far=1.6*t_tr_far;       
twin_near=[tini_tr_near,tfin_tr_near];
twin_far=[tini_tr_far,tfin_tr_far];
win='kaiserbessel';             % supported: 'blackman','hamming','blackmanharris','kaiserbessel','rectwin'

%% load data
thisdir='C:\Users\usuario\Documents\MATLAB\Thesis\Measurements\3D_cactus_01_11_2017';
[Spars_ref_near,freq_ref]=readSpars([thisdir, '\HornAntennaWR187G20dB_x2.s2p']);
[Spars_ref_far,freq_ref]=readSpars([thisdir, '\HornAntennaWR187G20dB_x2_far.s2p']);

%% parameters for the pulse 
nfreqs=length(freq_ref);
Nper=5;                             % Frequency domain zero padding multiplier (time domain interpoation factor)
pulsekind='normal';
fc=(max(freq_ref)+min(freq_ref))/2;	% in case of using a 'modulated' pulse
bw=2*(max(freq_ref)-min(freq_ref));	% bandwidth (assumed samplig frequency)
fre=[bw,fc];

%% Apply time gating to eliminate (most) reflection effects
S21_ref_near=squeeze(Spars_ref_near(2,1,:));
S21_ref_far=squeeze(Spars_ref_far(2,1,:));
S21_ref_near_tg=tdSpars(freq_ref,S21_ref_near,{pulsekind,fre},{win,twin_near},Nper,tp_est_near);
S21_ref_far_tg=tdSpars(freq_ref,S21_ref_far,{pulsekind,fre},{win,twin_far},Nper,tp_est_far);
%figure, plot(freq_ref,20*log10(abs(S21_ref_near)),'b',freq_ref,20*log10(abs(S21_ref_near_tg)),'r','linewidth',2); legend('meas','TGS')

%% Reference antenna calibration, assuming both antennas are equal in gain
Gtr_n_tg = abs(S21_ref_near_tg).*freq_ref*4*pi*r_near/3e8;
Gtr_f_tg = abs(S21_ref_far_tg).*freq_ref*4*pi*r_far/3e8;

%% DUT antenna gain calculation
%Gdut = (S21_dut*4*pi*r./lambda).^2./Gtr;

figure,plot(freq_ref,10*log10(abs(Gtr_n_tg)),'b',freq_ref,10*log10(abs(Gtr_f_tg)),'r','linewidth',2);
title('Reference antenna Gain'); legend('3m','10m')
xlabel('frequencies'); ylabel('Gain (dB)'); grid on;
% plot(freq,10*log10(Gdut),'linewidth',2);
% title('Cactus 3D monopole antenna Gain'); 
% xlabel('frequencies'); ylabel('Gain (dB)'); grid on;
% compare to the difference in dB: S21ref-S21dut