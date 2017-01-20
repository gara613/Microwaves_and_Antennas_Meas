%% Example: time gating routines 
clc,clear,close all,

%% Physical parameters 
% For further details on the example recreated herein see sec 6.2  http://cp.literature.agilent.com/litweb/pdf/5989-5723EN.pdf
% simulation data is obtained from CST design studio simulation
c=3e8;                                  % speed of light
d1=100e-3;                              % distance to the gated first capacitor (from port plane)
d=300e-3;                               % distance to second port (from port 1 to port 2 reference planes)
epsr=1;                                 % substrate's permittivity
t_rf=2*d1/(c/sqrt(epsr));               % estimated time of first reflection
t_tr=d/(c/sqrt(epsr));                  % estimated traveling time from port to port
tmax=40e-9;                             % maximum plotting time for impulse responses (doesn't affect the signal support)
tp_est=[tmax,t_rf; tmax,t_tr; tmax,t_tr; tmax,t_rf];

%% Parameters for the definition of the truncation window in time domain
tini_wrf=0.3e-9;                        % reflection initial time for truncation window (adjust using the presented graphics)
tfin_wrf=1.1e-9;                        % reflection final time for truncation window (adjust using the presented graphics)
tini_wtr=0.7e-9;                        % transmission initial time for truncation window (adjust using the presented graphics)
tfin_wtr=1.4e-9;                        % transmission final time for truncation window (adjust using the presented graphics)
twin=[tini_wrf,tfin_wrf; tini_wtr,tfin_wtr; tini_wtr,tfin_wtr; tini_wrf,tfin_wrf];
win='kaiserbessel';                          % truncation window name (supported: 'blackman', 'hamming','blackmanharris','kaiserbessel','rectwin')...

%% Data loading
curdir='C:\Users\usuario\Documents\MATLAB\Thesis\Calibraciones\';
[Spars,freqs,R]=readSpars([curdir 'tg3Lines_2Caps.s2p']);
plotSpars(freqs,Spars,'fase','','','linea desacoplada','CST Simulation');
nfreqs=length(freqs);
Nper=10;                                             % Frequency domain zero padding multiplier (time domain interpoation factor, use with caution)

%% parameters for the pulse 
pulsekind='normal';
fc=(max(freqs)+min(freqs))/2;	% in case of using a modulated pulse
bw=2*(max(freqs)-min(freqs));	% bandwidth (assumed samplig frequency)
fre=[bw,fc];

%% Time gating of Spars
newS=zeros(2,2,nfreqs);
cont=0;
for cont1=1:2
    for cont2=1:2
        cont=cont+1;
        newS(cont1,cont2,:)=tdSpars(freqs,squeeze(Spars(cont1,cont2,:)),{pulsekind,fre},{win,twin(cont,:)},Nper,tp_est(cont,:),false);
    end
end

% newS=timeGating(tdS,tini,tfin,win);                                           % Not sure if using one or two steps for this
% plotSpars(freqs,newS, 'fase');           % plot corrected S parameters
% 
[Sparsideal,freqs,R]=readSpars([curdir 'tg3Lines_2Caps_2ndRemoved.s2p']);           % correct value from simulation
% plotSpars(freqs,{newS,Sparsideal},'fase','cmp','comparison with ideal data','',{'Corrected','Ideal'});           % plot corrected S parameters

%% time gating for reflection only
Ref=squeeze(Spars(1,1,:));
newR=tdSpars(freqs,Ref,{pulsekind,fre},{win,twin(1,:)},Nper,tp_est(1,:),true);
normalizer=tdSpars(freqs,ones(size(Ref)).*exp(-1i*2*pi*freqs*0.7056e-9),{pulsekind,fre},{win,twin(1,:)},Nper,tp_est(1,:),true);
newRef(1,1,:)=newR./normalizer;
plotSpars(freqs,{Spars(1,1,:),newS(1,1,:),newRef,Sparsideal(1,1,:)},'fase','cmp','','comparison with ideal data',{'Measured','Corrected','Corrected normalized','Ideal'});           % plot corrected S parameters