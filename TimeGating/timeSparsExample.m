%% Example: time gating routines 
% The S parameters belong to the CST simulation of an uncoupled microstrip line consisting of three segments of different widths and lengths, 
% first and third segments of 3cm and Z_0=50\Omega, second segment of 4cm and Z_0=100\Omega (see "unCoupled_uStrip.cst" simulation file for further details).
clc,clear,close all,

%% Physical parameters 
c=3e8;                                              % speed of light
d1=41.5e-3;                                         % distance to first obstacle (from port plane)
d=123e-3;                                           % distance to second port (from port 1 to port 2 reference planes)
epsr=3;                                             % substrate's permittivity
he=1.52;                                            % substrate height
w1=3.78;                                            % first and third segment widths
w=0.6*3.78+0.4*1;                                   % ponderated average total line width
eps_1=(epsr+1)/2+(epsr-1)/2*(1+12*he/w1)^(-1/2);	% microstrip approximate effective permittivity for reflection segment 
eps_tr=(epsr+1)/2+(epsr-1)/2*(1+12*he/w)^(-1/2);	% microstrip approximate effective permittivity for transmission
t_1ref=2*d1/(c/sqrt(eps_1));                        % estimated time of first reflection
t_tr=d/(c/sqrt(eps_tr));                            % estimated traveling time from port to port
tmax=5e-9;                                          % maximum plotting time for impulse responses (doesn't affect the signal support)
tp_est=[tmax,t_1ref;tmax,t_tr;tmax,t_tr;tmax,t_1ref];

%% Parameters for the definition of the truncation window in time domain
tini_wrf=0.2e-9;                                    % reflection initial time for truncation window (adjust using the presented graphics)
tfin_wrf=0.6e-9;                                    % reflection final time for truncation window (adjust using the presented graphics)
tini_wtr=0.35e-9;                                   % transmission initial time for truncation window (adjust using the presented graphics)
tfin_wtr=0.85e-9;                                   % transmission final time for truncation window (adjust using the presented graphics)
twin=[tini_wrf,tfin_wrf;tini_wtr,tfin_wtr;tini_wtr,tfin_wtr;tini_wrf,tfin_wrf];
win='rectwin';                                 % truncation window name (supported: 'blackman', 'hamming','blackmanharris','kaiserbessel','rectwin')...

%% Data loading
curdir='C:\Users\usuario\Documents\MATLAB\Thesis\Measurements\';
[Spars,freqs,R]=readSpars([curdir 'unCopled_uStrip_FD.s2p']);         
plotSpars(freqs,Spars,'fase','','','linea desacoplada','CST Simulation');
% Spars_f=cleanSpars(Spars, freqs, 'RationalS', 'canonical', 11);                                   % fitting is not really benefical in this case
% plotSpars(freqs, {Spars(1,1,:),Spars_f(1,1,:)},'fase','cmp','S_{11} ','linea desacoplada',{'Simulated','Fitting'}); 
nfreqs=length(freqs);
Nper=10;                                             % Frequency domain zero padding multiplier (time domain interpoation factor, use with caution)

%% parameters for the pulse 
pulsekind='normal';
fc=(max(freqs)+min(freqs))/2;	% in case of using a modulated pulse
bw=2*(max(freqs)-min(freqs));	% bandwidth (assumed samplig frequency)
fre=[bw,fc];

%% Naive (direct) time gating. This is for comparison purposes only. 
x11=ifft(squeeze(Spars(1,1,:)));                                                                    % take IDFT of S11 as is
% figure, subplot(2,1,1);
% plot(abs(x11)); title('IFFT of S11'); xlabel('sample number'); ylabel('|S11[n]|');                % plot and determine sample number for cutting    
x21=ifft(squeeze(Spars(2,1,:)));                                                                    
% subplot(2,1,2);
% plot(abs(x21)); title('IFFT of S21'); xlabel('sample number'); ylabel('|S21[n]|');                  
x12=ifft(squeeze(Spars(1,2,:)));
x22=ifft(squeeze(Spars(2,2,:)));

Nfft=nfreqs; k=1:Nfft; %Nfft=1024;                                                          
nini_rf=3; nfin_rf=6; nini_tr=4; nfin_tr=7;	                             

X(1,1,:)=fft(x11(nini_rf:nfin_rf),Nfft); X(2,1,:)=fft(x21(nini_tr:nfin_tr),Nfft);
X(1,2,:)=fft(x12(nini_rf:nfin_rf),Nfft); X(2,2,:)=fft(x22(nini_tr:nfin_tr),Nfft);

% X11magdB=20*log10(abs(X11)); X21magdB=20*log10(abs(X21)); 
% X11phase=unwrap(angle(X11)); X21phase=unwrap(angle(X21));
% % figure, subplot(2,1,1);                                                        % Cut, take DFT (of few samples, interpolating to a big number of samples) 
% plot(k,X11magdB, k,X21magdB, [1 Nfft],[eS11 eS11], [1 Nfft],[eS21 eS21]); legend('S11', 'S21','idealS11','idealS21');	
% title('FFT of corrected S11 and S21'); xlabel('sample number'); ylabel('|S_{11}[k]|, |S_{21}[k]| (dB)'); 
% subplot(2,1,2);
% plot(k,X11phase, k,X21phase); legend('S11', 'S21'); 
% title('FFT of corrected S11 and S21'); xlabel('sample number'); ylabel('< S_{11}[k], < S_{21}[k]'); 

%% Time gating of Spars
newS=zeros(2,2,nfreqs);
cont=0;
for cont1=1:2
    for cont2=1:2
        cont=cont+1;
        newS(cont1,cont2,:)=tdSpars(freqs,squeeze(Spars(cont1,cont2,:)),{pulsekind,fre},{win,twin(cont,:)},Nper,tp_est(cont,:),false);
    end
end
plotSpars(linspace(min(freqs),max(freqs),size(newS,3)),newS, 'fase');           % plot corrected S parameters

Spars_fit=cleanSpars(newS, freqs, 'RationalS', 'canonical', 11);
idealSpars=ones(2,2,nfreqs);
idealSpars(1,1,:)=1/3*idealSpars(1,1,:); idealSpars(2,1,:)=(1-1/9)*idealSpars(2,1,:); % Squared due to the two impedance transitions 
idealSpars(1,2,:)=(1-1/9)*idealSpars(1,2,:); idealSpars(2,2,:)=1/3*idealSpars(2,2,:);

%plotSpars(freqs, {Spars,X,newS,Spars_fit,idealSpars},'fase','cmp','','linea desacoplada',{'Original','Naive','Timegated','Fitting','Ideal'}); 
% 
plotSpars(freqs, {Spars,newS},'fase','cmp','','linea desacoplada',{'Original','Timegated'}); 
plotSpars(freqs, {newS,idealSpars},'fase','cmp','','linea desacoplada',{'Timegated','Ideal'}); 
plotSpars(freqs, {newS,Spars_fit},'fase','cmp','','linea desacoplada',{'Timegated','Fitting'}); 
plotSpars(freqs, {X,idealSpars},'fase','cmp','','linea desacoplada',{'Naive','Ideal'}); 
plotSpars(freqs, {X,newS},'fase','cmp','','linea desacoplada',{'Naive','Timegated'}); 