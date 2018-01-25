% TL characterization
% routine for the extraction of the effective permitivity of a transmission
% line from Complex Reflection and Transmission Coefficients
close all; clc; clear all;

%% physical parameters
l=5e-2;     %line length
dl=-1.7e-2; %-0.86e-2;      % de-embedding length (negative to account for movement inwards the line from the port plane) 
mur_coax=1;                 % (0.78cm = dielectric physical length, 1.15cm inner conductor length)
epsr_coax=2.1; 
mur=1;      
c=3e8;
w=3.78e-3;	%microstrip width
h=1.52e-3;	%substrate thickness
eps_r=3;    %lower substrate relative permittivity
dc=3e-3;	%cover substrate thickness (2.75mm sample1, 3.4mm sample 2, 3mm CST simulations)
%w=w/h;     %normalized microstrip width

%% load data
baseroute='C:\Users\usuario\Documents\MATLAB\Thesis\Measurements\Mediciones_lineaUstrip_2Subs\';
%[ref, freqsr, Zref] = readSpars([baseroute 'line5cmAir.txt']);  % S parameters matrix of reference line
%[mat1, freqs1, ~] = readSpars([baseroute 'line5cmMat1.txt']);   % first and second material dielectric cover
%[mat2, freqs2, ~] = readSpars([baseroute 'line5cmMat2.txt']);
%plotSpars(freqsr,{ref,mat1,mat2},'','cmp','','Línea',{'ref','mat1','mat2'});

[ref, freqsr, Zref] = readSpars([baseroute 'CST_simulation\TwoSubsLine_0.03_8GHz_1.s2p']);  % S parameters matrix of reference line
[mat1, freqs1, ~] = readSpars([baseroute 'CST_simulation\TwoSubsLine_0.03_8GHz_34.s2p']);   % only S11 and S21 are of interest, hence S12=S22=0
[mat2, freqs1, ~] = readSpars([baseroute 'CST_simulation\TwoSubsLine_0.03_8GHz_50.s2p']);   
plotSpars(freqsr,{ref,mat1,mat2},'','cmp','','Línea',{'air','Dk=3','Dk=5'});

%% Time gating of S parameters
pulsekind='normal';
nfreqs=length(freqsr);
fc=(max(freqsr)+min(freqsr))/2;	% in case of using a modulated pulse
bw=2*(max(freqsr)-min(freqsr));	% bandwidth (assumed samplig frequency)
fre=[bw,fc];
t_1ref=abs(dl)*sqrt(epsr_coax)/c;
t_tr=l*sqrt(eps_r)/c + 2*abs(dl)*sqrt(epsr_coax)/c;
tini_wrf=-10*t_1ref;	% reflection initial time for truncation window
tfin_wrf=10*t_1ref;     % reflection final time 
tini_wtr=-10*t_tr;      % transmission initial time 
tfin_wtr=10*t_tr;       % transmission final time 
twin=[tini_wrf,tfin_wrf; tini_wtr,tfin_wtr; tini_wtr,tfin_wtr; tini_wrf,tfin_wrf];
tmax=50*t_tr;
tp_est=[tmax,t_1ref; tmax,t_tr; tmax,t_tr; tmax,t_1ref];
win='kaiserbessel'; 
Nper=5;
ref_tg=zeros(2,2,nfreqs);
mat1_tg=zeros(2,2,nfreqs);
mat2_tg=zeros(2,2,nfreqs);

% cont=0;
% for cont1=1:2
%     for cont2=1:2
%         cont=cont+1;
%         ref_tg(cont1,cont2,:)=tdSpars(freqsr,squeeze(ref(cont1,cont2,:)),{pulsekind,fre},{win,twin(cont,:)},Nper,tp_est(cont,:),false);
%         mat1_tg(cont1,cont2,:)=tdSpars(freqsr,squeeze(mat1(cont1,cont2,:)),{pulsekind,fre},{win,twin(cont,:)},Nper,tp_est(cont,:),false);
%         mat2_tg(cont1,cont2,:)=tdSpars(freqsr,squeeze(mat2(cont1,cont2,:)),{pulsekind,fre},{win,twin(cont,:)},Nper,tp_est(cont,:),false);
%     end
% end
% plotSpars(freqsr,{ref,ref_tg},'fase','cmp','','uStrip Tx line for material characterization',{'meas','tg corr'})
% plotSpars(freqsr,{mat1,mat1_tg},'fase','cmp','','uStrip Tx line for material characterization',{'meas','tg corr'})


%% retrieve effectivel Tx line parameters
u=w/h;
effepsr_model_line_air=(eps_r+1)/2+(eps_r-1)/2*(1+12/u)^(-1/2);	% "classic" microstrip approximate effective permittivity
% A = 1 + 1/49*log((u^4+(u/52)^2)/(u^4+0.432)) + 1/18.7*log(1+(u/18.1)^3);
% B = 0.564*((eps_r-0.9)/(eps_r+3))^0.053;
% effepsr_model_line_air = 0.5*(eps_r+1) + 0.5*(eps_r-1)*(1+10*h/w)^(-A*B) % microstrip Hammerstad and Jensen approximate effective permittivity
fu=6+(2*pi-6)*exp(-(30.666*h/w)^0.7528);
%Z=50/sqrt(eps_r)*1/(w/h+2.42-0.44*h/w+(1-h/w)^6)
Zf0=377;
Zl=Zf0/(2*pi*sqrt(eps_r))*log(fu*h/w + sqrt(1+(2*h/w)^2));
G=pi^2/12*(eps_r-1)/effepsr_model_line_air*sqrt(2*pi*Zl/Zf0);
fp=Zl/(2*4*pi*1e-7*h);
effepsr_model_line_air=eps_r-(eps_r-effepsr_model_line_air)./(1+G.*(freqsr/fp).^2);
        
effPars_line_air=txLineParsfromS(ref,freqsr,1,mur,l,Zref,dl,epsr_coax,mur_coax);
effPars_line_air_tg=txLineParsfromS(ref_tg,freqsr,1,mur,l,Zref,dl,epsr_coax,mur_coax);
effPars_line_mat1=txLineParsfromS(mat1,freqsr,1,mur,l,Zref,dl,epsr_coax,mur_coax);
effPars_line_mat1_tg=txLineParsfromS(mat1_tg,freqsr,1,mur,l,Zref,dl,epsr_coax,mur_coax);
effPars_line_mat2=txLineParsfromS(mat2,freqsr,1,mur,l,Zref,dl,epsr_coax,mur_coax);
effPars_line_mat2_tg=txLineParsfromS(mat2_tg,freqsr,1,mur,l,Zref,dl,epsr_coax,mur_coax);

%% plots of effective permittiviy, tandelta and characteristic impedance
figure, 
subplot(2,1,1)
plot(freqsr,real(effPars_line_air.gamma{1}),freqsr,real(effPars_line_air.gamma{2}),freqsr,real(effPars_line_air.gamma{3})); 
title('\gamma'); legend('JA','cmp_cosh','cmp_exp'); ylabel('real');xlabel('freq');
subplot(2,1,2)
plot(freqsr,imag(effPars_line_air.gamma{1}),freqsr,imag(effPars_line_air.gamma{2}),freqsr,imag(effPars_line_air.gamma{3})); 
title('\gamma'); legend('JA','cmp_cosh','cmp_exp'); ylabel('imag');xlabel('freq');

figure
subplot(2,1,1)
plot(freqsr,effepsr_model_line_air.*ones(size(freqsr)),freqsr,effPars_line_air.epsreff,freqsr,effPars_line_air_tg.epsreff,...
    freqsr,effPars_line_mat1.epsreff,freqsr,effPars_line_mat1_tg.epsreff,...
    freqsr,effPars_line_mat2.epsreff,freqsr,effPars_line_mat2_tg.epsreff,'linewidth',2); 
title('Effective permittivity'); legend('air model','air','air tg','mat1','mat1 tg','mat2','mat2 tg'); 
ylabel('\epsilon_{reff}');xlabel('freq'); axis([min(freqsr), max(freqsr), 0, 2*max(effepsr_model_line_air)])
subplot(2,1,2)
plot(freqsr,abs(effPars_line_air.tandelta),freqsr,abs(effPars_line_mat1.tandelta),freqsr,abs(effPars_line_mat2.tandelta));
title('loss tangent');legend('air','mat1','mat2'); ylabel('tan\delta');xlabel('freq'); axis([min(freqsr), max(freqsr), 0, 1])

figure
subplot(2,1,1)
plot(freqsr,real(effPars_line_air.Z),freqsr,real(effPars_line_mat1.Z),freqsr,real(effPars_line_mat2.Z)); 
title('Effective impedance');legend('air','mat1','mat2');ylabel('real');xlabel('freq');
subplot(2,1,2)
plot(freqsr,imag(effPars_line_air.Z),freqsr,imag(effPars_line_mat1.Z),freqsr,imag(effPars_line_mat2.Z)); 
title('Effective impedance');legend('air','mat1','mat2');ylabel('imag');xlabel('freq');
 
%% Comparison to the characteristic impedance of a covered microstrip line with a dielectric overlay (Barbuto, Alú, Bilotti, Toscano, Vegni)
