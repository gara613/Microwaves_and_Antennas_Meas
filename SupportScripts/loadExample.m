% Results reading and processing (based on the R&S FSH8 measurements, change accordingly in case of CST simulations)
% First step (not avoided yet) is to open the results in FSHview software,
% then right click on figure -> copy data -> paste on txt file, save and exit
%
% all routines can be downloaded from: https://github.com/gara613

%% first step is to load data
baseroute = 'C:\Users\usuario\Documents\MATLAB\Thesis\Measurements\AntennasDavid\';
[antenna1, freqsr, Zref] = readSpars([baseroute 'antenna1.txt']);	
[antenna2, ~, ~] = readSpars([baseroute 'antenna2.txt']);   % assume same sampling and reference impedance in all measurements
[antenna3, ~, ~] = readSpars([baseroute 'antenna3.txt']);   %
plotSpars(freqsr,{antenna1,antenna2,antenna3},'','cmp','','Línea',{'Ant_1','Ant_2','Ant_3'});

%% Time gating of S parameters
pulsekind='modulated';
nfreqs=length(freqsr);
fc=(max(freqsr)+min(freqsr))/2;	% in case of using a modulated pulse
bw=2*(max(freqsr));%-min(freqsr));	% pulse bandwidth
fre=[bw,fc];

c=3e8;
dl=1.2e-2;          %coax length
l=7e-2;             %distance between ports
eps_r=4.4;
epsr_coax=2.1;      
t_1ref=abs(dl)*sqrt(epsr_coax)/c;
t_tr=l*sqrt((eps_r+1)/2)/c + 2*abs(dl)*sqrt(epsr_coax)/c; % only an estimate as effective travel distance is not known beforehand

est_pulDur = 1/(max(freqsr)+min(freqsr));
tini_wrf=-200*t_1ref;           % reflection initial time for truncation window
tfin_wrf=200*t_1ref;            % reflection final time 
tini_wtr=-15*t_tr;              % transmission initial time 
tfin_wtr=25*t_tr;               % transmission final time 
twin=[tini_wrf,tfin_wrf; tini_wtr,tfin_wtr; tini_wtr,tfin_wtr; tini_wrf,tfin_wrf];
tmax=50*t_tr;%50*t_tr;
tp_est=[tmax,t_1ref; tmax,t_tr; tmax,t_tr; tmax,t_1ref];
win='kaiserbessel'; 
Nper=5;

% tfin=1/(10*(max(freqsr)+min(freqsr)))+tini;

antenna1_tg=zeros(2,2,nfreqs);
antenna2_tg=zeros(2,2,nfreqs);
antenna3_tg=zeros(2,2,nfreqs);
antenna1_JA=zeros(2,2,nfreqs);
antenna2_JA=zeros(2,2,nfreqs);
antenna3_JA=zeros(2,2,nfreqs);

cont=0;
for cont1=1:2
    for cont2=1:2
        cont=cont+1;
        normalizer=tdSpars(freqsr,ones(length(freqsr),1),{pulsekind,fre},{win,twin(cont,:)},Nper,tp_est(cont,:),false);
        antenna1_tg(cont1,cont2,:)=tdSpars(freqsr,squeeze(antenna1(cont1,cont2,:)),{pulsekind,fre},{win,twin(cont,:)},Nper,tp_est(cont,:),true)./normalizer;
        antenna2_tg(cont1,cont2,:)=tdSpars(freqsr,squeeze(antenna2(cont1,cont2,:)),{pulsekind,fre},{win,twin(cont,:)},Nper,tp_est(cont,:),false)./normalizer;
        antenna3_tg(cont1,cont2,:)=tdSpars(freqsr,squeeze(antenna3(cont1,cont2,:)),{pulsekind,fre},{win,twin(cont,:)},Nper,tp_est(cont,:),false)./normalizer;
        antenna1_JA(cont1,cont2,:) = timeGatedS11(freqsr,squeeze(antenna1(cont1,cont2,:)),twin(cont,1),twin(cont,2));       
        antenna2_JA(cont1,cont2,:) = timeGatedS11(freqsr,squeeze(antenna2(cont1,cont2,:)),twin(cont,1),twin(cont,2));       
        antenna3_JA(cont1,cont2,:) = timeGatedS11(freqsr,squeeze(antenna3(cont1,cont2,:)),twin(cont,1),twin(cont,2));       
    end
end
plotSpars(freqsr,{antenna1,antenna1_tg,antenna1_JA},'','cmp','','Antenna 1',{'meas','tg corr GR','tg corr JA'})
plotSpars(freqsr,{antenna2,antenna2_tg,antenna2_JA},'','cmp','','Antenna 2',{'meas','tg corr GR','tg corr JA'})
plotSpars(freqsr,{antenna3,antenna3_tg,antenna3_JA},'','cmp','','Antenna 3',{'meas','tg corr GR','tg corr JA'})

plotSpars(freqsr,{antenna1_tg,antenna2_tg,antenna3_tg},'','cmp','','Diversity Antenna Comparison',{'Ant1','Ant2','Ant3'})