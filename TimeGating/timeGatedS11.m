% S11tg = timeGatedS11(f,S11,tmin,tmax)
%
% (c) Javier Araque, november 2017

function S11tg = timeGatedS11(f,S11,tmin,tmax)

if true
  %warning('testing code');
  nn = ceil(0.1*length(f));
  
  deltaf = diff(f(1:2));
  
  f = [f(1)+(-nn:-1)*deltaf,f(:)',f(end)+(1:nn)*deltaf];
  S11 = S11(:).';
  
  S11fill1 = fliplr(S11(1:nn));
  S11fill1 = abs(S11fill1).*exp(-1j*angle(S11fill1))*exp(2j*angle(S11fill1(end)));
  
  S11fill2 = fliplr(S11(end-nn+1:end));
  S11fill2 = abs(S11fill2).*exp(-1j*angle(S11fill2))*exp(2j*angle(S11fill2(1)));
  
  S11 = [S11fill1,S11,S11fill2];
  
  filling = true;
else  
  filling = false;
end

wintype = 'rect';

if nargin == 0
  [f,S11,tmin,tmax] = testS11tg;
  test = true;
else
  test = false;
end

%check that f samples are uniformly spaced and f step
Df = diff(f);

if (max(Df) - min(Df))/mean(Df) > 1e-6
  error('f samples must be uniformly spaced');
else
  Df = mean(Df);
end

nf = length(f);
fmax = max(f);
fmin = min(f);
fmean = (fmax+fmin)/2;
hBW = (fmax-fmin)/2;

%check that time window specified is compatible with frequency sweep
if (tmax-tmin)*(fmax-fmin) < 10
  warning('Time window requested is shorter than pulse duration, results will be inaccurate');
  %test = true;
end

if rem(f(end),Df) > (Df*1e-3)
  %warning('Extended f grid does not include frequency 0, this causes some error with this implementation');
end

%complete grids using conjugate symmetry for S11
fext = -f(end):Df:f(end);
S11ext = zeros(size(fext));
S11ext(1:nf) = fliplr(conj(S11));
S11ext((end-nf+1):end) = S11;

next  = length(fext);

%compute exciting signal, a modulated Gaussian pulse, place +-3*sigma at band edges
%signalF = e.^(-linspace(-4.5,4.5,nf).^2);
%signalFext = zeros(size(fext));%with sine modulation to have zero mean
%signalFext(1:nf) = 0.5j*signalF;
%signalFext((end-nf+1):end) = -0.5*j*signalF;
signalFext = -0.5j*exp(-(fext-fmean).^2/(2*(hBW/3)^2)) + ...
      + 0.5j*exp(-(fext+fmean).^2/(2*(hBW/3)^2)); 

Dt = 1/(2*fmax);
Tmax = 1/Df;
tbase = linspace(0,Tmax-Dt,next);

timeSignal = ifft(ifftshift(S11ext.*signalFext));

if test
  subplot(211);
  timepulse = real(ifft(ifftshift(signalFext)));
  plot(tbase,timepulse/max(abs(timepulse)),'x-');
  hold on
  plot(tbase,timeSignal/max(abs(timeSignal)),'o-r');
  hold on
  legend('Normalized Input pulse','Normalized Time response');
  keyboard
end

%positions corresponding to the gated time

%negative times are measured w.r.t. Tmax
if tmin<0
  tmin = Tmax+tmin;
end

if tmax<0
  tmax = Tmax + tmax;
end

if tmin<=tmax
  tgpos = find((tbase>=tmin) & (tbase <=tmax));
else
  tgpos = find((tbase>=tmin) | (tbase <=tmax));  
end

%do actual time gating, use raised cosine window
%timeSignal((tbase<tmin) | (tbase>tmax)) = 0;
switch wintype
  case 'Gauss'%Gaussian window
    gatelength = tgpos(end)-tgpos(1);
    wf = GaussWindow(next,gatelength/3,(tgpos(end)+tgpos(1))/2);
    %wf = leoSigmoid(wf);
  case 'Kaiser'%Kaiser window with nulls at the borders of the gating interval
    tsamples = length(tgpos);
    wf = KaiserWindow(next,pi*sqrt((tsamples/2)^2-1),(tgpos(end)-tgpos(1))/2);
    %displace window to center of time gating interval
  case 'Tukey'%Tukey window, convolution of rectangle with raised cosine
    tmp = TukeyWindow(length(tgpos),0.1);
    wf = zeros(size(tbase));
    wf(tgpos) = tmp;
  case 'raisedCos' %raised cosine window
    wf = sin(pi*(tbase-tmin)/(tmax-tmin)).^2;
    wf((tbase<tmin) | (tbase>tmax)) = 0;
  case 'rect'%rectangular window
    wf = zeros(size(tbase));
    wf(tgpos) = 1;
end

S11tg = fftshift(fft(timeSignal.*wf));
S11tg = S11tg./signalFext;%eliminate effect of pulse shape
S11tg = S11tg((end-nf+1):end);

if test
  subplot(211);
  %plot(tbase,20*log10(abs(timeSignal)),'x-r');
  plot(tbase,wf,'-b','linewidth',2);legend('time response','window');
  
  hold on
  subplot(212);
  plot(f,20*log10(abs(S11)),'linewidth',2);
  hold on
  plot(f,20*log10(abs(S11tg)),'r--','linewidth',2);
  legend('Initial','With time gating');
  ylabel('|S11| (dB)');
  axis([f(1) f(end) -60 0]);
end

if filling
  S11tg = S11tg((nn+1):(end-nn));
end

function [f,S11,tmin,tmax] = testS11tg

nf = 631;
fmin = 100e6;
fmax = 5e9;
f = linspace(0,fmax,nf);
f = f(ceil(nf*fmin/fmax):end);

%our system is a TL segment of length 45cm with Z0, followed by a short 
%segment of 2cm with Z0*x, in cascade with a similar setup, which is finally terminated
%in z Z0 load.

Z0 = 50;
Lcon = 0.3;
Lelcon = Lcon/(3e8/f(1))*360;
Zcon = 0.9*Z0;
Lelcable = 0.43/(3e8/f(1))*360;

Zin = ZlineTerm(Zcon,Z0,f/f(1),Lelcon);
Zin = ZlineTerm(Z0,Zin,f/f(1),Lelcable);
Zin = ZlineTerm(Zcon,Zin,f/f(1),Lelcon);
Zin = ZlineTerm(Z0,Zin,f/f(1),Lelcable);

S11 = (Zin-Z0)./(Zin+Z0);

tmin = 2e-9;
tmax = 4e-9;


%Build a Tukey window of length N with relative transition length as given
function out = TukeyWindow(N,rtl)

Ntrans = ceil(N*rtl);

trans = sin(0.5*pi*(1:Ntrans)/Ntrans).^2;

out = ones(N,1);
out(1:Ntrans) = trans;
out((end-Ntrans+1):end) = fliplr(trans);


%build a Kaiser window with length N and width parameter beta
function out = KaiserWindow(N,beta,n0)

nbase = 0:(N-1);
out = besseli(0,beta*sqrt(1-(2*nbase/(N-1)-1).^2))./besseli(0,beta);
offset = round(n0-(N-1)/2);

out = shift(out,offset);

%build a Gauss window with length N, standard deviation sigma (in units) and centered at n0
function out = GaussWindow(N,sigma,n0)

nbase = 0:(N-1);
nbase = nbase - (N-1)/2;
out = exp(-nbase.^2/(2*sigma^2));

offset = round(n0-(N-1)/2);

out = shift(out,offset);

function out = leoSigmoid(in)

out = zeros(size(in));

pos1 = in>0 & in<=0.5;
out(pos1) = 0.5-sqrt(0.25-in(pos1).^2);

pos2 = in>0.5 & in<=1;
out(pos2) = 0.5+sqrt(0.25-(in(pos2)-1).^2);

pos3 = in>1;
out(pos3) = 1;