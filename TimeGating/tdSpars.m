% function h=tdSpars(freqs,Spar,pulse,varargin)
% Plots the time domain S parameter pulse response and allows to apply a
% window for time gating, then returns the corrected S parameter
% Inputs:
%   - freqs: frequency vector 
%   - Spar: S parameter vector Nfreqs x 1 
%   - pulse: cell array with
%       - string to indicate the pulse type: 'normal' (baseband), 'modulated' (gaussian), 'sinc'
%       - bandiwdth of the pulse (inverse of variance)
%       - central frequency of the pulse (only relevant for modulated pulse)
%   - varargin: possible fields include:
%       - {win,[tini,tfin]}: cell with: window kind (string), initial and final times for windowing
%       - Npad: multiplier for zero padding (time domain interpolation factor)
%       - [tmax,t_1ref,t_tr]: maximum time for display, estimated time of first reflection/arrival from port 1 to 1/2
%       - shouwpuls: boolean flag to enable time domain graphics
% Outputs: structure with fields
%       Smod: Gated S parameter
%
% Caution (Masking effect): using gating about a response that occurs after a substantial reflection can lead to quite erroneous results 
% Germán Augusto Ramírez Arroyave 
% CMUN - Universidad Nacional de Colombia 2016

function Smod=tdSpars(freqs,Spars,pulsepar,varargin)
    %% input reading
    r=length(Spars);
   	N=length(freqs); 
    if r~=N,	error('freqs and S pars sampling does not coincide');   end
 
    pulse=pulsepar{1};
    bw=pulsepar{2}(1);
    fc=pulsepar{2}(2);
   
    % time delay from S pars, could be used as delay time for the pulse (\tau_d(f))
%	tau_f=-unwrap(angle(Spars))./(2*pi*freqs);
%   td=5*abs(mean(tau_f(~isnan(tau_f))));
%   td=abs(max(tau_f));     
%     if false % used for debug only
%         figure; plot(freqs,tau_f,'b',freqs,ones(size(freqs))*td,'r','linewidth',2); grid on;
%         title('delay as a function of frequency'); xlabel('frequency (Hz)'); ylabel('delay (s)');
%     end
    
    %% optional arguments
    % truncation window
	win='rectwin';                          
    tini_w=0;
	tfin_w=(N-1)/(max(freqs)-min(freqs));     % time span (non periodic representable time)
    if length(varargin)>=1
        win=varargin{1}(1);
        tini_w=varargin{1}{2}(1);
        tfin_w=varargin{1}{2}(2);
    end
    % multiplier to get a new number of samples 
    Nper=1;                                 
	if length(varargin)>=2
        Nper=varargin{2}(1);
	end
    % variables used only for plotting: final time and estimated time of arrival for the first reflection/transmission
	if length(varargin)>=3
        tmax=varargin{3}(1);
        t_est=varargin{3}(2);
	end    
	showpuls=false;
    if length(varargin)>=4
        showpuls=varargin{4};
    end   
   
    %% extraction of signal parameters
    fs=2*max(freqs);                        % sampling frequency (assumed as twice the maximum representable frequency)
	deltaF=(max(freqs)-min(freqs))/(N-1);	% frequency step
    Nz=round(min(freqs)/deltaF);            % number of zeros for padding the initial part of the spectrum
    Nzp=(Nper-1)*(N+Nz);                    % number of zeros for padding at the end (to increase time resolution)
	Ns=2*(N+Nz+Nzp);                        % number of samples for the extended signal in frequency and interpolated time domain
    t=-(Ns-1)/(2*Nper*fs):1/(Nper*fs):(Ns-1)/(2*Nper*fs);       % time range is related to span as: (N-1)/span = 1/deltaF
    Nfft=Ns;%2^nextpow2(Ns);                 
    
    %% create pulse
    x=gaussPulse(t,fc,bw,Nfft,pulse);	% pulse time vector contains Nper times the original number of samples 
    f=linspace(-Nper*fs/2,Nper*fs/2,Ns);
    if showpuls
        figure, indplot=abs(t)<tmax;
        subplot(2,1,1);plot(t(indplot),x.time(indplot),'linewidth',2); grid on; title('time domain pulse'); xlabel('time (s)'); ylabel('amplitude');     
        subplot(2,1,2);plot(f,10*log10(abs(fftshift(x.freq))),'linewidth',2); grid on; title('frequency domain pulse'); xlabel('frequency (Hz)'); ylabel('magnitude (dB)');
    end

    %% pulse response for current S parameter (reflection/transmission)
	y=pulseRespfromSpar(Nz,Spars,Nzp,x);
   
    %% window    
    ind_w=find(t>tini_w & t<tfin_w);
    w=retrieveWin(Ns,ind_w,win);
    
    %% Time gated S parameters in the specified window
	hw=y.*w;
	Smod=fft(hw,Ns)./x.freq;        % get back S parameters removing the pulse effect              
	Smod=Smod(Nz+1:N+Nz);           % recover S parameters selecting the correct N samples              

    %% Efficient implementation completely in the frequency domain (not tested yet)
%     Smod_f = freqDomTimeGating(Nz,Spars,Nzp,x,w);
%     idx=fix(length(Smod)/2+1) : fix(length(Smod)/2)+N;
%     figure, plot(freqs,20*log10(abs(Spars)), freqs,20*log10(abs(Smod_f(idx)))); % ill behaved... 
%     Smod_f = Smod_f(Nz+1:N+Nz);

    %% Time domain plotting 
	if showpuls
        indplot=(abs(t)<tmax);%-td;
        tplot=t(indplot);

        figure; subplot(2,1,1);
        plot(tplot,real(y(indplot)),'linewidth',2); hold on;
        
        plot(tplot,max(real(y))*w(indplot),'g','linewidth',2); 
        legstr={'(scl) window','measured'}; 
        if exist('t_est','var')
            hold on; plot([t_est t_est],[min(real(y)) max(real(y))],'r','linewidth',2); 
            legstr={'(scl) window','measured','estimated'};
        end
        grid on; legend(legstr); 
        title('time domain pulse response (real part)'); xlabel('time (s)'); ylabel('amplitude');

        subplot(2,1,2); 
        plot(tplot,real(hw(indplot)),'linewidth',2); grid on;
        title('time domain gated pulse response (real part)'); xlabel('time (s)'); ylabel('amplitude');
	end
end

function w=retrieveWin(Ns,ind_win,win)
    w=zeros(Ns,1);
    Nwin=length(ind_win);  
	if strcmp(win,'blackman')
        w(ind_win)=blackman(Nwin);        
    elseif strcmp(win,'hamming')
        w(ind_win)=hamming(Nwin);
    elseif strcmp(win,'blackmanharris')
        w(ind_win)=blackmanharris(Nwin);
    elseif strcmp(win,'kaiserbessel')
        wi=kaiser(round(Nwin/2),32); % alpha is set to a high value to ensure a stepped slope 
        wkb=zeros(1,Nwin);
        for cont=1:Nwin 
            if cont<=round(Nwin/2)
                wkb(cont)=sqrt(sum(wi(1:cont))/sum(wi));
            else
                wkb(cont)=sqrt(sum(wi(1:Nwin-cont))/sum(wi));
            end
        end
        w(ind_win)=wkb;
    elseif strcmp(win,'rectwin')
        w(ind_win)=rectwin(Nwin);
	end
end

function y = pulseRespfromSpar(Nz,Spar,Nzp,x)
	posSij=[zeros(Nz,1); Spar; zeros(Nzp,1)];	% positive part of the spectrum (direct S parameter plus zero padding at the beginning and at the end)
    negSij=conj(posSij(end:-1:1));              % negative part of the spectrum (real signal is assumed)
    Hij=[posSij; negSij];                       % calculate the transfer function 
    Yij=Hij.*x.freq;                            % calculate output with a known analysis pulse input
    y=ifft(Yij);                                % fast convolution to calculate the time domain pulse response
    %[hij(Ns/2+1:end); hij(1:Ns/2)];            % negative time should be considered properly and causality ensured after delay correction      
end

function H = freqDomTimeGating(Nz,Spar,Nzp,x,w) 
	posSij=[zeros(Nz,1); Spar; zeros(Nzp,1)];	
    negSij=conj(posSij(end:-1:1));              
    Hij=[posSij; negSij];
    Yij=Hij.*x.freq;     
    W=fft(w);
    H=conv(Yij,W)./conv(x.freq,W);              % direct implementation of normalized time gating in the frequency domain
end