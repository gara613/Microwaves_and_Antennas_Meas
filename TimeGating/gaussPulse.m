% function x=gaussPulse(t,fc,bw,Nfft,option)
% creates a (un)modulated gaussian pulse, and returns a vector with the time and frequency domain values. 
% Inputs:
%   - t: time support of the intended signal 
%   - fc: central frequency of the signal
%   - bw: bandwidth of the signal (inverse of the std deviation)
%   - option: 'normal' or 'modulated'
% Outputs: 
%   - x: structure with fields: time and freq with the signal vector in time and frequency domain
% Germán Augusto Ramírez Arroyave
% CMUN - Universidad Nacional de Colombia 2016

function x=gaussPulse(t,fc,bw,Nfft,option)
    if strcmp(option,'impulse');
        x.time=zeros(length(t),1); 
        [~,ind]=min(abs(t));
        x.time(ind)=1;
    elseif strcmp(option,'normal');
        x.time=exp(-0.5*bw*bw*t.^2)'; %bw/sqrt(2*pi)*
    elseif strcmp(option,'modulated');
        x.time=(exp(-0.5*bw*bw*t.^2).*sin(2*pi*fc*t))';
	elseif strcmp(option,'sinc');
        x.time=sinc(bw*t)';
    end
	x.freq=fft(x.time,Nfft);
end