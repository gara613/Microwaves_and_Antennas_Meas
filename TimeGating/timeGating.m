% function Smod=timeGating(h,tini,tfin,win)
% Returns the time gated version of an input S parameters matrix, given a 
% cutoff time window
% Inputs:
%   - h: pulse response from S parameters (as calculated by tdSpars)
%   - tini: initial time for truncation
%   - tfin: final time for truncation
%   - win: string to indicate the kind of window, suppported options include ('gaussian','hamming','blackman')
% Outputs:
%   - Smod: Modified S parameters matrix
%
% Currently this processing is carried out by tdSpars.m
% Germán Augusto Ramírez Arroyave
% CMUN - Universidad Nacional de Colombia 2016

function Smod=timeGating(tdS,tini,tfin,win)
    t=tdS.t;
    ind=(t>tini & t<tfin);
    n=t(ind);
    N=length(n);
    
    if strcmp(win,'blackman')
        w=blackman(N);
    elseif strcmp(win,'hamming')
        w=hamming(N);
    elseif strcmp(win,'blackmanharris')
        w=blackmanharris(N);
    elseif strcmp(win,'rectwin')
        w=rectwin(N);
    end

	figure, plot(t,abs(squeeze(tdS.h(1,1,:))),'linewidth',2); 
    hold on; plot(n,w,'r','linewidth',2); legend('measured','window');   
    title('time domain impulse response reflection'); xlabel('time (s)'); ylabel('amplitude');

    figure, plot(t,abs(squeeze(tdS.h(2,1,:))),'linewidth',2); 
	hold on; plot(n,w,'r','linewidth',2); legend('measured','window'); 
    title('time domain impulse response transmission'); xlabel('time (s)'); ylabel('amplitude'); 
    
    for cont1=1:size(tdS.h(1))
        for cont2=1:size(tdS.h(1))
            H_tr=fft(squeeze(tdS.h(cont1,cont2,ind)).*w);
            Smod(cont1,cont2,:)=H_tr;%./x.f;
        end
    end  
end