%% function [pkVals, pkFreqs, pkInds, pkBW, pkQ] = my_pksFinder(X_dB,freq,varargin)
%
% Basic peak finder routine: Searches the dataset (assumed as S21 magnitude in dB) 
% for peaks along frequency sweep and determines the amplitude, frequency of peak location
% as well as the 3 dB bandwidth and the quality factor. 
% 
% - Inputs
%	* X_dB: Real Matrix (N_f x N_exps) with the data to process arranged in
%	columns (dB), (N_f is the number of frequency samples, N_exps is the number of experiments)
%   * freq: Frequency vector (N_f x 1). Uniform sampling assumed but not validated
%   * clipVal: (optional) (dB) Clipping value below the maximum of each vector,
%               Caution: values below threshold are discarded in the processing. 
%    
% - Outputs
%   * pkVals: (N_pks x N_exps): X_dB values at the peaks
%	* pkFreqs: (N_pks x N_exps): Frequencies of the peaks occurrence
%   * pkInds: (N_pks x N_exps): Indexes of the peaks within the frequencies vector
%   * pkBW: (N_pks x N_exps x 3) % 3: f_ini, f_fin, and BW. (Absolute BW).
%       If pkBW = 0, only two samples found in the peak's 3 dB BW and cannot be properly processed
%   * pkQ: (N_pks x N_exps) Calculated as the reciprocal of fractional BW
%       If pkQ = inf, only two samples found in the peak's 3 dB BW and cannot be properly processed
% 
% Ex:
% freq = linspace(1e9,5e9,1001);
% Z01 = 50;
% Z02 = 50;
% R = 1e-3;
% L = 50e-9; 
% C = 56.4e-15;
% Z = R + 1i*2*pi*freq*L + 1./(1i*2*pi*freq*C);
% Ds = Z + Z01 + Z02;
% S11 = (Z01 + Z02 - Z)./Ds;
% S21 = 2*sqrt(Z01*Z02)./Ds;
% S21_dB = 20*log10(abs(S21)).';
% plot(freq,S21_dB,'linewidth',2); grid on,
% [pkVals, pkFreqs, pkInds, pkBW, pkQ] = my_pksFinder(S21_dB,freq,80)
% 
% Ex: 
% freq = linspace(1e9,5e9,1001);
% Z01 = 50;
% Z02 = 50;
% R1A = 1e-3;
% L1A = 40e-9; 
% C1A = 106.4e-15;
% R1B = 1e-3;
% L1B = 60e-9; 
% C1B = 36.4e-15;
% R2 = 1e-3;
% L2 = 50e-9; 
% C2 = 56.4e-15;
% Z1A = R1A + 1i*2*pi*freq*L1A + 1./(1i*2*pi*freq*C1A);
% Z1B = R1B + 1i*2*pi*freq*L1B + 1./(1i*2*pi*freq*C1B);
% Z2A = R2 + 1i*2*pi*freq*L2 + 1./(1i*2*pi*freq*C2);
% Z1 = Z1A.*Z1B./(Z1A+Z1B);  
% Z2 = Z2A; 
% Ds_1 = Z1 + Z01 + Z02;
% Ds_2 = Z2 + Z01 + Z02;
% S21_1 = 2*sqrt(Z01*Z02)./Ds_1;
% S21_2 = 2*sqrt(Z01*Z02)./Ds_2;
% S21_dB(:,1) = 20*log10(abs(S21_1)).';
% S21_dB(:,2) = 20*log10(abs(S21_2)).';
% [pkVals, pkFreqs, pkInds, pkBW, pkQ] = my_pksFinder(S21_dB,freq,80)
% inds2plot = [1,2];
% Nplots = length(inds2plot);
% pkVals3dB = pkVals(:,inds2plot);
% pkVals3dB = [pkVals3dB(:),pkVals3dB(:)]-3;
% xx = pkBW(:,inds2plot,1:2);
% bw2plot = reshape([xx(:,:,1),xx(:,:,2)],Nplots*size(pkInds,1),2);
% figure,
% plot(freq, S21_dB(:,inds2plot),'linewidth',2), grid on; 
% hold on
% plot(pkFreqs(:,inds2plot), pkVals(:,inds2plot),'r*');
% plot(bw2plot',pkVals3dB','k','linewidth',2);
% xlabel('frequency (GHz)'); ylabel('S_{21} (dB)');
% axis([1e9,5e9,-30,0])
% 
% Germán Ramírez, 
% EPFL - MAG, July 2024

function [pkVals, pkFreqs, pkInds, pkBW, pkQ] = my_pksFinder(X_dB,freq,varargin)
    if ~isempty(varargin)
        clipVal = varargin{1};
    else 
        clipVal = -inf; 
    end
    freq = freq(:);

    [N_fpts, N_exps] = size(X_dB); 
    indVec = 1:N_fpts;
    delta_f = mean(unique(diff(freq)));  
    if length(unique(diff(freq))) > 1 % abs(delta_f/min(unique(diff(freq)))) > 1.01
        warning('Nonuniform frequency sampling detected, results may be inaccurate for BW and Q');
    end

    maxVal = max(X_dB,[],1); 
    X_dB(X_dB<maxVal-clipVal)=-inf;

    % Derivative of measured data. (Can be really noisy, pre-processing should be required)
    X_dB_der = [zeros(1,N_exps);diff(X_dB,1,1)/delta_f];    % Zero padding at the beginning
    % Peaks are found as the change of slope sign 
    pksFound = [diff(sign(X_dB_der));zeros(1,N_exps)] < -1;	% Zero padding at the end
    N_pks = sum(pksFound);                                  % Number of peaks for each experiment
    pkVals = zeros(max(N_pks),N_exps);
    pkFreqs = zeros(max(N_pks),N_exps);
    pkInds = zeros(max(N_pks),N_exps);
	pkBW = zeros(max(N_pks),N_exps,3);
    pkQ= zeros(max(N_pks),N_exps);
        
    for cont = 1:N_exps                 % Process the peaks for each of the signals provided 
        pkVals(1:N_pks(cont),cont) = X_dB(pksFound(:,cont),cont);
        pkFreqs(1:N_pks(cont),cont) = freq(pksFound(:,cont));
        pkInds(1:N_pks(cont),cont) = indVec(pksFound(:,cont));
        for conpk = 1:N_pks(cont)       % Process each peak for the current signal 
            xx = freq(X_dB(:,cont)>=pkVals(conpk,cont)-3);              % Find the samples 3dB below the peak
            [~,yy] = min(abs( pkFreqs(conpk,cont) - xx ));              % Index of the peak in the 3dB BW
            aa = [1; find(diff(xx) > 1.01*mean(diff(xx))); length(xx)]; % Check for peaks in discontinuous frequency ranges
            if length(xx) > 2             
                fLim = sum(aa<=yy); 
                if aa(fLim)+fLim-1 <= length(xx)
%                pkBW(conpk,cont,:) = deal([xx(aa(fLim)+1), xx(aa(fLim+1)), xx(aa(fLim+1)) - xx(aa(fLim)+1)]);                    
                    pkBW(conpk,cont,:) = deal([xx(aa(fLim)+fLim-1), xx(aa(fLim+1)), xx(aa(fLim+1)) - xx(aa(fLim)+fLim-1)]);
                    pkQ(conpk,cont) = pkFreqs(conpk,cont)/pkBW(conpk,cont,3);
                else
                    pkBW(conpk,cont,:) = 0; % This case needs to be dealt with. Detected when there are nonuniform samples within the 3dB BW of a single peak
                    pkQ(conpk,cont) = inf;
                end
            else                            % Peak found, but it's too narrow to characterize with the current signal sampling
                pkBW(conpk,cont,:) = 0;     % No interpolation is done, just "anormal" values returned
                pkQ(conpk,cont) = inf;  
            end
        end
    end
end