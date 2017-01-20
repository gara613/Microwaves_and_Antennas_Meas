% Example of the TRL calibration method procedure
 
%% Configuration parameters
Order = 8;% Max order of the rational approximation used
cleaning = true;
 
%% Measurements reading sub-section
% S parameter matrix for Thru standard will be called T, R for Reflect standard, and L for Line
l = 2.22e-2;   % line standard length in meters 
[T, freqst, R_T] = readSpars('Thru_30MHz_3GHz.s2p');    
[R, freqsr, R_R] = readSpars('Reflect_30MHz_3GHz.s2p'); 
[L, freqsl, R_L] = readSpars('Line_30MHz_3GHz.s2p');    
[Smeas, freqsm, R_M] = readSpars('TestLine75Ohm.s2p');  
 
% validate measurements before calling the correction routine
if (freqst ~= freqsr) | (freqst ~= freqsl) | (freqst ~= freqsm)
    disp('Las frecuecias de las mediciones no coinciden');
end
if (R_T ~= R_R) | (R_R ~= R_L) | (R_L ~= R_M) | ~all([R_L,R_T,R_R,R_M])
    disp('La impedancia característica de las mediciones no coincide');
else
    Z0 = R_M;    % Characteristic impedance of the line
end
 
%% Data cleaning
% Rational model is a physical model approximation to the S-parameters,
% this approximation is prone to error for SNR below 23dB
if cleaning
    T = cleanSpars(T, freqst, 'RationalS', 'canonical', Order); 
    R = cleanSpars(R, freqsr, 'RationalS', 'canonical', Order); 
    L = cleanSpars(L, freqsl, 'RationalS', 'canonical', Order); 
    Smeas = cleanSpars(Smeas, freqsm, 'RationalS', 'canonical', Order); 
end 
  
%% Call to the data correction routine
% Even if only five measurements are needed (and used) all twelve are passed to the correctiion routine
Sdut = correctionTRL(T,R,L,Smeas,freqsm);
 
%% Plot results
plotSpars(freqsm, {Smeas, Sdut}, 'fase', 'cmp','','Uncoupled line', {'Meas', 'TRL Corr'});