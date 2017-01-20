% TRL calibration method procedure, from three measurements of S parameters
% for Thru, Reflect and Line standards [T], [R], [L] matrixes notation used
% as in Pozar. i.e. Lij equals Sij for the Line network.
%
% Universidad Nacional de Colombia
% Germán Augusto Ramírez Arroyave
% Grupo de Investigación en Telecomunicaciones - CMUN (2014)

clc, clear, %close all % Remove this after validation of results!

%% Physical parameters involved
Order = 8;  % This is for the order of the rational approximation or the bandwidth of the Gaussian kernel, this should be chosen automatically

%% Configuration
validationData = 'VNA_fisica'; % CST_circuital, 'CST_FullWave', 'VNA_EMC', 'VNA_fisica'
criterio = 'sectSlope'; % 'minJump', 'sectSlope'
cleaning = true;
noise = true;     % use true only when simulation data is used in order to verify noise immunity of calibration
SNR = 30; % dB 
show_interp = false;
cmp_gaussKer = false;
reference = true;

%% Measurements reading sub-section
% S parameter matrix for Thru standard will be called T, R for Reflect standard, and L for Line

switch validationData
    case 'CST_circuital'
        l = 3e-2;   % line length in meters (FWEM & circuit simulation)
        [T, freqst, R_T] = readSpars('Thru_ckt_wError.s2p');	% Thru is done by a pair of line segments of 9.53mm of 49 Ohm and eps_r=2.1 substrate
        [R, freqsr, R_R] = readSpars('Reflect_ckt_wError.s2p');	% Reflect is done by a pair of open lines of 2cm
        [L, freqsl, R_L] = readSpars('Line_ckt_wError.s2p');	% Line is made using 3cm long, 50 Ohm line, using an eps_r=3 substrate, and the same pair of line segments used in the thru as connectors
        [Smeas, freqsm, R_M] = readSpars('TestLine_75Ohm_ckt_wErrort.s2p'); % Test_ckt75_wError_Design 

    case 'CST_FullWave'
        l = 2.22e-2;   % line standard length in meters 
        [T, freqst, R_T] = readSpars('Thru_30MHz_3GHz.s2p');	
        [R, freqsr, R_R] = readSpars('Reflect_30MHz_3GHz.s2p');	
        [L, freqsl, R_L] = readSpars('Line_30MHz_3GHz.s2p');    
        [Smeas, freqsm, R_M] = readSpars('TestLine75Ohm.s2p');  

    case 'VNA_EMC'
        l = 5e-2;   % line length in meters (Measurement)
        [T, freqst, R_T] = readSpars('thru.s2p');
        [R, freqsr, R_R] = readSpars('reflect.s2p');
        [L, freqsl, R_L] = readSpars('line.s2p');
        [Smeas, freqsm, R_M] = readSpars('line.s2p'); % S parameter matrix for the measurement to be corrected
    case 'VNA_fisica'
        l = 5e-2;   % line length in meters (Measurement)
        [T, freqst, R_T] = readSpars('Thru.txt');
        [R, freqsr, R_R] = readSpars('Reflect.txt');
        [L, freqsl, R_L] = readSpars('Line.txt');
        [Smeas, freqsm, R_M] = readSpars('Line.txt'); % S parameter matrix for the measurement to be corrected
end
Smeas_orig = Smeas; % Enable if want to compare simulation result with synthetic noise addition, rational, and gaussian kernel interpolation

% validate measurements before calling the correction routine
if (freqst ~= freqsr) | (freqst ~= freqsl) | (freqst ~= freqsm)
    disp('Las frecuecias de las mediciones no coinciden');
end
if (R_T ~= R_R) | (R_R ~= R_L) | (R_L ~= R_M)
    disp('La impedancia característica de las mediciones no coincide');
else
    Z0 = R_M;    % Characteristic impedance of the line
end

%% Synthetic noise addition (Only for test "measurements" with simulation data)
if noise
    snr = 10^(SNR/10); % lineal
     
    pS = powerSpars(T);
    an = sqrt(pS/(2*snr*size(T,3)));
    N = an*randn(size(T)) + 1i*an*randn(size(T));
    T = T + N; 
    
    pS = powerSpars(R);
    an = sqrt(pS/(2*snr*size(R,3)));
    N = an*randn(size(R))+ 1i*an*randn(size(R));
    R = R + N;
    
    pS = powerSpars(L);
    an = sqrt(pS/(2*snr*size(L,3)));
    N = an*randn(size(L)) + 1i*an*randn(size(L));
    L = L + N;
    
    pS = powerSpars(Smeas);
    an = sqrt(pS/(2*snr*size(Smeas,3)));
    N = an*randn(size(Smeas))+ 1i*an*randn(size(Smeas));
    Smeas = Smeas + N;
end 
Smeas_noi = Smeas; % Enable if want to compare simulation result with synthetic noise addition, rational, and gaussian kernel interpolation

%% Data cleaning
% Rational model is a physical model approximation to the S-parameters,
% this approximation is prone to error for SNR below 23dB
if cleaning
    T = cleanSpars(T, freqst, 'RationalS', 'canonical', Order); 
    R = cleanSpars(R, freqsr, 'RationalS', 'canonical', Order); 
    L = cleanSpars(L, freqsl, 'RationalS', 'canonical', Order); 
    Smeas = cleanSpars(Smeas, freqsm, 'RationalS', 'canonical', Order); 
    Smeas_rat = Smeas; % Enable if you want to compare simulation result with synthetic noise addition, rational, and gaussian kernel interpolation
end 

if show_interp % Enable if want to compare measurement with rational interpolation
    figure, plotSpars(freqsm, {Smeas_orig, Smeas_noi, Smeas}, 'fase', 'cmp', {'Original', ['plus SNR=' num2str(SNR) 'dB Noise'], 'Rational'});
end

% This method is very effective for very low SNR measurements, altough discouraged in normal situations where SNR is high 
if cmp_gaussKer
    Tgk = cleanSpars(T, freqst, 'KernelS', 'gaussi', 0.5, Order); 
    Rgk = cleanSpars(R, freqsr, 'KernelS', 'gaussi', 0.5, Order); 
    Lgk = cleanSpars(L, freqsl, 'KernelS', 'gaussi', 0.5, Order); 
    Smeas_gk = cleanSpars(Smeas_noi, freqsm, 'KernelS', 'gaussi', 0.5, Order); 
    %Smeas=Smeas_gk; 

    % Enable if want to compare measurement with synthetic noise addition, rational, and gaussian kernel interpolation
    figure, plotSpars(freqsm, {Smeas_orig, Smeas_noi, Smeas_rat, Smeas_gk}, 'fase', 'cmp', {'Original', 'Ruido sint', 'Rational', 'G Kernel'});
end

%% Call to the data correction routine
% Even if only five measurements are needed (and used) all twelve are passed to the correctiion routine (all of them should have the same frequency sampling)
Sdut = correctionTRL(T,R,L,Smeas,freqsm, criterio, Z0, l);

%% Plot results
figure, plotSpars(freqsm, {Smeas, Sdut}, 'fase', 'cmp', {'Meas', 'TRL Corr'});
 
% If there is a reference for comparison, load it here (freqsex must be equal to freqsm!)
if reference
    switch validationData
        case 'CST_circuital'
            [Sexpected, freqsex] = readSpars('TestLine_75Ohm_ckt.s2p'); % Ideal_ckt_Line75
        case 'CST_FullWave'
            [Sexpected, freqsex] = readSpars('TestLine75OhmIdeal.s2p'); 
    end
	figure, plotSpars(freqsm, {Smeas, Sdut, Sexpected}, 'fase', 'cmp', {'Meas', 'TRL Corr', 'Expected'});
end