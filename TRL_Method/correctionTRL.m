% TRL calibration method procedure, from three measurements of S parameters
% for Thru, Reflect and Line networks: [T], [R], [L]. Reference: Pozar.
%
% sDUT = correctionTRL(T,R,L,Smeas,freqsm, varargin)
% varargin should specify: 
%   option: Method of choice for the selection of the root in the e^{gamma l} equation (Slope criterion is taken by default)
%       - minJump
%       - sectSlope
%   Z0: characteristic impedance of the measurements (50 Ohm is taken by default)  
%   l:  length of the Line standard (1cm is taken by default)
%
% Universidad Nacional de Colombia
% Germán Augusto Ramírez Arroyave
% Grupo de Investigación en Telecomunicaciones - CMUN (2014)

function sDUT = correctionTRL(T, R, L, Smeas, freqsm, varargin)

    % Default values for the optional parameters involved
    option = 'sectSlope'; 
    Z0 = 50;    % line's characteristic impedance
    l = 1e-2;   % line length in meters, maybe some formulation should make use of this value!?
    
    if ~isempty(varargin)
        switch length(varargin)
            case 1
                option = varargin{1};     
            case 2
                option = varargin{1};     
                Z0 = varargin{2};     
            case 3
                option  = varargin{1};     
                Z0 = varargin{2};     
                l = varargin{3};
            otherwise
                disp('Wrong number of optional input parameters');
        end
    end

    % Data cleaning and frequency sampling validation are assumed and should be performed before calling this routine 
    T11 = squeeze(T(1,1,:));
    T12 = squeeze(T(1,2,:)); 
    R11 = squeeze(R(1,1,:));
    L11 = squeeze(L(1,1,:));
    L12 = squeeze(L(1,2,:));
    ABCDmeas = myS2ABCD(Smeas,Z0);

    %% Determine the correction factors
    % Propagation coefficient of the line
    expgammal_mas   = ( L12.^2 + T12.^2 - (T11 - L11).^2 + sqrt( (L12.^2 + T12.^2 - (T11 - L11).^2).^2 - 4*L12.^2.*T12.^2 ) ) ./ (2*L12.*T12);
    expgammal_menos = ( L12.^2 + T12.^2 - (T11 - L11).^2 - sqrt( (L12.^2 + T12.^2 - (T11 - L11).^2).^2 - 4*L12.^2.*T12.^2 ) ) ./ (2*L12.*T12);
              
    %% Different algorithms to choose the right root for e^{\gamma l}
	% Minimum jump from one value in the root to the next
    if strcmp(option,'minJump') 
    	howFar = 1;           
        if sign(angle(expgammal_mas(1+howFar)) - angle(expgammal_mas(1)))>=0   % Choose as starting root the one with increasing phase
            expgammal = expgammal_mas;
            for cont = 1:length(expgammal)-howFar % Minimum jump from initial root to create a 'continuous' root
                if abs(expgammal(cont) - expgammal_mas(cont+howFar)) > abs(expgammal(cont) - expgammal_menos(cont+howFar))
                    expgammal(cont+howFar) = expgammal_menos(cont+howFar);
                % In case of equality use a slope criterion (curves can cross each other in frequencies where no sample exists and hence cross is not detected by this condition)
                elseif abs(expgammal(cont) - expgammal_mas(cont+howFar)) == abs(expgammal(cont) - expgammal_menos(cont+howFar))
                    trend = sign(angle(expgammal(cont)) - angle(expgammal(cont-howFar))); % assumes that there are no crosses in the first sample
                    if sign(angle(expgammal_menos(cont+howFar)) - angle(expgammal_menos(cont))) == trend
                        expgammal(cont+howFar) = expgammal_menos(cont+howFar); 
                    end
                end
            end 
        else % Choose as starting root the one with increasing phase
            expgammal = expgammal_menos;  
            for cont = 1:length(expgammal)-howFar 
                if abs(expgammal(cont) - expgammal_menos(cont+howFar)) > abs(expgammal(cont) - expgammal_mas(cont+howFar))
                    expgammal(cont+howFar) = expgammal_mas(cont+howFar);
                elseif abs(expgammal(cont) - expgammal_menos(cont+howFar)) == abs(expgammal(cont) - expgammal_mas(cont+howFar))
                    trend = sign(angle(expgammal(cont)) - angle(expgammal(cont-howFar)));
                    if sign(angle(expgammal_mas(cont+howFar)) - angle(expgammal_mas(cont))) == trend
                        expgammal(cont+howFar) = expgammal_mas(cont+howFar); 
                    end
                end
            end
       end
    % creates a root with a continuous positive phase
    elseif strcmp(option,'sectSlope') % Can have problems due to sharp (noisy) changes in the curves if those are not previously fitted
        ang_g_menos = (angle(expgammal_menos));
        ang_g_mas = (angle(expgammal_mas)); 
        tol = 1e6*eps;                % Already tried with a for loop and with the 'diff' command ang got similar results...
        expgammal = expgammal_menos.*( ([ang_g_menos(1); ang_g_menos(2:end)-ang_g_menos(1:end-1)]) > tol ) +...
                     expgammal_mas.* ( ([ang_g_mas(1);   ang_g_mas(2:end)-ang_g_mas(1:end-1)]) > tol ); % both roots shouldn't have positive phase at the same time
    end
     
%% Debug: These plots must be erased once I figure out how to choose the correct root!
    buggy = false;
    if buggy
        figure; plot(freqsm, abs(expgammal), freqsm, unwrap(angle(expgammal)*180/pi));
        legend('abs e^{\gamma l}', 'angle e^{\gamma l}'); title('e^{\gamma l} tomado'); 
        figure; plot(freqsm, real(log(expgammal)/l), freqsm, imag(log(expgammal)/l));
        legend('\alpha', '\beta'); title('\gamma');
    end
    
    %% S parameters of the error boxes (Symmetry is assumed for S12 and S21)
    S22 = (T11 - L11) ./ (T12 - L12./expgammal);
    S11 = T11 - S22.*T12;
    S12 = sqrt(T12.*(1 - S22.^2));
    % Reflection coefficient of the reflect standard
    GammaL = (R11 - S11) ./ (S12.^2 + S22.*(R11 - S11)); % this may be (un)used...
      
    % Error boxes are assumed identical
    S(1,1,:) = deal(S11); S(1,2,:) = deal(S12);
    S(2,1,:) = deal(S12); S(2,2,:) = deal(S22);
       
    %% Retrieve the ABCD parameters of error boxes
    ABCDerror = myS2ABCD(S, Z0);
    
    % ABCD matrix for the second error box due to connection inversion
    DBCAerror(1,1,:) = ABCDerror(2,2,:); DBCAerror(1,2,:) = ABCDerror(1,2,:);
    DBCAerror(2,1,:) = ABCDerror(2,1,:); DBCAerror(2,2,:) = ABCDerror(1,1,:);

    %% Find the corrected S parameters
    ABCDdut = zeros(size(ABCDmeas));
     
    for cont = 1:length(freqsm)
        ABCDdut(:,:,cont) = ABCDerror(:,:,cont) \ ABCDmeas(:,:,cont) / DBCAerror(:,:,cont);
    end
    
    sDUT = myABCD2S(ABCDdut, Z0);