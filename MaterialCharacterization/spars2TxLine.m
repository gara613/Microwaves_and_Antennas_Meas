%% function [x,y] = spars2TxLine(Spars)
% Obtain the transmission line propagation constant and characteristic impedance from the port S-parameters 
% 
% Inputs: 
%   * Spars: 2x2xN_freqs, S parameters matrix for a two port network
% Outputs:
%   * x = exp(-\gamma*l), where l is the length of the line, not known to this function
%   * y = Z_Line/Z_ref, normalized line impedance, Z_ref is the reference impedance
%
% Ex: 
% fre = linspace(1,10,101)'*1e9;
% Zref = 50; 
% c = 3e8;
% lambda = c./fre;
% l = 5e-2;             %
% ZL = 53; YL = 1/ZL;
% beta = (1-1i*0.005)*2*pi./lambda;
% ABCD(1,1,:) = cos(beta.*l);       % A
% ABCD(1,2,:) = 1i*ZL*sin(beta*l);  % B
% ABCD(2,1,:) = 1i*YL*sin(beta*l);  % C
% ABCD(2,2,:) = cos(beta*l);        % D
% S = myABCD2S(ABCD,Zref);
% [x,y] = spars2TxLine(S);
% figure,
% plot(fre/1e9,1/l*real(log(x.minus)), fre/1e9,1/l*imag(log(x.plus)), 'linewidth',2); grid on
% hold on, 
% plot(fre/1e9,real(beta), '--', fre/1e9,imag(beta), '--', 'linewidth',2); 
% xlabel('Frequency (GHz)'); ylabel('\gamma'); 
% legend('Re_{est}','Im_{est}','Re_{act}','Im_{act}')
% figure,
% plot(fre/1e9,real(y*Zref), fre/1e9,ZL*ones(size(fre)),'--', 'linewidth',2); grid on
% xlabel('Frequency (GHz)'); ylabel('Z_0 (\Omega)'); 
% legend('Estimated','Actual');
% 
% Germán A. Ramírez, 
% MAG, EPFL, Dic 2023

function [x,y] = spars2TxLine(Spars,varargin)
    option = 'symm';
    if exist('varargin','var') & ~isempty(varargin)
        option = varargin{1};
    end

    if strcmpi(option,'symm')
        S11 = squeeze(Spars(1,1,:));
        S21 = squeeze(Spars(2,1,:));        

        x.plus =  (1-S11.^2+S21.^2 + sqrt((S11.^2-S21.^2-1).^2-4*S21.^2))./(2*S21);
        x.minus = (1-S11.^2+S21.^2 - sqrt((S11.^2-S21.^2-1).^2-4*S21.^2))./(2*S21);   
        y = sqrt( ((1+S11).^2-S21.^2)./((1-S11).^2-S21.^2) );
    else 
        S11 = squeeze(Spars(1,1,:));
        S21 = squeeze(Spars(2,1,:));
        S12 = squeeze(Spars(1,2,:));
        S22 = squeeze(Spars(2,2,:));
        
        A = (1-S11.^2).*(1-S22.^2) + 2*S21.*S12.*(1-S22.*S11) + S21.^2.*S12.^2;
        B = (1-S11.^2).*(1-S22.^2) - 2*S21.*S12.*(1+S22.*S11) + S21.^2.*S12.^2;
        x.plus = ( sqrt(A)+sqrt(B) )./ (2*S21);
        x.minus = ( sqrt(A)-sqrt(B) )./ (2*S21);
        y = sqrt( ((1+S11).*(1+S22)-S21.*S12) ./ ((1-S11).*(1-S22)-S21.*S12) );
    end
end