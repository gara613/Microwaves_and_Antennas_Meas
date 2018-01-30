% Function that returns the effective transmission line parameters from its
% S pars measurements
% inputs:
%   - Smat: S parameters matrix (as returned by readSpars) 
%   - freqs: frequency sample vector
%   - mur: effective relative permeability of TEM line
%   - Zref: port's reference impedance 
%   - varargin: with fields 
%       * dl, epsr_pt, mur_pt (used for de-embedding port length)
% outputs:
%   - effPars: structure with fields
%       * Z, gamma, epsreff, tandeltaeff, R,L,G,C
%
% Germán Augusto Ramírez - CMUN 2017

function effPars=txLineParsfromS(Smat,freqsr,mur,l,Zref,varargin)
    c=3e8;
    dl=0; epsr_pt=1; mur_pt=1;
    if ~isempty(varargin)
        dl=varargin{1};         %de-embedding length    
        epsr_pt=varargin{2};    %port material relative permittivity
        mur_pt=varargin{3};     %port material
    end
    
    % Symmetry is assumed for S parameters, hence only S11 and S21 are used
    S11=squeeze(Smat(1,1,:));
    S21=squeeze(Smat(2,1,:));
    % De-embedding of port length
    S11=S11.*exp(-1i*2*2*pi*freqsr/c*sqrt(mur_pt*epsr_pt)*dl);
    S21=S21.*exp(-1i*2*2*pi*freqsr/c*sqrt(mur_pt*epsr_pt)*dl);      % assumes connector symmetry

    A = sqrt(-(S11.^2+2*S11-S21.^2+1)./(-S11.^2+2*S11+S21.^2-1));
    effPars.gamma{1} = -1/l*log(S21.*(A/2+0.5)-(S11-1).*(S11+A.*S11-A+1)./(2*S21));                 % JA derivation from S pars measurements

    effPars.Z = Zref.*sqrt(((1+S11).^2-S21.^2)./((1-S11).^2-S21.^2));
    effPars.gamma{2} =  1/l*acosh((1-S11.^2+S21.^2)./(2*S21));                                      % from S pars to ABCD conversion (naive)
    effPars.gamma{3} = -1/l*log((1-S11.^2+S21.^2-sqrt((1+S11.^2-S21.^2).^2-4*S11.^2))./(2*S21));	% from S pars to ABCD conversion (stable)

    effPars.epsreff=(unwrap(imag(effPars.gamma{1}))*c./(2*pi*freqsr)).^2/mur;
    
    effPars.R = real(effPars.gamma{1}.*Zref);
    effPars.L = imag(effPars.gamma{1}.*Zref)./(2*pi*freqsr);
    effPars.G = real(effPars.gamma{1}./Zref);
    effPars.C = imag(effPars.gamma{1}./Zref)./(2*pi*freqsr);
    effPars.tandelta = effPars.G./(2*pi*freqsr.*effPars.C);
end