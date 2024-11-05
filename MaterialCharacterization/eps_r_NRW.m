%% function [eps_est,mu_est] = eps_r_NRW(S11,S21,fre, MUT_L)
% Following the equations in the papers from Nicolson, Ross, and Weir 
% for determining the permittivity of a Material Under Test filling the
% transversal section of a Waveguide
%
% Inputs: 
%   - S11: N_meas x N_freqs (complex)
%   - S21: N_meas x N_freqs (complex)
%   - fre: N_freqs x 1
%   - MUT_L:
%   - f_cut:
%   - method (optional). Use of direct NRW method by default. 
%       Iterative solution via Baker-Jarvis procedure is proposed as an extension but not currently implemented. 
% Outputs: 
%   - eps_est: N_meas x N_freqs (complex)
%   - mu_est: N_meas x N_freqs (complex)
%
% Ex: R&S Measurement of Dielectric Material Properties Application Note
% S11 = 0.856*exp(1i*163.2*pi/180);     (eps_r = 5.7 + j7.2)
% S21 = 0.609*exp(-1i*140.5*pi/180);
% sL = 0.4e-2;
% f_0 = 8e9;
% f_c = 5.26e9;
% [eps_est,mu_est] = eps_r_NRW(S11,S21,f_0,sL,f_c,'NRW') 
% [eps_est,mu_est] = eps_r_NRW(S11,S21,f_0,sL,f_c,'NonIter') 
%
% Ex: 
% See: 'eps_r_NRW_test' script
%
% Germán A. Ramírez
% EPFL - MAG, July 2023

function [eps_est,mu_est] = eps_r_NRW(S11,S21,fre,MUT_L,f_cut,varargin)
    c = 3e8;
    if exist('varargin','var') & ~isempty(varargin)
        method = varargin{1};
    else
        method = 'NRW';
    end
 
    fre = fre(:).';
    lambda_0 = c./fre;
    lambda_c = c/f_cut;
    k_0 = 2*pi./lambda_0; 
    k_c = 2*pi/lambda_c; 
    beta_g = sqrt(k_0.^2-k_c.^2);
 
    K = (S11.^2-S21.^2+1)./(2*S11);
    Gamma_plus =  K + sqrt(K.^2-1);
    Gamma_minus = K - sqrt(K.^2-1);
    Gamma = Gamma_plus.*(abs(Gamma_plus)<=1) + Gamma_minus.*(abs(Gamma_minus)<=1);

    P = (S11+S21-Gamma)./(1-(S11+S21).*Gamma);  
    phi = log(1./P); 

	% m = 0*ones(1,length(fre));
    % phi = phi + 1i*m*2*pi;

    beta_est = -1i*phi/MUT_L;
    mu_est = bsxfun(@rdivide,beta_est,beta_g).*(1+Gamma)./(1-Gamma);
    eps_eff = bsxfun(@rdivide,beta_est,beta_g).*(1-Gamma)./(1+Gamma);
    if strcmpi(method,'NRW')
        eps_est = (beta_est.^2 + k_c.^2) ./ (k_0.^2 .*mu_est);
    elseif strcmpi(method,'NonIter')
        if abs(mean(mu_est)-1) <= 1e-2
            eps_est = (beta_est./beta_g).^2;
        else
            eps_est = (1-lambda_0.^2./lambda_c.^2).*eps_eff + lambda_0.^2./lambda_c.^2.*1./mu_est;
        end        
    elseif strcmpi(method,'Iter') % -> Baker-Jarvis paper
        error('Not implemented yet...')
    end
end