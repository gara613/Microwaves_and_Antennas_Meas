% function yhat = mySimpleFit(npartes, x, y, ponderador)
%
% Function to fit a data-set by a linear interpolation based on a direct least 
% squared error formulation, and a sliding window to fit npartes of the vector x
%
% output:
%   - yhat  = vector with the estimated value 
% inputs: 
%   - x     = values of the independent variable
%   - y     = values of the dependent variable (dataset)
%   - vargarin:
%       - npartes    = number of partitions to the x interval
%       - ponderador = sliding window shape, can be 'lineal' (default) or 'triangle' recommended
%
% Germán Augusto Ramírez Arroyave
% Universidad Nacional de Colombia
% CMUN 2014

function [yhat, rmse] = mySimpleFit(x, y, varargin)
    x = x.'; % due to the formulation used
    y = y.';
    if ~isempty(varargin)
        ponderacion = varargin{1};
        Npartes = varargin{2};
    else 
        Npartes = 1; % not a big effort
        ponderacion = 'lineal';
    end
    
    M = length(x);              % Number of samples in the input signal
    k = floor(M/Npartes);       % Number of samples in the interpolation window
    k_2 = floor(k/2);           % in case of triangle weigthing
    yhat = zeros(1,length(x));

    % weigthing vectors, 'triangle' or 'lineal' (default)
    if strcmp(ponderacion,'triangle')
        if mod(k,2)==0 % even values of k
            ponderador(1:k_2) = 1/x(k_2)*x(1:k_2); ponderador(k_2+1:k) = 1-1/x(k/2)*x(1:k_2);
        else % odd values of k
            ponderador(1:k_2+1) = 1/x(k_2+1)*x(1:k_2+1); ponderador(k_2+1:k) = 1-1/x(k_2+1)*x(1:k_2+1);
        end            
    else    % only sqare window so far
        ponderador = 0.5*ones(1,k);
    end
 
    % initialize the lower end
    ind = 1:k_2;
    A = [sum(x(ind)) k_2; x(ind)*x(ind).' sum(x(ind))];
	b = [sum(y(ind)); sum(x(ind).*y(ind))];
	coefs = A\b;
	yhat(ind) = yhat(ind) + (coefs(1)*x(ind) + coefs(2)).*ponderador(k_2+1:end);

    % Linear interpolation algorithm using a direct least squared error formulation
   	ind = 1:k;
    for cont = 1:2*Npartes-1
        A = [sum(x(ind)) k; x(ind)*x(ind).' sum(x(ind))];
        b = [sum(y(ind)); sum(x(ind).*y(ind))];
        coefs = A\b;
        yhat(ind) = yhat(ind) + (coefs(1)*x(ind) + coefs(2)).*ponderador;       %x = (A.^T A)^-1 A^T b
        ind = ind+k_2;
    end
    
    % initialize the upper end
    ind = min(ind):M;
    N = length(ind);
    A = [sum(x(ind)) N; x(ind)*x(ind).' sum(x(ind))];
    b = [sum(y(ind)); sum(x(ind).*y(ind))];
    coefs = A\b;
    yhat(ind) = yhat(ind) + (coefs(1)*x(ind) + coefs(2)).*ponderador(1:N);
        
    ind = k*Npartes:M;
    N = length(ind);
    A = [sum(x(ind)) N; x(ind)*x(ind).' sum(x(ind))];
    b = [sum(y(ind)); sum(x(ind).*y(ind))];
    coefs = A\b;
    yhat(ind) = yhat(ind) + (coefs(1)*x(ind) + coefs(2)).*ponderador(1:N);

    rmse = sqrt(1/length(x)*norm(y-yhat));
    yhat = yhat.';
return 