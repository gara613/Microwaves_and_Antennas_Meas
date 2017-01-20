% function [yhat, rmse] = kernelSpar(x, y, option)
% Function to fit a data-set using kernel interpolation
%
% outputs: estimated vector, rmse
% inputs : 
%   - x      = vector (nx1) with the values of the independent variable
%   - y      = vector (nx1) with the values of the dependent variable (dataset)
%   - option = {'linReg', 'linear', 'polynom', 'gaussi'} 
%   - Suboptions for each kernel kind:
%       - linear: lambda=regularization parameter to control complexity of the model
%       - polynom: lambda, d=degree of the polynomial
%       - gaussi: lambda, gammac=bandwidth of the RBFs (gammac = 1/2simga^2)
%
% Germán Augusto Ramírez Arroyave
% Universidad Nacional de Colombia
% CMUN 2015

function [yhat, rmse] = kernelSpar(x, y, varargin)
    % In case inputs are not vectors
	magOrderX = floor(log10(max(x)));   % Another kind of normalization  
    if size(x,1)==1 && size(y,1)==1
        x = x.'/10^(magOrderX);
        y = y.';
    else
        x = x/10^(magOrderX);          % Removes the exponent
    end
    
    noptions = size(varargin, 2);
    if ~isempty(varargin)
        option = varargin{1};
    else
        option = 'linReg';
    end
    
    %% Fit the dataset using standard linear regression
    if strcmp(option,'linReg')
        w = (x.'*x)\x.'*y;          % weights matrix
        yhat = w'*x;                % estimation of data

    %% Fit the dataset using the regularized kernel expression (Ridge regression)
    elseif strcmp(option,'linear')
        if noptions>=2
            lambda = varargin{2};   % controls the complexity of the model
        else
            lambda = 0.01;          
        end       
        G = x*x.';              
        alpha = (G + lambda*eye(size(x,1)))\y;
        yhat = alpha.'*G;           % estimation

    %% Fit the dataset using nonlinear (polynomial) kernel
    elseif strcmp(option,'polynom')
        if  noptions>=2
            lambda = varargin{2};   % controls the complexity of the model
        else
            lambda = 0.5;           % controls the complexity of the model
        end
        if noptions>=3
            d = varargin{3};   
        else
            d = 7;                  % degree of the model
        end
        c = 1;                      % This parameter usually is fixed to one
        % Fit the dataset using the regularized kernel expression... Ridge regression, poly kernel may suffer from instability        
        G = x*x.';
        K = (G + c).^d; 
        alpha = (K + lambda*eye(size(x,1)))\y;
        yhat = alpha.'*K;

    %% Fit the dataset using Gaussian kernel K(x,z) = exp(-gamma*norm(x-z)^2)
    elseif strcmp(option,'gaussi')
        if  noptions>=2
            lambda = varargin{2};   % controls the complexity of the model
        else
            lambda = 0.05;          
        end
        if  noptions>=3
            gammac = varargin{3};   % a.k.a as the bandwidth of the RBFs \gamma = 1/2simga^2 
        else
            gammac = 0.1; 
        end
%       Straight (naive) implementation of Gaussian Kernel as interpreted by GAR
%        K = zeros(length(x),length(x));
%         for cont1 = 1:length(x) % This is a really time-expensive routine, try to avoid its use until proper fast implementation
%             for cont2 = 1:length(x)
%                 K(cont1,cont2) = exp( -gammac*((x(cont1,:)-x(cont2,:))'*(x(cont1,:)-x(cont2,:))));
%             end
%         end  
        %	Alternative (efficient) matricial version of RBF Kernel as implemented by Taylor & Christiannini
        n = size(x,1);
        K = gammac*(x*x.');
        d = diag(K);
        K = K - ones(n,1)*d.'/2 - d*ones(1,n)/2; % Ensures zeros in main diag keeping symmetry
        K = exp(K);     

        alpha = (K + lambda*eye(size(x,1)))\y;
        yhat = alpha.'*K;

    else 
        disp(['Input: ' option ' is not implemented yet']);
    end
    
    yhat = yhat.';
    rmse = sqrt(1/length(x)*norm(y-yhat));
return