% function [yhat, ord, rmse] = rationalSpar(x, y)
% Function to fit a data-set using rational interpolation
%
% outputs: 
%   - yhat   = estimated vector, 
%   - ord    = order of numerator (order of denominator=m-1), 
%   - rmse
% inputs: 
%   - x      = values of the independent variable
%   - y      = values of the dependent variable (dataset)
%   - option = {'barycentric', 'canonical'} representation of the rational function
% 
% Germán Augusto Ramírez Arroyave
% Universidad Nacional de Colombia
% CMUN 2014

function [yhat, orNum, rmse] = rationalSpar(x, y, varargin)
    %% Default values and initialization of parameters
    option = 'canonical';
    maxOrd = 1;
    % normalize sampling points 
	% xnorm = x./max(x);                % Normalizes to [0 1]
    magOrderX = floor(log10(max(x)));   % Another kind of normalization  
    xnorm = x./10^(magOrderX);          % Removes the exponent
    yhat = zeros(size(y));

	tol = 1e-2;             % Maximum approximation error (high to avoid overfitting)
    xtam = length(x);
    delta = 10*eps; 
    num = zeros(xtam,1);
    den = zeros(xtam,1);
    orNum = 1; orDen = 0;   % Assume orNum = orDen+1, increase sequentially in steps of one
    if ~isempty(varargin)
        option = varargin{1}; 
        maxOrd = varargin{2};
    end
    %% Berrut & Mittelmann. "Matrices for the direct determination of the barycentric weights of rational interpolation" 
    % Journal of Computational and Applied Mathematics, 78 355-370, 1997.
    if strcmp(option, 'barycentric')
        while (norm(y-yhat)^2/norm(y)^2 > tol) && (orDen < maxOrd)          
            % Sample the input vector with Nsamp samples
            Nsamp = orNum + orDen + 1; % Limit number of samples (condition for existence of nullspace of A)        
            % sampInd = randperm(xtam,Nsamp); % random sampling doesn't seem to work very well
            % newIdx = sort(sampInd);         
            newIdx = floor(linspace(1,xtam,Nsamp)); % uniform sampling behaves better than random
            xsamp = xnorm(newIdx);              
            ysamp = y(newIdx);              % this downsampled representation is enough for reconstruction

            % Create matrix 
            A = ones(size(xsamp)).';        
            for cont = 1:orNum-1
                A = [A; (xsamp.^cont).'];
            end
            for cont = 0:orDen-1
                A = [A; (ysamp.*xsamp.^cont).'];
            end
            
            u = null(A); % Values of parametrization in barycentric representation
            % Reconstruction of data with barycentric representation       
            for cont = 1:xtam
                num(cont) = sum(ysamp.*u(:,1)./(xnorm(cont)+delta - xsamp));
                den(cont) = sum(u(:,1)./(xnorm(cont)+delta - xsamp));
            end
            yhat = num./den;
            orNum = orNum+1;
            orDen = orDen+1;
        end
    %% Canonical representation of rational functions
   	elseif strcmp(option, 'canonical')
        % Similar idea can be found on: "Automated Fitting and Rational Modeling Algorithm for EM-Based S-Parameter Data" by Dhaene
        while (norm(y-yhat)^2/norm(y)^2 > tol) && (orDen < maxOrd) % Order must be limited to avoid overfitting and poles
            xsamp = xnorm;  % no down-sampling is performed
            ysamp = y;
            
            A = ones(size(xsamp));
            
            for cont = 1:orNum
                A = [A xsamp.^cont];
            end
        	for cont = 1:orDen
                A = [A -ysamp.*xsamp.^cont];
        	end
            
            coefs = A\ysamp; % Values of parametrization in canonical representation fixing b0 = 1
            % Reconstruction of data with canonical representation
            num = polyval(coefs(orNum+1:-1:1),xnorm);
        	den = polyval([coefs(end:-1:orNum+2); 1],xnorm);
        	yhat = num./den;
            orNum = orNum+1;
            orDen = orDen+1;
        end
    end
    rmse = sqrt(1/xtam*norm(y-yhat));
end