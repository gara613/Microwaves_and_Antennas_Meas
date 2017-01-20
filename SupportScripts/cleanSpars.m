% function Spars = cleanSpars(Smatrix, freqs, method, parameters)
% 
% Inputs: 
%   - Smatrix: S matrix of dimensions m,n,l, where l is the number of frequency samples
%   - freqs: The frequencies sampling vector
%   - method: 
%       - 'Picewise'
%       - 'RationalS'
%       - 'KernelS'
%   - parameters: configuration vector for the chosen method
%       - 'Picewise'
%           - Weigth: 'lineal' or 'triangle'
%           - Nparts: number of segments in the freqs sweep
%       - 'RationalS'
%           - Representation: ('canonical' or 'barycentric')
%           - Maxorder: (to avoid overfitting)
%       - 'KernelS'
%           - Representation: ('polynom' or 'gaussi')
%           - lambda: regularization to avoid overfitting
%           - Maxorder: (for polynomial), for Gaussian kernel parameter gamma 
%                       is calculated as gammac = 1/Maxorder
%
% Germán Augusto Ramírez Arroyave 
% Universidad Nacional de Colombia
% CMUN 2014

function Spars = cleanSpars(Smatrix, freqs, method, varargin)
    % Use varargin to determine aditional parameters for each method
    [m,n,r] = size(Smatrix);
    Spars = zeros(m,n,r);
    if r ~= length(freqs)
        disp('error: incoherent frequency sampling'); 
        return
    end
    % This option calls to the Nparts picewise interpolation script mySimpleFit
    if strcmp(method,'Picewise')
        if ~isempty(varargin)
            weight = varargin{1};
            Nparts = varargin{2};
        else
            weight = 'triangle';
            Nparts = 5;
        end
        for cont1 = 1:m % Every S_{ij} is passed as an independet vector to the fitting function
            for cont2 = 1:n
                %Sij = squeeze(Smatrix(cont1,cont2,:)); 
                Spars(cont1,cont2,:) = mySimpleFit(freqs, squeeze(Smatrix(cont1,cont2,:)), weight, Nparts);
            end
        end
    % This option calls to the rational interpolation function rationalSpar
    elseif strcmp(method,'RationalS')
        if ~isempty(varargin)
            Representation = varargin{1};
            Maxorder = varargin{2};
        else
            Representation = 'canonical';
            Maxorder = 5;
        end
        for cont1 = 1:m 
            for cont2 = 1:n
                Spars(cont1,cont2,:) = rationalSpar(freqs, squeeze(Smatrix(cont1,cont2,:)), Representation, Maxorder);
            end
        end
	% This option calls to the Kernel interpolation function kernelSpar
    elseif strcmp(method,'KernelS')
        if ~isempty(varargin)
            Representation = varargin{1};
            lambda = varargin{2};
            if strcmp(Representation,'polynom')
                Maxorder = varargin{3};
            else%if strcmp(Representation,'gaussi')
                Maxorder = varargin{3}; % in case there is a different way to establish gamma
            end
        else
            Representation = 'polynom';
            lambda = 0.5;
            Maxorder = 7;
        end
        for cont1 = 1:m 
            for cont2 = 1:n
                Spars(cont1,cont2,:) = kernelSpar(freqs, squeeze(Smatrix(cont1,cont2,:)), Representation, lambda, Maxorder);
            end
        end
    else
        disp('error: Unimplemented interpolation method');
    end
end