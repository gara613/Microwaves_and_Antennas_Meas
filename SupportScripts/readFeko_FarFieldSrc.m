%% [ret_field, thetaGrid, phiGrid, sim_freq, type] = readFeko_FarFieldSrc(varargin)
%
% Routine to read the field pattern in a Feko FarField Source (*.ffe) file
% Format: 8
%   Theta, Phi, Re(E_Theta), Im(E_Theta), Re(E_Phi), Im(E_Phi), Gain(Theta), Gain(Phi), Gain(Total)
%
% Inputs: 
%   - path to file
%   - (optional) scalar: string to request only the scalar field magnitude 
% Outputs: 
%	- ret_field: Returned field, structure with fields:
%       * E_theta
%       * E_phi
%       * E_mag
%       * abs_radPat (dB): 10log10( (|E_th|^2 + |E_ph|^2)/(60*Prad) ),
%       + NOTE: Theta and Phi components of abs_radPat (in the file) are not returned
%   - thetaGrid: theta values arranged as a proper matlab grid 
%   - phiGrid: phi values arranged as a proper matlab grid 
%   - sim_freq: Frequency specified in the file
%   - type: of the returned magnitude (Directivity, Gain, Realized Gain)
%       + NOTE: ret_field is the electric field 
% 
% EPFL, MAG, Feb 2023
% Germán Augusto Ramírez Arroyave

function [ret_field, thetaGrid, phiGrid, sim_freq, type] = readFeko_FarFieldSrc(varargin)
    narginchk(1,2);
    filename = char(varargin{1});
    if nargin == 2	% kind is currently unimplemented, use to limit the size of the returned structure
        scalar = varargin{2};
    else
        scalar = false;
    end
    % Initialize default return values
    if scalar
        ret_field = struct('abs_radPat',[]);
    else 
        ret_field = struct('E_theta',[],'E_phi',[],'E_mag',[],'abs_radPat',[]);
    end
       
    %% Read datafile  
    fileId = fopen(filename, 'rt');
	[sim_freq, Npts, type] = readHdrFFE(fileId);        % reads first rows of file
    radPat = fscanf(fileId,'%f', [9, inf]).';           % reads data columnwise into a (9,:) matrix, transposed to restore file appearance
	fclose(fileId);

    %% Return field and angular sweep grid
    thetaGrid = reshape(radPat(:,1)*pi/180, Npts.theta, Npts.phi);
    phiGrid = reshape(radPat(:,2)*pi/180, Npts.theta, Npts.phi);

    if scalar
        ret_field.abs_radPat = reshape(radPat(:,9), Npts.theta, Npts.phi);
    else
        e_theta = radPat(:,3) + 1i*radPat(:,4);
        e_phi = radPat(:,5) + 1i*radPat(:,6);
        ret_field.E_theta = reshape(e_theta, Npts.theta, Npts.phi);
        ret_field.E_phi = reshape(e_phi, Npts.theta, Npts.phi);
        ret_field.E_mag = sqrt(abs(ret_field.E_theta).^2 + abs(ret_field.E_phi).^2);
        ret_field.abs_radPat = reshape(radPat(:,9), Npts.theta, Npts.phi); 
    end
end

% File header reading
function [sim_freq, Npts, type] = readHdrFFE(fileId)
	while feof(fileId) == 0
        linea_actual = fgets(fileId);
        if ~isempty(regexp(linea_actual,'\w*Etheta\w*', 'match', 'ignorecase'))% ~isempty(regexp(linea_actual,{'\w*Etheta\w*','\w*Ephi\w*'}, 'match', 'ignorecase'))
            break
        end    
        % Manually read the values...
        linea_actual = regexp(linea_actual,':','split');
        if strcmp(linea_actual{1},'#Frequency')
            sim_freq = str2double(linea_actual{2});
        end
        if strcmp(linea_actual{1},'#No. of Theta Samples')
            Npts.theta = str2double(linea_actual{2});
        end
        if strcmp(linea_actual{1},'#No. of Phi Samples')
            Npts.phi = str2double(linea_actual{2});
        end
        if strcmp(linea_actual{1},'#Result Type')
            type = linea_actual{2};
        end
    end       
end