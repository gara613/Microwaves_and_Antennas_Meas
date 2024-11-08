%% [ret_field, thetaGrid, phiGrid, P_ant, sim_freq] = readCST_FarFieldSrc(varargin)
%
% Routine to read the field pattern in a CST FarField Source (*.ffs) file
% Format: V3.0 
%   Phi, Theta, Re(E_Theta), Im(E_Theta), Re(E_Phi), Im(E_Phi)
%
% Inputs: 
%   - path to file
%   - (optional) kind: string to request raw E field components or just the Power pattern
% Outputs: 
%	- ret_field: Returned field, structure with fields:
%       * E_Field.E_theta
%       * E_Field.E_phi
%       * abs_radPat (dB): 10log10( (|E_th|^2 + |E_ph|^2)/(60*Prad) )
%   - thetaGrid: theta values arranged as a proper matlab grid 
%   - phiGrid: phi values arranged as a proper matlab grid 
%   - P_ant: Radiated/Accepted/Stimulated Power/radiation Efficiency/total Efficiency/max realized Gain
%   - sim_freq: Frequency specified in the file
% 
% Uiversidad Nacional de Colombia
% Germán Augusto Ramírez Arroyave
% CMUN, 2018

function [ret_field, thetaGrid, phiGrid, P_ant, sim_freq] = readCST_FarFieldSrc(varargin)
    narginchk(1,2);
    filename = char(varargin{1});
    if nargin == 2	% kind is currently unimplemented, use to limit the size of the returned structure
        kind = varargin{2};
    else
        kind = false;
    end
    % Initialize default return values
    if kind
        ret_field = struct('E_Field',struct(kind,[]));
    else 
        ret_field = struct('E_Field',struct('E_theta',[],'E_phi',[]),'abs_radPat',[]);
    end
    thetaGrid = []; 
    phiGrid = [];
    P_ant = [];
    sim_freq = [];
    
    %% Read datafile  
    fileId = fopen(filename, 'rt');
	[P_ant, sim_freq, Npts] = readHdrFFS(fileId);	% reads first rows of file
    radPat = fscanf(fileId,'%f', [6, inf]).';       % reads data columnwise into a (6,:) matrix, transposed to restore file appearance
	fclose(fileId);

    %% Return field and angular sweep grid
    phiGrid = reshape(radPat(:,1)*pi/180, Npts.theta, Npts.phi);
    thetaGrid = reshape(radPat(:,2)*pi/180, Npts.theta, Npts.phi);

	multTerm = (1/60)*(1/P_ant.Psource);    %  = 0.5*(sphFac/eta_0), (eta_0 = 120*pi, sphFac = 4*pi)

    e_theta = radPat(:,3) + 1i*radPat(:,4);
    e_phi = radPat(:,5) + 1i*radPat(:,6);
%     ret_field.E_Field.E_theta = reshape(e_theta, Npts.theta, Npts.phi);
%     ret_field.E_Field.E_phi = reshape(e_phi, Npts.theta, Npts.phi);
    ret_field.E_Field.E_theta = sqrt(multTerm)*reshape(e_theta, Npts.theta, Npts.phi);
    ret_field.E_Field.E_phi = sqrt(multTerm)*reshape(e_phi, Npts.theta, Npts.phi);

    abs_radPat = multTerm*(e_theta.*conj(e_theta) + e_phi.*conj(e_phi));
    abs_radPat = 10*log10(abs_radPat);
    ret_field.abs_radPat = reshape(abs_radPat, Npts.theta, Npts.phi);

    P_ant.maxRdGain = max(max(abs_radPat));
end

% File header reading
function [P_ant, sim_freq, Npts] = readHdrFFS(fileId)
	while feof(fileId) == 0
        linea_actual = fgets(fileId);
        if strcmp(linea_actual(1:end-1),'// >> Phi, Theta, Re(E_Theta), Im(E_Theta), Re(E_Phi), Im(E_Phi): ')
            break
        end    
        % Manually read the values...
        linea_actual = regexprep(linea_actual, ' +\s', ' ');

        if strcmp(linea_actual,'// Radiated/Accepted/Stimulated Power , Frequency ')
            P_ant.Prad = str2double(fgets(fileId));
            P_ant.Pin = str2double(fgets(fileId));
            P_ant.Psource = str2double(fgets(fileId));
            sim_freq = str2double(fgets(fileId));

            P_ant.radEff = 10*log10(P_ant.Prad/P_ant.Pin);
            P_ant.totEff = 10*log10(P_ant.Prad/P_ant.Psource);
        end

        if strcmp(linea_actual(1:end-1),'// >> Total #phi samples, total #theta samples')
%             linea_actual = fgets(fileId);
%             nnpts = str2num(linea_actual); 
            nnpts = str2num(fgets(fileId));            
            Npts.phi = nnpts(1);
            Npts.theta = nnpts(2);
        end
	end
end