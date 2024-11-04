% function val = writeSnp(filePath,fileName,rawData,format,varargin)
%
% Trivial function to write either s1p or s2p data to a file, pending
% support for higher order snp, and for noise data
%
% Inputs: 
%   * filePath: a valid path for storing the file, if folder doe snot exist, it is created
%   * fileName: a valid name for the file
%   * rawData: Can be either
%       - An array with s(1/2)p format (i.e., <freq> <S11_re> <S11_im> <S21_re> <S21_im> <S12_re> <S12_im> <S22_re> <S22_im>) or
%       - A cell with two arrays, the first one is the frequency, while the second is the conventional 3D S-parameters matrix (i.e., N_pts x N_pts x N_freqs)
%   * format: RI/MA, it is assumed that the rawData provided is a complex
%   number, real/imaginary in the columns case, and complex in the 3D array
%   case
%   * R_car (Optional): characteristic impedance of the measurement system 
%
% Outputs:
%   * val: Flag indicating if everything was fine
% 
% GermÃ¡n Augusto RamÃ­rez 
% Universitat PolitÃ¨cnica de Catalunya, 2020

function val = writeSnp(filePath,fileName,rawData,format,varargin)
    if ~isempty(varargin)
        if length(varargin)>2
            R_car = varargin{1};
            comments = varargin{2};       
            silentWrite = varargin{3};
        elseif length(varargin)>1
            R_car = varargin{1};
            comments = varargin{2};       
        else
            R_car = varargin{1};
        end
    else
        R_car = 50; 
        comments = ''; 
        silentWrite = false;
    end
    
%% Arrange data into a matrix if not already provided
	if isa(rawData,'cell')
        freqs = rawData{1};
        nFreqs = size(freqs,1);
        sPars = rawData{2};
        nPorts = size(sPars,1);        
        rawData = zeros(nFreqs,2*nPorts^2+1);
        rawData(:,1) = freqs;
        rawData(:,2:2:end) = reshape(real(sPars),nPorts^2,nFreqs).';
        rawData(:,3:2:end) = reshape(imag(sPars),nPorts^2,nFreqs).';     
        
	else
        nPorts = sqrt((size(rawData,2)-1)/2);
	end
    
    %% Convert data to specified format
    if strcmpi(format,'ma')	% Convert RI columns to MA
        mag_lin = abs(rawData(:,2:2:end).^2+rawData(:,3:2:end).^2);
        pha_rad = atan2(rawData(:,3:2:end),rawData(:,2:2:end));
        rawData(:,2:2:end) = mag_lin;
        rawData(:,3:2:end) = pha_rad*180/pi;
    elseif strcmpi(format,'db')	% Convert RI columns to dB
        mag_dB = 10*log10(abs(rawData(:,2:2:end).^2+rawData(:,3:2:end).^2));
        pha_rad = atan2(rawData(:,3:2:end),rawData(:,2:2:end));
        rawData(:,2:2:end) = mag_dB;
        rawData(:,3:2:end) = pha_rad*180/pi;
    end
    
    %% open file
	if ~exist(filePath, 'dir')
        mkdir(filePath)
    end
    if ~exist([fullfile(filePath,fileName),'.s',num2str(nPorts),'p'], 'file')
        fileId = fopen([fullfile(filePath,fileName),'.s',num2str(nPorts),'p'],'wt');        
    else       
        apnd_str = char(datetime('now','format','yyyyMMdd_hhmm')); % Optional: Recursive call to validate if xx_1/2/3.cff already exists...
        fprintf('%s.cff,  already exists, saving with a different name.\n', fullfile(filePath,fileName))
        fileId = fopen([fullfile(filePath,fileName),apnd_str,'.s',num2str(nPorts),'p'],'wt');
	end

    %% Write file header
    fprintf(fileId,'%s \n','! Microwaves and Antennas Group - MAG');
    fprintf(fileId,'%s \n','! École polytechnique fédérale de Lausanne');
    fprintf(fileId,'%s \n',['! Date: ' char(datetime('now','format','yyyy:MM:dd hh:mm'))]);
    fprintf(fileId,'%s \n',['! S',num2str(nPorts),'p file : Measurements: ' ]);
    fprintf(fileId,'%s \n',['! ',fileName, ' created with writeSnp.m']); 
	if iscell(comments)
        for cont = 1:length(comments)
            fprintf(fileId,'%s \n',['! ',comments{cont}]);
        end
    else
        fprintf(fileId,'%s \n',['! ',comments]);
    end   
    fprintf(fileId,'%s \n',['# Hz  S ', upper(format), ' R ' num2str(R_car)]);
    
    %% Write data
    fomrStr = repmat('%4.6g \t',1,2*nPorts^2+1);
	val = fprintf(fileId ,[fomrStr '\n'],rawData.');
    
    if ~silentWrite
        if val
            disp(['Successfully created: ' fullfile(filePath,fileName),'.s',num2str(nPorts),'p']);
        else
            error(['Something went wrong with file creation, error: ', num2str(val)]);
        end
    end
    fclose(fileId);
end