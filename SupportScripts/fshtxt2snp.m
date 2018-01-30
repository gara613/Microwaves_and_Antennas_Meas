% function fshtxt2snp(route,option)
% function to convert a *.set (.txt) FSH data file to a *.snp file
% inputs:
%   - route: path to the *.txt file (for two ports the following parameters
%            order is assumed: S11, S22, S21, S12, please ensure this in the VNA)
%   - option: text string to indicate if one ('one') or two ('two') ports measurements are reported
% outputs:
%   - A *.snp is created in the same directory where the txt file is located
%   
% Germán Augusto Ramírez Arroyave
% CMUN - Universidad Nacional de Colombia 2016
function fshtxt2snp(route,option)

%     %% second option seems a better alternative-...
%     comma2point_overwrite( filespec )
%     file    = memmapfile( filespec, 'writable', true );
%     comma   = uint8(',');
%     point   = uint8('.');
%     file.Data( transpose( file.Data==comma) ) = point;

end