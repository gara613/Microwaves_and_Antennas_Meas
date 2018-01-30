% function Scorr=normalCorr(Smeas,Snorm,Sref)
% performs the reflection coefficient correction according to the relation:
% Scorr = (Smeas-Sref)./(Snorm-Sref);
% Inputs:
%   - Smeas: S parameters of the DUT
%   - Snorm: S parameters of the Normalization case
%   - Sref: Sparameters of the reference case (Free of interference)
% Outputs:
%   - Scorr: Corrected S parameters
function Scorr=normalCorr(Smeas,Sref,Snorm)
    if size(Smeas)~=size(Snorm)||size(Smeas)~=size(Sref)
        error('Incongruent data sizes');
    end 
	Scorr = (Smeas-Sref)./(Snorm-Sref);
end