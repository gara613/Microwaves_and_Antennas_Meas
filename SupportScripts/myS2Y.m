% [Z] = function myS2Z(S, Z0) 
% Rutina para convertir parámetros S a Z
% Recibe la matriz de parámetros S (Spars entregada por la rutina readSpars) para cada una de las frecuencias 
% Se asume que la impedancia característica es independiente de la frecuencia
%
% Germán Augusto Ramírez Arroyave
% CMUN - 2018
function [Y] = myS2Y(S, Z0)

    narginchk(1,2)
    if nargin < 2
        Z0 = 50;
    end
    if any(Z0 == 0)
        error('The characteristic impedance cannot be zero');
    end
    Y0 = inv(Z0);
	nports = size(S,1);
    Y = zeros(size(S));
%    Z0 = Z_ref*eye(totN_ports);
%    G0 = 1/Z_ref*eye(totN_ports);

    if nports == 2
        S11 = deal(S(1,1,:));	S12 = deal(S(1,2,:));
        S21 = deal(S(2,1,:));   S22 = deal(S(2,2,:));

        den = (1+S11).*(1+S22)-S12.*S21;              
        
        Y(1,1,:) = Y0*( (1-S11).*(1+S22) + S12.*S21 ) ./den;
        Y(1,2,:) = -Y0*2*S12 ./den;
        Y(2,1,:) = -Y0*2*S21 ./den;
        Y(2,2,:) = Y0*( (1+S11).*(1-S22) + S12.*S21 ) ./den;
    else
        if isscalar(Z0)
            Y0 = 1/Z0*eye(nports);
        end
        for cont = 1:size(S,3)
            %Y(:,:,cont) = G0\ (S(:,:,cont)*Z0 + conj(Z0))\ (eye(nports) - S(:,:,cont)) *G0;    
            Y(:,:,cont) = Y0* inv(S(:,:,cont) + eye(nports))* (eye(nports) - S(:,:,cont));
        end
    end
return 