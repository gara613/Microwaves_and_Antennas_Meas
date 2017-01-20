% [ABCD] = function myS2ABCD(S, Z0) 
% Rutina para convertir parámetros S a ABCD
% Recibe la matriz de parámetros S (Spars entregada por la rutina readSpars) para cada una de las frecuencias 
% Se asume que la impedancia característica es independiente de la frecuencia
%
% Germán Augusto Ramírez Arroyave
% CMUN - 2014
function [ABCD] = myS2ABCD(S, Z0)

    narginchk(1,2)
    if nargin < 2
        Z0 = 50;
    end
    if any(Z0 == 0)
        error('The characteristic impedance cannot be zero');
    end

    % Implementación poco eficiente en memoria en favor de un código más claro, 
    % en caso de que las matrices sean realmente muy grandes habrá que ajustarlo
    S11 = deal(S(1,1,:));	S12 = deal(S(1,2,:));
    S21 = deal(S(2,1,:));   S22 = deal(S(2,2,:));

    den = 2*S21;

    ABCD(1,1,:) = ((1 + S11).*(1 - S22) + S12.*S21) ./den;
    ABCD(1,2,:) = Z0 .* ((1 + S11).*(1 + S22) - S12.*S21) ./den;
    ABCD(2,1,:) = 1./Z0 * ((1 - S11).*(1 - S22) - S12.*S21) ./den;
    ABCD(2,2,:) = ((1 - S11).*(1 + S22) + S12.*S21) ./den;
return 