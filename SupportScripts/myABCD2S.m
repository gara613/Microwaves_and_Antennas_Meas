% [S] = function myABCD2S(ABCD, Z0) 
% Rutina para convertir parámetros ABCD a S
% Debe recibir la matriz de parámetros ABCD para cada una de las frecuencias 
%
% Germán Augusto Ramírez Arroyave
% CMUN - 2014
function [S] = myABCD2S(ABCD, Z0)

    narginchk(1,2)
    if nargin < 2
        Z0 = 50;
    end
    if any(Z0 == 0)
        disp('Error: the characteristic impedance can not be zero');
    end

    % Implementación poco eficiente en memoria en favor de un código más
    % claro, ajustar en caso de que las matrices sean realmente muy grandes
    A = deal(ABCD(1,1,:));  B = deal(ABCD(1,2,:));
    C = deal(ABCD(2,1,:));  D = deal(ABCD(2,2,:));

    den = (A + B./Z0 + C.*Z0 + D);

    S(1,1,:) = (A + B./Z0 - C.*Z0 - D) ./ den;
    S(1,2,:) = 2*(A.*D - B.*C) ./ den;
    S(2,1,:) = 2./ den;
    S(2,2,:) = (-A + B./Z0 - C.*Z0 + D) ./ den;
return 