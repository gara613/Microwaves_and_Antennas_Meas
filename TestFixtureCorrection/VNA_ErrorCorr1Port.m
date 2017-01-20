% Rutina para corregir errores de medición de un puerto a partir del modelo de error de tres términos 
% descrito en la nota de aplicación de Agilent: "Network Analyzer Error Models and Calibration Methods".

%% Mediciones de referencia: corto, abierto y carga de 50 Ohm.
[GammaM(1,:), freq_Open, R_Open] = readSpars('Open_Cal.s1p');
[GammaM(2,:), freq_Short, R_Short] = readSpars('Short_Cal.s1p');
[GammaM(3,:), freqs, Z0] = readSpars('Load_Cal.s1p');

if freq_Open ~= freq_Short | freq_Short ~= freqs
    error('las frecuencias de las mediciones no coinciden');
end
if R_Open ~= R_Short || R_Short ~= Z0
    warning('las mediciones no se tomaron con la misma impedancia de referencia');
end

%% Corrección de mediciones
% Mediciones que serán corregidas
listaind = {'Induc6p8nH.s1p', 'Induc10nH.s1p', 'Induc68nH.s1p'};
listacap = {'Cap330Taiyo.s1p', 'CapJD220pF.s1p', 'CapJT10pF.s1p', 'CapJT22pF.s1p', 'CapPanasonic560pF.s1p', 'CapTaiyo10pF.s1p', 'CapTDK220pF.s1p'};
% Llamar a la rutina de corrección
Z = correction1Port(GammaM,freqs,Z0,listaind,1,'L');
Z = correction1Port(GammaM,freqs,Z0,listacap,0,'C');