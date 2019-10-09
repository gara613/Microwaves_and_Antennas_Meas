% Homework: Antenna Phased Arrays
% 
% Germán Augusto Ramírez Arroyave
% UPC - 2019

function D = arrayDir(In,kd,alpha,varargin)
    N = length(In);
    dubleterms = zeros(size(alpha));
    for n = 1:N-1 
        for q = n+1:N
            dubleterms = dubleterms + In(n)*In(q) * bsxfun(@times, sin(kd*(n-q))./(kd*(n-q)), cos(alpha*(n-q)));
        end
    end
    den = sum(In.^2) + 2*dubleterms;
    D = sum(In)^2./den;
end