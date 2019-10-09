% Design considerations of a microstrip transmission line for non-linear
% diodes characterization

c = 3e8;
epsr = [2.2; 3.0; 3.0; 3.0; 3.55; 3.55; 3.55; 3.66; 3.66; 11.2; 11.2; 11.2];
tand = 0.0009;%0.0027;
h = 1e-3*[0.125; 0.25; 0.75; 1.52; 0.508; 0.813; 1.524; 0.254; 0.762; 0.25; 0.64; 1.28];
f_max = 50e9;
Zair = 120*pi;
Z0 = 50; % Intended characteristic impedance

A = Z0/60 * sqrt((epsr+1)/2) + (epsr-1)./(epsr+1).*(0.23+0.11./epsr);
B = Zair*pi./(2*Z0*sqrt(epsr));

w = zeros(size(epsr));
w = h.*8.*exp(A)./(exp(2*A)-2) .* (8*exp(A)./(exp(2*A)-2) < 2);   % w/h < 2
w = w + h.*2./pi.* ( B-1-log(2*B-1) +(epsr-1)./(2*epsr).*(log(B-1) +0.39 -0.61./epsr) ) .* (8*exp(A)./(exp(2*A)-2) > 2);    % w/h > 2
u = w./h;

epsr_eff = 0.5*(epsr+1) + 0.5*(epsr-1)./sqrt(1+12./u); % Pozar's book, use HJ expressions if more accuracy needed 
lambda_eff = c./(sqrt(epsr_eff)*f_max);

fu = 6 + (2*pi-6)*exp(-(30.666./u).^0.7528);
Z_L = 60./sqrt(epsr_eff).* log(fu./u + sqrt(1+(2./u).^2)); % HJ accurate expresssion   

h_rec = 0.1*lambda_eff;  % Max recommended substrate height

f_s = c*atan(epsr)./(sqrt(2)*pi*h.*sqrt(epsr-1));  % Cut frequency for surface waves occurrence
f_c = c./(sqrt(epsr).*(2*w+0.8*h));  % Cut frequency for higher order modes emergence