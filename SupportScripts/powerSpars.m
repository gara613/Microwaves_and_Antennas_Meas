function P = powerSpars(Smat)
% Returns the power of the signal with the minimum power in the S parameters matrix Smat

P = zeros(size(Smat,1)^2,1);
cont=1;
for i = 1:size(Smat,1)
    for j = 1:size(Smat,2)
        % power in the spectrum of a signal a: Pa = a'*a (operator ' = *T = Hermitian or conjugate transpose)
        P(cont) = norm(squeeze(Smat(i,j,:)))^2;
        cont=cont+1;
    end
end
P = min(P); % power could be defined for each parameter and hence returned as a vector

% or even given as a mean value, but this masks the power of the smaller signal
% P_s = 0;
% for i = 1:size(Smat,1)
%     for j = 1:size(Smat,2)
%         P_s = P_s + norm(squeeze(Smat(i,j,:)))^2;
%     end
% end
% N = size(Smat,1)*size(Smat,2)*size(Smat,3);
% P_s = P_s/N;