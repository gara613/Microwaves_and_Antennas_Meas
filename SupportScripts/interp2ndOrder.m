% function out = interp2ndOrder(X,Y)
% 
% Provide a multivariate second degree polynomial fit to the data Y
% given the m features X specified. 
% Y = a_{0} + a_{1}*x_{1} +...+ a_{m}*x_{m} + a_{12}*x_{1}x_{2} +...+ a_{1m}*x_{1}x_{m} +...+ a_{m1}*x_{m-1}x_{m} + a_{m+1}*x_{1}^2 +...+ a_{r}*x_{m}^2
% where r = 1+m*(m+3)/2 is the resulting number of coefficients
%
% Inputs:
%   x1,...,xm: m input vectors specifying the l evaluation points (X \in \mathbb{R}^{l\times m}) (m: number of features, l: number of samples)
%   y1,...,yn: array of function samples evaluated at the input locations
%              specified by y1,...,yn (Y \in \mathbb{R}^{l\times n}) (n: number of outputs, l: number of samples)
%
% Output:
%   out: data structure with fields:
%       - coefs: r x n array specifying the weight of terms 
%       - relerror: relative fitting error
%       - expr: string expression for the fitting polynomial using the variables 'x1',...,'xm'
%
% Germán Augusto Ramírez
% Universidad Nacional de Colombia, 2016

function out = interp2ndOrder(X,Y)
	[l,m]=size(X);
	[ll,n]=size(Y);
    r=m*(m+3)/2+1;
    % validate input data is consistent
	if l~=ll
        error('Input grid data must be conformable to output or have length(xn) = size(fval,n)');
	end

    %% create features powers matrix
	A=ones(l,r);
	cont=1+m;
	for cont1=1:m 
		A(:,cont1+1)=X(:,cont1);
	    A(:,end-m+cont1)=X(:,cont1).^2;
        for cont2=cont1+1:m
            cont=cont+1;
            A(:,cont)=X(:,cont1).*X(:,cont2);
        end
	end
    
    %% Calculate coefficients
%     alpha=1e-3;
%     G=alpha*eye(r); % regularization matrix (Nfeatures x Nfeatures)
%     out.coefsN=(A'*A+G'*G)\A'*Y; %regularized normal equation solution (Ordinary Least Squares) this may have problems given rank(X)<r
%     [Q,R]=qr(A,0);
%     out.coefsQ=R\Q'*Y; % normal equation solution (QR factorization) may have problems given rank(X)<r
    out.coefs=A\Y; %Straight approach (Matlab solves the linear problem by the most suitable way)

    out.relerror=norm(A*out.coefs-Y)/norm(Y);
    
    %% build literal expression for easy interpretation
    [cont1,cont2]=find(tril(ones(m,m),-1));
	out.expr = cell(n,1);
	for countr=1:n
        for cont=1:r
            if out.coefs(cont,countr) ~= 0 % put sign and coefs into expression
                if out.coefs(cont,countr) < 0
                    pref=sprintf(' - %g',abs(out.coefs(cont,countr)));
                elseif isempty(out.expr) %positive, first term, insert without sign
                    pref=sprintf('%g',abs(out.coefs(cont,countr)));
                else
                    pref=sprintf(' + %g',abs(out.coefs(cont,countr)));
                end
            else
                pref = '';
            end

            if ~isempty(pref) % put variables into expression
                if cont>1&&cont<m+2
                    pref = [pref,'*X',num2str(cont-1)];
                elseif cont>m+1&&cont<=r-m
                    pref = [pref,'*X',num2str(cont2(cont-m-1)),'*X',num2str(cont1(cont-m-1))];
                elseif cont>r-m
                    pref = [pref,'*X',num2str(cont-r+m),'^2'];
                end

                out.expr{countr} = [out.expr{countr},pref];
            end
        end   
	end
return 