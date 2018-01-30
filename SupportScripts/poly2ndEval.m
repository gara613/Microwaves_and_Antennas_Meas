% out = poly2ndEval(coefs,[x1,...,xm])
%
% Evaluate an m-variable second order polynomial at the given abscissa values.
%
% Inputs: 
%     - coefs: m-array polynomial specification as output by interp2Order.
%     - X (x1,..., xm): conformable arrays of abscissa points for polynomial evaluation.
%
% Germán Augusto Ramírez
% Universidad Nacional de Colombia, 2016

function out = poly2ndEval(coefs,X)
	[l,m]=size(X);
    r=m*(m+3)/2+1;
    % validate input data is consistent
	if r~=size(coefs,1)
        error('Input grid data must be conformable to coefs');
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
    %% evaluate hypothesis
    out.eval=A*coefs;
return