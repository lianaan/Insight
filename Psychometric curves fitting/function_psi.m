function f = function_psi(s,mu,sigma,lambda)
%PSYCHOFUN Psychometric function cumulative density of Gaussian

f = bsxfun(@plus, lambda/2, ...
    bsxfun(@times,1-lambda,0.5*(1+erf(bsxfun(@rdivide,bsxfun(@minus,s,mu),sqrt(2)*sigma)))));

end
