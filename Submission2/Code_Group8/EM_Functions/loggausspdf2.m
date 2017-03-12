function y=loggausspdf2(X,mu,cov) % log guassian, mu(mean) cov(covariance)
d = size(X, 1);
X = bsxfun(@minus, X, mu);
[U, p]= chol(cov);
if p ~= 0
    error('ERROR: cov is not SPD.');
end
Q = U'\X;
q = dot(Q, Q, 1);
c = d * log(2 * pi) + 2 * sum(log(diag(U)));
y = -(c + q) / 2;
