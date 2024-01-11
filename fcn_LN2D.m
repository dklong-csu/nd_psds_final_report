%   Probability density function for log-normal
function pdf = fcn_LN2D(D,L,mu,sigma)
    X = log(D);
    Y = log(L);
    x = X(:);
    y = Y(:);
    R = [x,y];
    psd = mvnpdf(R,mu,sigma)./exp(x)./exp(y);
    [r,c] = size(X);
    pdf = reshape(psd,[r,c]);
end