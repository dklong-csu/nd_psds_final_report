%   Probability density function for normal
function pdf = fcn_MVN2D(D,L,mu,sigma)
    x = D(:);
    y = L(:);
    R = [x,y];
    psd = mvnpdf(R,mu,sigma);
    [r,c] = size(D);
    pdf = reshape(psd,[r,c]);
end