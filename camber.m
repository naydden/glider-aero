% Function camber
% It calculates the z-coordinates of the mean camber line

function z = camber(x, m, p)
if(x<p)
    z = m*(2*p*x-x^2)/p^2;
else
    z = m*(1-2*p+2*p*x-x^2)/(1-p)^2;
end
end