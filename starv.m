function [ S ] = starv( x, dmin, dmax, alpha2 )


S = dmax - (( dmax - dmin)*exp(-(alpha2*x)));
end

