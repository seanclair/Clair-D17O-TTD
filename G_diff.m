function G = G_diff(gamma,delta,t)
%G_DIFF Calculates the diffusive age spectrum from Waugh and Hall 2002 eqn
%9
%   Inputs: 
%       gamma: mean age
%       delta: spectral width
%       t: elapsed time
%   Output:
%       G: age spectrum with lenght(t) elements

that = t ./ gamma;


G = (1 ./ (2 .* delta .* sqrt(pi .* that.^3))) .* exp((-1 .* gamma.^2 .* (that - 1).^2) ./ (4 .* delta.^2 .* that));



end