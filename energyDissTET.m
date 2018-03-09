function [E,A] = energyDissTET(xi_na,Z00,kappa,p)

% energyDissTET returns the reduction of the main system's dimenionless
%   energy  and amplitude as defined in the research paper
% "Performance Measures for Targeted Energy Transfer and Resonance
% Capture Cascading in Nonlinear Energy Sinks" by Kevin Dekemele et al,
% equation (28) and (29). DOI: 10.1007/s11071-018-4190-5
%
%
%       [E,A] = energyDissTET(xi_na,Z00,kappa)
%
%               E      Reduction of main system dimensionless energy [fraction]
%               A      Reduction of main system amplitude  [fraction]
%
%               xi_na   Dimensionless damping of absorber
%               Z00     Initial dimensionless 'energy'
%               kappa   Dimensionless linear stiffness of absorber
%               p       Power of absorber
%
% Copyright (c) 2018, Kevin Dekemele (kevindekemele@gmail.com)
% Source code available at 
% Licensed under GNU GPLv3 license.
%

b = nchoosek(p,(p-1)/2);
%% maximum and minumum SIM

ZnaM = nthroot(2^(p-2)*((1-kappa)*(p+1)-sqrt((p-1)^2*(1-kappa)^2-4*p*xi_na^2))/(b*p),(p-1)/2);
ZnaP = nthroot(2^(p-2)*((1-kappa)*(p+1)+sqrt((p-1)^2*(1-kappa)^2-4*p*xi_na^2))/(b*p),(p-1)/2);

Z0P = (xi_na^2+(1-kappa-b/(2^(p-1))*ZnaM^((p-1)/2))^2).*ZnaM;
Z0M = (xi_na^2+(1-kappa-b/(2^(p-1))*ZnaP^((p-1)/2))^2).*ZnaP;
%% Calculate E_diss and A_Reduction
E = 1-Z0M/Z00;

A = 1-sqrt(Z0M/Z00);

end

