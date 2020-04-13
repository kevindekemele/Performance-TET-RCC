function [T_pump] = pumpingTime(xi_na,Z00,kappa,p)
% pumpingTime returns pumping time as defined in the research paper
% "Performance Measures for Targeted Energy Transfer and Resonance
% Capture Cascading in Nonlinear Energy Sinks" by Kevin Dekemele et al,
% equation (33)
%
%       T_pump = pumpingTime(xi_na,Z00,kappa)
%
%               T_pump      Duration of TET [epsilon*s]
%        
%               xi_na   Dimensionless damping of absorber
%               Z00     Initial dimensionless 'energy'
%               kappa   Dimensionless linear stiffness of absorber
%               p       Power of nonlinear absorber
%
% Copyright (c) 2018-2020, Kevin Dekemele (kevindekemele@gmail.com)
% Source code available at 
% Licensed under GNU GPLv3 license.
%

b = nchoosek(p,(p-1)/2);
%% maximum and minumum SIM
ZnaM = nthroot(2^(p-2)*((1-kappa)*(p+1)-sqrt((p-1)^2*(1-kappa)^2-4*p*xi_na^2))/(b*p),(p-1)/2);
ZnaP = nthroot(2^(p-2)*((1-kappa)*(p+1)+sqrt((p-1)^2*(1-kappa)^2-4*p*xi_na^2))/(b*p),(p-1)/2);

Z0P = (xi_na^2+(1-kappa-b/(2^(p-1))*ZnaM^((p-1)/2))^2).*ZnaM;
Z0M = (xi_na^2+(1-kappa-b/(2^(p-1))*ZnaP^((p-1)/2))^2).*ZnaP;


%% Find Zna(0) that corresponds with Z0(0)

r = [];
r(1)=b^2/2^(2*p-2);
r((p+1)/2)=-2*(1-kappa)*b/(2^(p-1));
r(p)=xi_na^2+(1-kappa)^2;
r(p+1)=-Z00;

rot = roots(r);
s_double =  rot(imag(rot)==0);
Zna0 = max(s_double) % Right most solution on SIM

%% Calculate the t for each state
t1 = p*b^2/((p-1)*2^(2*(p-1)))*Zna0^(p-1)-b*((1-kappa)*(p+1))/((p-1)/2*2^(p-1))*Zna0^((p-1)/2) + ((1-kappa)^2+xi_na^2)*log(Zna0);
t2 = p*b^2/((p-1)*2^(2*(p-1)))*ZnaP^(p-1)-b*((1-kappa)*(p+1))/((p-1)/2*2^(p-1))*ZnaP^((p-1)/2) + ((1-kappa)^2+xi_na^2)*log(ZnaP);

T_pump = (t1-t2)/(xi_na);
end

