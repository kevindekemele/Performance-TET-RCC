function [Z0,Zna,T] = slowManifold(xi,xi_na,Z00,kappa,p,T_end,right)

% 
% slowManifold returns the evolution on the slow time manifold of Zna and  Z0
% as defined in the research paper
% "Performance Measures for Targeted Energy Transfer and Resonance
% Capture Cascading in Nonlinear Energy Sinks" by Kevin Dekemele et al,
% equation (21). DOI: 10.1007/s11071-018-4190-5
%
%       [ Z0,Zna,T] = slowManifold(xi,xi_na,Z00,kappa,p,T_end)
%
%               Z0      Evolution of main system dimensionless energy
%               Zna     Evolution of absorber dimensionless energy
%               T       Time vector in slow time scale [epsilon*s]
%
%               xi      Dimensionless damping of main system (physical/modal)
%               xi_na   Dimensionless damping of absorber (physical)
%               Z00     Initial dimensionless 'energy'
%               kappa   Dimensionless linear stifness of absorber
%               p       Power of nonlinear spring
%               T_end   The final time of the simulations, mutiplied by the mode vibration frequency
%                       and epsilon
%               right   'true' if the right branch of the SIM should always
%                        be chosen
%
% Copyright (c) 2018, Kevin Dekemele (kevindekemele@gmail.com)
% Source code available at 
% Licensed under GNU GPLv3 license.
%

% Binomial Coefficent
b = nchoosek(p,(p-1)/2);
%% maximum and minumum SIM

ZnaM = nthroot(2^(p-2)*((1-kappa)*(p+1)-sqrt((p-1)^2*(1-kappa)^2-4*p*xi_na^2))/(b*p),(p-1)/2);
ZnaP = nthroot(2^(p-2)*((1-kappa)*(p+1)+sqrt((p-1)^2*(1-kappa)^2-4*p*xi_na^2))/(b*p),(p-1)/2);

Z0P = (xi_na^2+(1-kappa-b/(2^(p-1))*ZnaM^((p-1)/2))^2).*ZnaM;
Z0M = (xi_na^2+(1-kappa-b/(2^(p-1))*ZnaP^((p-1)/2))^2).*ZnaP;


%% Find Zna(0) that corresponds with Z0(0)

syms x;
s = solve((xi_na^2+(1-kappa-b/(2^(p-1))*x^((p-1)/2))^2)*x-Z00==0,x);
s_double = double(vpa(s));
cond = abs(real(s_double)) > 1e6 * abs(imag(s_double)); % Condition to Remove imaginary round-offs errors
s_double(cond) = real(s_double(cond));

%% Choose Zna(0) corresponding with Z0(0)
if(Z00 < Z0P && ~right)
    Zna0 = min(s_double);
else
    Zna0 = max(s_double);
end

%% Find Zna after jump in SIM

s = solve((xi_na^2+(1-kappa-b/(2^(p-1))*x^((p-1)/2))^2)*x-Z0M==0,x);
s_double = double(vpa(s));
cond = abs(real(s_double)) > 1e6 * abs(imag(s_double)); % Condition to Remove imaginary round-offs errors
s_double(cond) = real(s_double(cond));
Znajump = min(s_double);


%% Simulate up to jump
f2 = @(t,y)[-xi_na*y(2)-xi*y(1);...
(-xi_na*y(2)-xi*y(1))/((p*b^2)/(2^(2*(p-1)))*y(2)^(p-1)-(b*(1-kappa)*(p+1))/(2^(p-1))*y(2)^((p-1)/2)+(1-kappa)^2+(xi_na)^2)];

Prec = 1e-16;
options = odeset('RelTol',Prec,'AbsTol',[Prec]);
 [T2,Y2] = ode45(f2,[0 T_end],[Z00 Zna0],options);
     
if(Z00 >Z0P || right)
    
    Zna0=0.999*Znajump;
    Y2(end,1);
    %Z00 after jump
    Z00=(xi_na^2+(1-kappa-b/(2^(p-1))*Zna0^((p-1)/2))^2)*Zna0;
%% Simulate After jump

    options = odeset('RelTol',Prec,'AbsTol',[Prec]);
    [T3,Y3] = ode45(f2,[0 T_end-T2(end)],[Z00 Zna0],options);

    T = [T2;T2(end)+T3];
    Z0 = [Y2(:,1);Y3(:,1)];
    Zna = [Y2(:,2);Y3(:,2)];


else
     T = [T2];
    Z0 = [Y2(:,1)];
    Zna = [Y2(:,2)];

end

