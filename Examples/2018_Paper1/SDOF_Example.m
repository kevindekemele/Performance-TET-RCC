close all
clear all

% The numerical simulations as performed in the research paper
% "Performance Measures for Targeted Energy Transfer and Resonance
% Capture Cascading in Nonlinear Energy Sinks" by Kevin Dekemele et al,
% section 6: Numerical study and comparison actual dynamics
% with slow flow. DOI: 10.1007/s11071-018-4190-5

%% Add functions
 addpath('../../');


%% SDOF system

% Main system
m = 1;
k = 1;
c = 0;
omega0 = sqrt(k/m);
x0_dot = 0.1;

% Absorber

mna =0.02;
cna = 0.1*omega0*mna;
ep = mna/m;
xi = c/m/omega0/ep;
xi_na = (cna/mna)/omega0 ;


% Tuning parameters, these were changed for the other simulations
p = 3;
kappa = 0;



%% Tuning
b = nchoosek(p,(p-1)/2);

ZnaM = nthroot(2^(p-2)*((1-kappa)*(p+1)-sqrt((p-1)^2*(1-kappa)^2-4*p*xi_na^2))/(b*p),(p-1)/2)
Z0P = (xi_na^2+(1-kappa-b/(2^(p-1))*ZnaM^((p-1)/2))^2).*ZnaM;

kna = 1*Z0P^((p-1)/2)*mna*omega0^(p+1)/(x0_dot^(p-1))*1.1; % Normal absorber
Om = kna/mna/omega0^(p+1);

klin = kappa*mna*omega0^2;


%% Simulations

% Simulation parameters
% x0_real is the actual applied initial condition, change to test
% robustness
x0_real= x0_dot;
Tl =200;

Zoo = nthroot(Om*x0_real^(p-1),(p-1)/2)
x0 = [0 0];
x0_dot = [x0_real 0];
fs = 100;
a=0;
t=0:1/fs:Tl;
 f1 = @(t,y)[y(3);y(4);...
    -(k*y(1) + c*y(3) + a*k*y(1)^3 + kna*(y(1)-y(2))^p+klin*(y(1)-y(2))+cna*(y(3)-y(4)))/m;...
    -(kna*(y(2)-y(1))^p + klin*(y(2) - y(1))+ cna*(y(4)-y(3)))/mna;...
    ];
    
Prec = 1e-12;
% Actual numerical simulation

options = odeset('RelTol',Prec,'AbsTol',[Prec Prec Prec Prec]);
      [T1,Y1] = ode45(f1,[0 Tl],[x0 x0_dot],options);
      
% Plot simulation result      
figure(1)
plot(T1,Y1(:,1))
ylabel('x [m]')
figure(2)
plot(T1,Y1(:,2) - Y1(:,1))
ylabel('x_{na}-x [m]')
xlabel('Time [s]')
% The E_Fast as defined in equation (35)
eno = ((omega0^2)*(Y1(:,1)+ep*Y1(:,2)).^2+(Y1(:,3)+ep*Y1(:,4)).^2);
ena = ((omega0^2)*(Y1(:,1)-Y1(:,2)).^2+(Y1(:,3)-Y1(:,4)).^2);
[peaks_en0,locs_en0] = findpeaks(((omega0^2)*(Y1(:,1)+ep*Y1(:,2)).^2+(Y1(:,3)+ep*Y1(:,4)).^2),T1);
[peaks_ena,locs_ena] = findpeaks(((omega0^2)*(Y1(:,1)-Y1(:,2)).^2+(Y1(:,3)-Y1(:,4)).^2),T1);
      
      
      
     
figure(3)
subplot(2,1,1)
plot(locs_en0,peaks_en0,'k','LineWidth',2);
hold on
subplot(2,1,2)
plot(locs_ena,peaks_ena,'k','LineWidth',2);
hold on
% Slow flow simulation
[Z0_slow, Zna_slow,T_slow]=slowManifold(xi,xi_na,Zoo,kappa,p,Tl*ep,false);

%Plot slow from and fast on eachother
subplot(2,1,1)
plot(T_slow/ep,Z0_slow/((Om)^(2/(p-1))),'k--','LineWidth',2);
legend('E0 fast','E0 slow')
subplot(2,1,2)
plot(T_slow/ep,Zna_slow/((Om)^(2/(p-1))),'k--','LineWidth',2);
legend('Ena fast','Ena slow')
xlabel('Time [s]')
figure
plot(Zna_slow,Z0_slow,(Om)^(2/(p-1))*((omega0^2)*(Y1(:,1)-Y1(:,2)).^2+(Y1(:,3)-Y1(:,4)).^2),((omega0^2)*(Y1(:,1)+ep*Y1(:,2)).^2+(Y1(:,3)+ep*Y1(:,4)).^2)*(Om)^(2/(p-1)))
xlabel('Zna')
ylabel('Z0')
%     

%% Performance metrics

T_pump_mode_kappa= pumpingTime(xi_na,Zoo,kappa,p)/ep/(omega0)*2*pi
E_TET_kappa = energyDissTET(xi_na,Zoo,kappa,p)

   