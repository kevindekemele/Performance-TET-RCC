close all
clear all
% The numerical simulations as performed in the research paper
% "Performance Measures for Targeted Energy Transfer and Resonance
% Capture Cascading in Nonlinear Energy Sinks" by Kevin Dekemele et al,
% section 7: Multi-modal vibrations and resonance capture cascading.
% DOI: 10.1007/s11071-018-4190-5

%% Add functions
 addpath('../../');
%% MDOF system

% Main system
m = 1;
k = 1;
M = [m 0 0; 0 m 0; 0 0 m];
K = [2*k -k 0;-k 2*k -k; 0 -k k];

%Eigenfrequences and vectors
[E,F] = eig(K,M);
EFREQ = sqrt(F);

%initial conditions and initial modal conditions
x0_dot = [0.1;0;0];
P0_dot = inv(E)*x0_dot;

% Absorber
mna = 0.06


% Mode to damp
j = 3;
% Attachment point of NES
l = 2;
n = length(E);

%Here xi = 0.1
cna = 0.1*EFREQ(j,j)*mna;
mshape = E(l,j);


%modal mass
m1=1/(E(l,j)^2);

c1=0*EFREQ(j,j)*m1;

%decomposing initial conditions
u0 = E(l,j)*P0_dot(j);

ep = mna/m1;
xi = c1/m1/EFREQ(j,j)/ep;
xi_na = (cna/mna)/(EFREQ(j,j));

%% Tuning

% Tuning parameters, uneven power and linear part
p = 3;
kappa = 0;

b = nchoosek(p,(p-1)/2);

ZnaM = nthroot(2^(p-2)*((1-kappa)*(p+1)-sqrt((p-1)^2*(1-kappa)^2-4*p*xi_na^2))/(b*p),(p-1)/2)
Z0P = (xi_na^2+(1-kappa-b/(2^(p-1))*ZnaM^((p-1)/2))^2).*ZnaM;
%Coefficient nonlinear part
kna = 1.1*Z0P^((p-1)/2)*mna*EFREQ(j,j)^(p+1)/(u0^(p-1)); % Normal absorber

%Omega
Om = kna/mna/EFREQ(j,j)^(p+1);

klin = kappa*mna*EFREQ(j,j)^2;



%% Sim of each mode seperataly


% Simulation parameters
u0_real= 1*u0; % Actual applied initial conditions, to verify robustness
Tl =2000;


for i=1:n
    
k1=EFREQ(i,i)^2/(E(l,i)^2);
m1=1/(E(l,i)^2);
c1=0*EFREQ(j,j)*m1;

u0 = E(l,i)*P0_dot(i);
% Dimensionless parameters
ep = mna/m1;
xi = c1/m1/EFREQ(i,i)/ep;
xi_na = (cna/mna)/(EFREQ(i,i));
kappa = klin/mna/EFREQ(i,i)^2;


ZnaM = nthroot(2^(p-2)*((p+1)-sqrt((p-1)^2-4*p*xi_na^2))/(b*p),(p-1)/2);
ZnaP = nthroot(2^(p-2)*((p+1)+sqrt((p-1)^2-4*p*xi_na^2))/(b*p),(p-1)/2);

Z0P = (xi_na^2+(1-b/(2^(p-1))*ZnaM^((p-1)/2))^2).*ZnaM;
Z0M = (xi_na^2+(1-b/(2^(p-1))*ZnaP^((p-1)/2))^2).*ZnaP;

u0 = E(l,i)*P0_dot(i);
u0_real = u0;
Om = kna/mna/EFREQ(i,i)^(p+1);

%% Single mode Simulations
Zoo = nthroot(Om*u0_real^(p-1),(p-1)/2);
u0 = [0 0];
u0_dot = [u0_real 0];
fs = 100;
t=0:1/fs:Tl;
 f1 = @(t,y)[y(3);y(4);...
    -(k1*y(1) + c1*y(3) + kna*(y(1)-y(2))^p + klin*(y(1)-y(2))+cna*(y(3)-y(4)))/m1;...
    -(kna*(y(2)-y(1))^p + klin*(y(2)-y(1))+ cna*(y(4)-y(3)))/mna;...
    ];
    
Prec = 1e-12;
options = odeset('RelTol',Prec,'AbsTol',[Prec Prec Prec Prec]);
[T1,Y1] = ode45(f1,[0 Tl],[u0 u0_dot],options);

figure
plot(T1,Y1(:,1))
ylabel("x")
xlabel("time [s]")
title("Mode " + i)
figure
plot(T1,Y1(:,2) - Y1(:,1))
ylabel("x_{na}-x")
xlabel("time [s]")
title("Mode " + i)
%E fast
[peaks_en0,locs_en0] = findpeaks(((EFREQ(i,i)^2)*(Y1(:,1)+mna/m1*Y1(:,2)).^2+(Y1(:,3)+mna/m1*Y1(:,4)).^2),T1);
[peaks_ena,locs_ena] = findpeaks(((EFREQ(i,i)^2)*(Y1(:,1)-Y1(:,2)).^2+(Y1(:,3)-Y1(:,4)).^2),T1);

figure
subplot(2,1,1)
plot(locs_en0,peaks_en0,'k','LineWidth',2);
hold on
subplot(2,1,2)
plot(locs_ena,peaks_ena,'k','LineWidth',2);
hold on
    
[Z0_slow, Zna_slow,T_slow]=slowManifold(xi,xi_na,Zoo,0,p,Tl*ep*(EFREQ(i,i)));
      
subplot(2,1,1)
plot(T_slow/ep/((EFREQ(i,i))),Z0_slow/((Om)^(2/(p-1))),'k--','LineWidth',2);
hold on
legend('E_{0}fast','E_{0}')
subplot(2,1,2)
plot(T_slow/ep/((EFREQ(i,i))),Zna_slow/((Om)^(2/(p-1))),'k--','LineWidth',2);
legend('E_{na}fast','E_{na}')
hold on
title("Comparison slow flow and real of single mode " + i)

%% Performance metrics
T_pump_mode(i) = pumpingTime(xi_na,Zoo,kappa,p)/ep/(EFREQ(i,i))*2*pi
E_TET_mode(i) = energyDissTET(xi_na,Zoo,kappa,p)


end
%Cascading time
T_cascade = sum(T_pump_mode)

%% Actual numerical simulations
f1 = @(t,y)[y(5);y(6);y(7);y(8);...
    -(k*y(1)  + k*(y(1)-y(2)) )/m;...
    -(k*(y(2)-y(1)) +klin*(y(2)-y(4)) + k*(y(2)-y(3))+kna*(y(2)-y(4))^p +cna*(y(6)-y(8)))/m;...
    -(k*(y(3)-y(2)))/m;...
    -(kna*(y(4)-y(2))^p +cna*(y(8)-y(6))+klin*(y(4)-y(2)))/mna;
    ];

x0 = [0 0 0 0]; 
x0_dot = [0.1 0 0 0];
t = 0:0.01:Tl;
options = odeset('RelTol',Prec,'AbsTol',[Prec Prec Prec Prec Prec Prec Prec Prec]);
[T4,Y4] = ode45(f1,t,[x0 x0_dot],options);
figure
plot(T4,Y4(:,4)-Y4(:,2))
ylabel("x_{na}-x_{2}")
title("NES vibration")
figure
plot(T4,Y4(:,1),T4,Y4(:,2),T4,Y4(:,3))
ylabel("x_{i}")
legend("x_{1}","x_{2}","x_{3}")
ylabel("Time [s]")
title("Main system's vibration")

   