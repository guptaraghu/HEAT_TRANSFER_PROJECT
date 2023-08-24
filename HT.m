clear all;
clc;
ti = 0;
tf = 0.10;%final thickness of TPS
dx = 0.010; % space step
gp = floor((tf - ti) / dx) + 1; % No of grid points: discretization
dt = 1; %time step
n = floor(720 / dt) + 1; %No of grid points in time space
K = 1.5;                %thermal_conductivity
den = 1300;             % density
c = 820;                %specific heat
alpha = K / (den * c);  %thermal diffusivity
cfl = alpha * dt / dx^2; %cfl no
x = linspace(ti, tf, gp); %partiation the thickness
t = linspace(0, 720, n);  % partiation of the time
sigma = 5.67e-8; % steafan Boltzmann constant
e = 0.9; % Emissivity
hf = zeros(1, n);%heat flux intital all value zero
for p = 1:201
    hf(p) = 4000 * t(p) + 200000;
end
for q = 202:551
    hf(q) = 1000000;
end
for r = 552:651
    hf(r) = 4960000 - 7200 * t(r);
end
for s = 652:n
    hf(s) = 2.8e5;
end

Temp = zeros(gp, n);
Temp(:, 1) = 300; % Initial condition
for j = 1:n-1
    for i = 2:gp-1
        Temp(i, j+1) = cfl*Temp(i+1, j) + (1-2*cfl)*Temp(i, j) + cfl*Temp(i-1, j);
    end
    Temp(gp, j+1) = Temp(gp-1, j+1); % Boundary condition
    Temp(1, j+1) = Temp(1, j) + dt / (den * c) * (hf(j) + K/dx* (Temp(2, j) - Temp(1, j)) - e * sigma * Temp(1, j)^4); % Radiative cooling
end
figure();
hold on 
xlabel('Distance (m)');
ylabel('Temperature (K)');
plot(x, Temp(:, 100), 'b'); 
plot(x, Temp(:, 200),'g' ); 
plot(x, Temp(:, 300),'r'); 
plot(x, Temp(:, 450),'y'); 
plot(x, Temp(:, 550),'k'); 
plot(x, Temp(:, 650),'m'); 
plot(x, Temp(:, end));
legend('100s','200s','300s','450s','550s','650s','720') 
title(" temperature distribution through the TPS material")
hold off 
figure(); 
plot(Temp(1,:)); 
title(" temperature variation with time at the exposed end")
xlabel('Time(sec)'); 
ylabel('Temperature(k)');
figure()
plot(Temp(11,:));
title(" temperature variation with time at the back end")
xlabel('Time(sec)'); 
ylabel('Temperature(k)');

