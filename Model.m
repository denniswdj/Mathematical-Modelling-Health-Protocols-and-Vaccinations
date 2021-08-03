clc
clear all

% Definition of Initial Conditions at time zero. 
% target vaksinasi DKI Jakarta: 8815157
initial.S = 5251843;        % S = target vaksinasi - V - D - R - Q - C
initial.C = 29839;          % Rata2 kasus positif aktif dalam 7 hari 
initial.Q = 15884 ; 
initial.R = 789261; 
initial.D = 12209; 
initial.V = 2716121;

% Definition of model parameters
a = 0.218;                              % penerapan physical distancing & testing Covid-19
b = 0.0883 * 0.4885 * 0.3937 * 0.2391;  % masker, cuci tangan, googgle, sarung tangan
c = 0.9786;                             % sanitasi rutin
d = 0.7 ;                               % penentuan durasi maksimal karyawan di kantor
e = 0.88;                               % ventilasi udara yang baik
    
param.betaa =  a * b * c * d * e * 0.125;   % laju transmisi asimptomatik
param.alpha = 0.5;                          % proporsi pasien asimptomatik
param.betap = a * b * c * d * e * 0.250;    % laju transmisi presimptomatik
param.tau = 0.192;                          % laju perkembangan gejala asimptomatik menjadi simptomatik
param.gamma = 0.0473;                       % laju kesembuhan
param.delta = 1.33 * (10^(-3));             % laju kematian
param.kappa = 1.525 * (10^(-3));            % laju vaksinasi

% Definition of the simulation time
end_time = 200;

initial_values = [];
variable_names = fieldnames(initial);
for i = 1:length(variable_names)
    initial_values = [initial_values; initial.(variable_names{i})];
end

[t, y] = ode45(@(t, x) ode_system(t, x, param), [0 end_time], initial_values);

% prepare legend texts
legend_texts = cell(length(variable_names), 1);
for i = 1:length(variable_names)
   text = [variable_names{i} '(t)'] ;
   legend_texts{i} = text;
end

% plot the results
subplot(3,2,1)
plot(t, y(:,1), '--k', 'LineWidth',2);
xlabel('Waktu (hari)');
ylabel('Jumlah (orang)');
title('Susceptible (S)');
legend(legend_texts(1));

subplot(3,2,2)
plot(t, y(:,2), 'cyan', 'LineWidth',2);
xlabel('Waktu (hari)');
ylabel('Jumlah (orang)');
title('Silent Carrier (C)');
legend(legend_texts(2));

subplot(3,2,3)
plot(t, y(:,3), 'yellow', 'LineWidth',2);
xlabel('Waktu (hari)');
ylabel('Jumlah (orang)');
title('Quarantined (Q)');
legend(legend_texts(3));

subplot(3,2,4)
plot(t, y(:,4), 'green', 'LineWidth',2);
xlabel('Waktu (hari)');
ylabel('Jumlah (orang)');
title('Recovered (R)');
legend(legend_texts(4));

subplot(3,2,5)
plot(t, y(:,5), 'red', 'LineWidth',2);
xlabel('Waktu (hari)');
ylabel('Jumlah (orang)');
title('Death (D)');
legend(legend_texts(5));

subplot(3,2,6)
plot(t, y(:,6), 'LineWidth', 2);
xlabel('Waktu (hari)');
ylabel('Jumlah (orang)');
title('Vaccinated (V)');
legend(legend_texts(6));

% Definition of the ODE system
function deriv = ode_system(t,x,param)
    % Input:
    %   t: Time (not used in this example because there is no
    %           explicit time dependence)
    %   x: Vector of the current values of all variables in the same
    %           order as defined the initial values.
    %   param: Used to pass parameter values.   
    % Output:
    %   deriv: Column vector of derivatives, must be the same order as the input vector x.
    S = x(1);
    C = x(2);
    Q = x(3);
    R = x(4);
    D = x(5);
    V = x(6);
    
    dS = -param.betaa * param.alpha * C * S / (S + C + R + V) - param.betap * (1 - param.alpha) * C * S / (S + C + R + V) - param.kappa * S;
    dC = param.betaa * param.alpha * C * S / (S + C + R + V) + param.betap * (1 - param.alpha) * C * S / (S + C + R + V) - param.tau * (1 - param.alpha) * C;
    dQ = param.tau * (1 - param.alpha) * C - param.gamma * Q - param.delta * Q;
    dR = param.gamma * Q;
    dD = param.delta * Q;
    dV = param.kappa * S;
    
    deriv = [dS; dC; dQ; dR; dD; dV];    
end