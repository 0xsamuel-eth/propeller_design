%% Sam Masten, Mickey Kerrigan, Hanna McDaniel
%  MAE451
%  Lab 5 Code

clear;
clc;

%% Analytical Calculations

%% Given Values

N = 2; %number of blades
D = 10; %in
R = D/2; %in
c = 1; %in

%Chosen Values
dr = 0.05;
RPM = 3200;
f = RPM / 60;

%Nondimensional radius
r_R = 0.1:dr:1;

%Radius increments
r = 0.5:0.25:5;

%Theta Values
theta = [45 41 38 35 33 31 29 27 25 23 21 19 17 15 13 11 9 7 5]; %deg

%% Importing CL and CD

foil = readmatrix('xf-naca2412-il-200000.xlsx');
foil_alpha = foil(:,1);
foil_CL = foil(:,2);
foil_CD = foil(:,3);

%Polyfit CL and CD to find the values at r increments
CL_fit = polyfit(foil_alpha,foil_CL,3);
CD_fit = polyfit(foil_alpha,foil_CD,3);

%% Calculations

% Calculate Propeller Solidity
sigma = N*c/(pi*R);

%Calculate omega
omega = 2 * pi * f;

%Calculate V2
V2 = omega .* r;

%For Each J Value
J = [0 0.6661 0.9381 1.1476];

for i = 1:length(J)

    % Calculate V0
    V0{i} = J(i) * f * D;

    %Calculate Angle Phi
    phi{i} = atand(V0{i} ./ V2);

    %Calculate Angle Alpha
    alpha{i} = theta - phi{i};

    %Fit the polyfit coefficients to 3rd order polynomial
    CL{i} = CL_fit(1).*alpha{i}.^3 + CL_fit(2).*alpha{i}.^2 + CL_fit(3).*alpha{i} + CL_fit(4);
    CD{i} = CD_fit(1).*alpha{i}.^3 + CD_fit(2).*alpha{i}.^2 + CD_fit(3).*alpha{i} + CD_fit(4);

    %Calculate dCp
    dCT{i} = (J(i)^2 + (pi^2) .* (r.^2)) * sigma .* (CL{i} .* cosd(phi{i})- CD{i} .* sind(phi{i}))*dr;
    sumCT{i} = sum(dCT{i})*dr;

    %Calculate dCP
    dCP{i} = (pi .* r) .* (J(i)^2 + (pi^2) .* (r.^2)) .* sigma .* (CL{i} .* sind(phi{i}) + CD{i} .* cosd(phi{i}))*dr;
    sumCP{i} = sum(dCP{i})*dr;

    %Calculate eta
    sum_eta{i} = sumCT{i} * J(i) / sumCP{i};

end

%Convert values back to vectors from cell arrays
sumCT = cell2mat(sumCT);
sumCP = cell2mat(sumCP);
sum_eta = cell2mat(sum_eta);



%% Experimental Calculations

%% Constants

R = 287; %ideal gas constant
RPM = 3200; %rotations per minute

%% Import Data

Data = readmatrix('Lab5Data.xlsx');

%% Tabling Values

diameter = Data(1,1); %diameter in in
q = Data(:,3); %psf
Temp = Data(:,4); %degrees F
p = Data(1,5); %psi
f = Data(1,6); %Hz
T = Data(:,7); %lb
Q = Data(:,8); %in-lb

%% Converting Values

diameter = diameter * 0.0254; %in to m
q = q.* 47.88; %psf to Pa
Temp = (Temp-32).*(5/9) + 273.15; %deg F to K;
p = p * 6894.76; %psi to Pa
T = T.* 4.4482216; %lb to N
Q = Q.* 0.11298482933333; %in-lb to N-m

%% Calculating rho

rho = p ./ (R*Temp);

%% Calculations

%V_inf
Vinf = real(sqrt(2*q./rho));

%J
J = Vinf / (f*diameter);

%C_T
CT = T ./ (rho*(f^2)*(diameter^4));

%C_Q
CQ = Q ./ (rho*(f^2)*(diameter^5));

%C_P
CP = 2 * pi * CQ;

%eta_pr
eta = CT .* J ./ CP;

%% Plots

% CT vs. J
figure(4)

subplot(3,1,1)
plot(J,sumCT, 'r');
hold on;
plot(J,CT, 'b');
hold on;
scatter(J,sumCT, 'r','filled');
hold on;
scatter(J,CT, 'b','filled');
xlabel('Advance Ratio');
ylabel('Thrust Coefficient');
title('C_{T} vs J');
grid on;
legend('Analytical','Experimental')

% CP vs. J
subplot(3,1,2)
plot(J,sumCP, 'r');
hold on;
plot(J,CP, 'b');
hold on;
scatter(J,sumCP, 'r','filled');
hold on;
scatter(J,CP, 'b','filled');
xlabel('Advance Ratio');
ylabel('Power Coefficient');
title('C_{P} vs J');
grid on;
legend('Analytical','Experimental')

% eta_pr vs. J
subplot(3,1,3)
plot(J,sum_eta, 'r');
hold on;
plot(J,eta, 'b');
hold on;
scatter(J,sum_eta, 'r','filled');
hold on;
scatter(J,eta, 'b','filled');
xlabel('Advance Ratio');
ylabel('Propeller Efficiency');
title('{\eta}_{pr} vs J');
grid on;
legend('Analytical','Experimental')
