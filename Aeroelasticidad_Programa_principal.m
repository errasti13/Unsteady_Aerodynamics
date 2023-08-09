% Aerodynamics Simulation Code

% Author: Jon Errasti Odriozola
% Date: May 2020

% Description: This MATLAB code performs aerodynamics simulations using
% various models and notations. It covers steady and non-steady cases,
% gusts, and different airfoil characteristics. The code also visualizes
% the results using plots and legends.

%%
% Clear the workspace and command window
clear all
clc

%% Input (Menu)

wingRootChord = 2;      % Wing root chord (m)
angleOfAttack = 1*pi/180; % Angle of attack (rad)
flightSpeed = 1;        % Flight speed (m/s)
numProfileDivisions = 10; % Number of divisions in the airfoil profile
flightHeight = 0;       % Flight height (m)
Qinf = flightSpeed * ones(numProfileDivisions, 1) * -angleOfAttack;
w = 1;

%% Air Parameters

ambientTemperature = 288.15 - 0.0065 * flightHeight; % Ambient temperature (K)
ambientPressure = 101325 * (ambientTemperature / 288.15)^(9.81 / (287.075 * 0.0065)); % Ambient pressure (Pa)
airDensity = ambientPressure / (287.074 * ambientTemperature); % Air density (kg/m^3)

%% NACA 2412 Airfoil

camberLine = 2 / 100;
locationOfMaximumCamber = 4 / 10;
thickness = 12 / 100;

%% Time Definition

initialTime = 0;
finalTime = 50 * wingRootChord * flightSpeed;
timeIncrement = 0.01;
numTimeSteps = (finalTime - initialTime) / timeIncrement;

%% Variables

iteration = 1;
alpha = 0;
height = 0.5;
Integral_r = angleOfAttack;

%% Discretize the airfoil profile

for i = 1:numProfileDivisions
    xPanel(i) = wingRootChord / numProfileDivisions * (i - 1);
    xpControl(i) = wingRootChord / numProfileDivisions * (i - 1 + 3/4);
    xCgamma(i) = wingRootChord / numProfileDivisions * (i - 1 + 1/4);
end

%% Discretize the wake

for i = 1:numTimeSteps
    xWake(i) = wingRootChord + timeIncrement * (i - 1);
end

for i = 1:numTimeSteps
    xCgammaWake(i) = wingRootChord + timeIncrement * 0.2 + timeIncrement * (i - 1);
end




%% Calculations

A = zeros(numProfileDivisions,numProfileDivisions);
for i = 1:numProfileDivisions
    for j = 1:numProfileDivisions
        A(i, j) = -1 / (2 * pi) * (1 / (xpControl(i) - xCgamma(j)));
    end
end

VectorCirculation = A \ Qinf;

DistrLe = airDensity * flightSpeed * VectorCirculation';
DistrPe = DistrLe / (wingRootChord / numProfileDivisions);

Le = sum(DistrLe);
Cle = Le / (0.5 * airDensity * flightSpeed^2 * wingRootChord);
Cle1 = sum(VectorCirculation);
MomLe = -sum(DistrLe .* xCgamma);

CMe = MomLe / (0.5 * airDensity * flightSpeed^2 * wingRootChord^2);

Cl_e = ones(1, numTimeSteps) * Cle;


%% Non-steady calculation

X0 = -flightSpeed * timeIncrement;
Xp0 = -flightSpeed;
theta = sin(w * timeIncrement);
thetap = w * cos(w * timeIncrement);

%% Procedure

timeVector = zeros(1, numTimeSteps);

%% Angle of Attack

% Angle of attack for non-steady case

alphaMatrix = zeros(numProfileDivisions, numTimeSteps);

for i = 1:numTimeSteps
    if i == 1
        timeVector(1, i) = timeIncrement;
    else
        timeVector(1, i) = timeVector(1, i - 1) + timeIncrement;
    end
    alphaMatrix(:, i) = -angleOfAttack;
end

% Wagner

timeVectorAn = linspace(0, numTimeSteps * timeIncrement, numTimeSteps);

Phi = zeros(1, numTimeSteps);
Cl_alpha = zeros(1, numTimeSteps);

for i = 1:length(timeVectorAn)
    Phi(1, i) = 1 - 0.165 * exp(-0.045 * timeVectorAn(1, i)) - 0.335 * exp(-0.3 * timeVectorAn(1, i));
    Cl_alpha(1, i) = 2 * pi * angleOfAttack * Phi(1, i);
end
%% Flexion

% Flexion for non-steady case

flexionMatrix = zeros(numProfileDivisions, numTimeSteps);
timeVector = zeros(1, numTimeSteps);

for i = 1:numTimeSteps
    if i == 1
        timeVector(1, i) = timeIncrement;
    else
        timeVector(1, i) = timeVector(1, i - 1) + timeIncrement;
    end
    flexionMatrix(:, i) = height * 1i * iteration * exp(1i * iteration * timeVector(1, i));
end
% Theodorsen

timeVectorAn = linspace(0, numTimeSteps * timeIncrement, numTimeSteps);
C = Funcion_Theodorsen(iteration);

Claf1 = zeros(1, numTimeSteps);
Claf2 = zeros(1, numTimeSteps);
Claf = zeros(1, numTimeSteps);
alpha = zeros(1, numTimeSteps);

for i = 1:length(timeVectorAn)
    Claf1(1, i) = -2 * pi * height * C * 1i * iteration * exp(1i * iteration * timeVectorAn(1, i));
    Claf2(1, i) = pi * height * iteration^2 * exp(1i * iteration * timeVectorAn(1, i));
    Claf(1, i) = Claf1(1, i) + Claf2(1, i);
    alpha(1, i) = -1i * iteration * height * exp(1i * iteration * timeVectorAn(1, i));
end


%% Step Gust

% Gust for non-steady case

gustMatrix = zeros(numProfileDivisions, numTimeSteps);
timeVector = zeros(1, numTimeSteps);
v_Int_r = ones(numProfileDivisions, 1) * -Integral_r;

for i = 1:numTimeSteps
    if i == 1
        timeVector(1, i) = timeIncrement;
    else
        timeVector(1, i) = timeVector(1, i - 1) + timeIncrement;
    end
    for j = 1:numProfileDivisions
        if timeVector(1, i) - xpControl(1, j) - 1 == 0
            gustMatrix(j, i) = -Int_r;
        else
            gustMatrix(j, i) = -Integral_r * heaviside(timeVector(1, i) - xpControl(1, j) - 1);
        end
    end
end

% Kussner

timeVectorAn = linspace(0, numTimeSteps * timeIncrement, numTimeSteps);
Psi = zeros(1, numTimeSteps);
Clar = zeros(1, numTimeSteps);

for i = 1:length(timeVectorAn)
    Psi(1, i) = 1 - 0.5 * exp(-0.13 * timeVectorAn(1, i)) - 0.5 * exp(-timeVectorAn(1, i));
    Clar(1, i) = 2 * pi * Integral_r * Psi(1, i);
end


[Cl_ane, Circulacion_ane, Circulacion_gradiente_ane, Circulacion_estela_ane] = CalculoNoEstacionario(numProfileDivisions, numTimeSteps, A, alphaMatrix, xpControl, xCgammaWake, timeIncrement, xCgammaWake);
[Cl_fne, Circulacion_fne, Circulacion_gradiente_fne, Circulacion_estela_fne] = CalculoNoEstacionario(numProfileDivisions, numTimeSteps, A, flexionMatrix, xpControl, xCgammaWake, timeIncrement, xCgammaWake);
[Cl_rne, Circulacion_rne, Circulacion_gradiente_rne, Circulacion_estela_rne] = CalculoNoEstacionario(numProfileDivisions, numTimeSteps, A, gustMatrix, xpControl, xCgammaWake, timeIncrement, xCgammaWake);
%% Section 2

% Initial Conditions

initialConditions = [1, 0, 0, 0];

alpha0 = 1 * pi / 180;
h0 = 0;
dhe = 0;
dalpha = 0;
dt = timeIncrement;

% Defined Parameters

wingSpan = 1;
ah = -0.5;
x_alpha = 0.555;
r_alfa = sqrt(0.8222);
walpha = 1;
dif_w = 1.83;  % wh/walpha
wh = dif_w * walpha;
nu = 172;
m = nu * pi * airDensity * wingSpan^2;
Ialpha = r_alfa * m * wingSpan^2;
kh = m * wh^2;
kalpha = Ialpha * walpha^2;

xG = 1/2 * wingRootChord + wingSpan * ah + wingSpan * x_alpha;


% Fuerzas y momentos no estacionarias

L_alpha_ne=0;

M_alpha_ne=0;

%for i=1:length(vector_tiempo)

%funcion=@(t,x)[x(2);(1/(r_alfa^2-x_alpha^2))*(M_alpha_ne/m+-x_alpha*L_alpha_ne/m+x_alpha*wh^2*x(3)-r_alfa^2*walpha^2*x(1));
%x(4);(1/(x_alpha^2-r_alfa^2))*(x_alpha*M_alpha_ne/m+r_alfa^2*L_alpha_ne/m-r_alfa^2*wh^2*x(3)-r_alfa^2*walpha^2*x_alpha*x(1))];

%[t,x]=ode45(funcion,[0,t_final],x0);

%alphapp=1/(Ialpha-m*x_alpha^2)*(M_alpha_ne+kh*h0*x_alpha-kalpha*alpha0);

%dalphanp1=dalpha+dt*alphapp;

%dhenp1=dhe+dt*1/m*(-L_alpha_ne-kh*h0-m*x_alpha*alphapp);

%henp1=dhe*dt+h0;

%alphanp1=dalpha*dt+alpha0;

%v_alphaprima=ones(Np,Ne)*-alpha0(1); 

%[Cl_anep,g_anep,gd_anep,gw_anep]=CalculoNoEstacionario(Np,Ne,A,v_alphaprima,xpcontrol,xcgamma_estela,inc_t,xcgamma);

%Distr_L_alpha_ne=1/2*ro*U_inf^2*c*Cl_anep;

%L_alpha_ne=-sum(Distr_L_alpha_ne);

%M_alpha_ne=sum(Distr_L_alpha_ne.*(x_alpha+ah));

%h0=henp1; dhe=dhenp1; alpha0=alphanp1; dalpha=dalphanp1;

%hev(i)=h0; alphav(i)=alpha0; L_alpha(i)=L_alpha_ne;
%dhev(i)=dhe; dalphav(i)=dalpha; M_alpha(i)=M_alpha_ne;

%end

%% Apartado 3- Metodo VG

k=1:-0.000001:0.000001;

for i=1:length(k)
    kred=k(i);
C=Funcion_Theodorsen(kred);
Clh=(kred^2-2*kred*1i*C);
Cla=(-kred*1i-ah*kred^2-2*C*(1+(1/2-ah)*kred*1i));    
CMah=(-ah*kred^2+2*C*(ah+1/2)*kred*1i);
CMaa=((-1/2+ah)*kred*1i+(1/8-ah^2)*kred^2+2*C*(ah+1/2)*(1+(1/2-ah)*kred*1i));
    
% Matriz rigidez
m11=1;
m12=x_alpha;
m21=m12;
m22=r_alfa^2;

M=[m11,m12;m21,m22];

% Matriz de disipacion

k11=dif_w^2;
k12=0;
k21=k12;
k22=r_alfa^2;

K=[k11,k12;k21,k22];

% Matriz de fuerzas


q11=Clh;
q12=Cla;
q21=CMah;
q22=CMaa;

Q=[q11,q12;q21,q22];

Matriz_A=(M+(1/(nu*kred^2))*Q);
Matriz_h_alpha=inv(K)*Matriz_A;

lambda=eig(Matriz_h_alpha);

for j=1:length(lambda)

g=imag(lambda(j)/real(lambda(j)));

if (real(lambda(j))>0) && (abs(g)<0.00001)
    wflameo=walpha/sqrt(real(lambda(j)));
    Uflameo=wflameo/kred
    break
else
    continue
end
end
end

%% Representacion grafica  

figure(1)

for i=1:length(xPanel)
plot([0,wingRootChord],0)
hold on
plot([0,wingRootChord],[0,0],'b-')
plot(xPanel(i),0,'-*')
plot(xpControl(i),0,'-o')
plot(xCgamma(i),0,'-+')
xlim([0 wingRootChord]) 
ylim([-1 1])
title('	Perfil discretizado')
xlabel('x'),ylabel('y')
grid on
end

figure(2)
plot(timeVector, Cl_e, 'b')
hold on
plot(timeVectorAn, Cl_alpha, 'r')
plot(timeVectorAn, Cl_ane, 'g')
title('Wagner')
legend({'Steady Cl', 'Wagner', 'Non-steady Cl'}, 'Location', 'Southeast')
xlabel('Time'), ylabel('Cl')
grid on


figure(3)
plot(timeVector, Cl_e, 'b')
hold on
plot(timeVector, Clar, 'r')
plot(timeVectorAn, Cl_rne, 'g')
title('Kussner')
legend({'Steady Cl', 'Kussner', 'Non-steady Cl'}, 'Location', 'Southeast')
xlabel('Time'), ylabel('Cl')
grid on


figure(4)
plot(timeVectorAn, Claf, 'b')
hold on
plot(timeVectorAn, Cl_fne, 'r')
title('Theodorsen')
legend({'Theodorsen', 'Non-steady Cl'}, 'Location', 'Southeast')
xlabel('Time'), ylabel('Cl')
grid on

