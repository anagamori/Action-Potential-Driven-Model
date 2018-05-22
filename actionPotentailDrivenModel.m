%--------------------------------------------------------------------------
% actionPotentialDrivenModel.m
% Author: Akira Nagamori
% Last update: 5/20/2018
% Implementation of the action potential-drive model by Kim et al. (2015)
%--------------------------------------------------------------------------

close all
clear all
clc

%--------------------------------------------------------------------------
Fs = 40000;
%--------------------------------------------------------------------------
% Module 1 parameters (Table 1 )
K1 = 3000; % [M-1*s-1]
K2 = 3; % [s-1]
K3 = 400; % [M-1*s-1]
K4 = 1; % [s-1]
K5 = 4e5; % [M-1*s-1]
K6 = 150; % [s-1]
K = 850; % [M-1] : site binding constant for Ca2+ to activate the pump in SR
R_max = 10; % [s-1]: maximum permeability of the SR
U_max = 2000; % [M*s-1]: maximum pump rate
tau_1 = 3/1000; % [ms]
tau_2 = 25/1000; % [ms]
phi_1 = 0.03; % [mm-1]
phi_2 = 1.23; 
phi_3 = 0.01; % [mm-1]
phi_4 = 1.08; % 
%--------------------------------------------------------------------------
% Initial condition
Ca_SR = 2.5/1000; % [mM]
CS_0 = 30/1000; % [mM]
B_0 = 0.43/1000; % [mM]
T_0 = 70*10^-3/1000; %[mM]

Ca_SRCS = 0;
Ca_SP = 0;
Ca_SPB = 0;
Ca_SPT = 0;


%--------------------------------------------------------------------------
% Module 2 parameters (Table 1 )
C1 = 0.128;
C2 = 0.093;
C3 = 61.206/1000; %[ms]
C4 = -13.116;
C5 = 5.095;
alpha = 2;
alpha_1 = 4.77;
alpha_2 = 400/1000; % [ms]
alpha_3 = 160/1000; % [ms]
beta = 0.47;
gamma = 0.001; 

A_tilda = 0;
%--------------------------------------------------------------------------
FR = 40;
time = 0:1/Fs:5;

diffusion1 = 1-exp(-time./tau_1);
diffusion2 = exp(-time./tau_2);

spikeTrain_temp = zeros(1,length(time));
spikeTrain_temp_2 = zeros(1,length(time));
spikeTrain = [zeros(1,1*Fs) spikeTrainGenerator(0:1/Fs:1,Fs,FR) zeros(1,length(time)-(2*Fs+1))];
% spikeTrain = zeros(1,length(time));
% spikeTrain(1*Fs) =1;
diffusion1_conv = zeros(1,length(time));
diffusion2_conv = zeros(1,length(time));
spike_time = 0;

R_vec = zeros(1,length(time));
U_vec = zeros(1,length(time));
Ca_SR_vec = zeros(1,length(time));
Ca_SRCS_vec = zeros(1,length(time));
Ca_SP_vec = zeros(1,length(time));
Ca_SP_dot_vec = zeros(1,length(time));
Ca_SPB_vec = zeros(1,length(time));
Ca_SPT_vec = zeros(1,length(time));
A_vec = zeros(1,length(time));

for t = 1:length(time)
    if spikeTrain(t) == 1
        spikeTrain_temp(t) = 1;
        spikeTrain_temp_2(t) = 1;
        diffusion1_conv = conv(spikeTrain_temp, diffusion1);
        diffusion2_conv = conv(spikeTrain_temp_2, diffusion2);
        
        spikeTrain_temp_2 = zeros(1,length(time));
    end
    R = Ca_SR*R_max*diffusion1_conv(t)*diffusion2_conv(t); % equation 3
    U = U_max*(Ca_SP^2*K^2/(1+Ca_SP*K+Ca_SP^2*K^2))^2; % equation 4
    
    Ca_SR_dot = -K1*CS_0*Ca_SR + (K1*Ca_SR+K2)*Ca_SRCS -R + U; % equation 1
    Ca_SRCS_dot = K1*CS_0*Ca_SR - (K1*Ca_SR+K2)*Ca_SRCS; % equation 2
    
    
    Ca_SR = Ca_SR_dot/Fs + Ca_SR;
    Ca_SRCS = Ca_SRCS_dot/Fs + Ca_SRCS;
    
    Ca_SP_dot = -(K3*B_0+K5*T_0)*Ca_SP + (K3*Ca_SP+K4)*Ca_SPB ...
        + (K5*Ca_SP+K6)*Ca_SPT + R - U; % equation 5
    Ca_SPB_dot = K3*B_0*Ca_SP - (K3*Ca_SP+K4)*Ca_SPB; % equation 6
    Ca_SPT_dot = K5*T_0*Ca_SP - (K5*Ca_SP+K6)*Ca_SPT; % equation 7
    
    Ca_SP = Ca_SP_dot/Fs + Ca_SP;
    Ca_SPB = Ca_SPB_dot/Fs + Ca_SPB;
    Ca_SPT = Ca_SPT_dot/Fs + Ca_SPT;
    
    Ca_SR_vec(t) = Ca_SR;
    Ca_SRCS_vec(t) = Ca_SRCS;
    Ca_SP_vec(t) = Ca_SP;
    Ca_SP_dot_vec(t) = Ca_SP_dot;
    Ca_SPB_vec(t) = Ca_SPB;
    Ca_SPT_vec(t) = Ca_SPT;
   
    A_tild_inf = 0.5*(1+tanh((Ca_SPT/T_0-C1)/C2));
    tau_A_tilda = C3*(cosh((Ca_SPT/T_0-C4)/(2*C5)))^-1;
    A_tilda_dot = (A_tild_inf - A_tilda)/tau_A_tilda;
    A_tilda = A_tilda_dot/Fs + A_tilda;
    A = A_tilda*alpha;
    A_vec(t) = A;
    
    R_vec(t) = R;
    U_vec(t) = U;
end

figure(1)
subplot(3,1,1)
plot(time,Ca_SR_vec)
subplot(3,1,2)
plot(time,Ca_SP_vec)
subplot(3,1,3)
plot(time,A_vec)

function spikeTrain = spikeTrainGenerator(t,Fs,freq)

spikeTrain = zeros(1,length(t));
ISI = round(1/freq*Fs);
numSpikes = round(length(t)/ISI);
index = [1:numSpikes]*ISI;
index(index>length(t)) = [];
spikeTrain(index) = 1;
spikeTrain(1) = 1;

end
