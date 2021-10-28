global Fom Lom lst A Dlmt MA gamma1 gamma2 d lt u Mmes phi_0 time e delay
%%THIS SCRIPT ESTIMATES WITH EXISTENCE OF GROUND TRUTH DATA%%

%% Load data for a movement
Data='FW_SRL_S1.mat';
%load(Data);
Samp_freq=s1.fs;%10240;
ix=find(s1.Data(:,7)==4.1);
iemg=s1.Data(ix,3);
force=s1.Data(ix,12);

%% Data Processing
fsample=500; % first sample to be included 
lsample=length(iemg); % last sample to be included 
sample_int=lsample-fsample+1; % sample interval to be processed
filtx11=DataProc(iemg,fsample,lsample,Samp_freq);
figure, plot(filtx11)
title('Processed EMG')
%% NMS MODEL
time=sample_int;
t=1:1:time;
gamma1=0.5;
gamma2=0.5;
delay=2;
d=ceil(delay*Samp_freq/1000);
%initial Muscle Parameters 
Fom=12; %N
Lom=6.3; %cm
lst=24.4;%
A=-3;
Dlmt=4.1; %cm;
Lmt=lst+1.2*Lom;
MA=1.5; %cm
lb=[1, 1,-3, 0.9*MA,0,0]; %lower bounds
ub=[2000, 30, -0.001,1.1*MA,1,1]; %Upper bounds

%Musc_model(Params)
D2R=2*pi/360;
R2D =180/pi;
phi_0=3*D2R; %rad
% r_ave=1.9; %cm
lambda=0.15;

%% Model parameters calibration
% Pass fixed parameters to objfun
objfun = @(x)NMS(x,filtx11,force,time,Samp_freq,phi_0);

% Set nondefault solver options
%options2 = optimoptions('fmincon','PlotFcn','optimplotfvalconstr');
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% Solve
    x0=(lb+ub)/ceil((10*rand(1)));
    J = NMS(x0,filtx11, force, time,Samp_freq, phi_0);
    [solution2,objectiveValue2] = fmincon(objfun,x0,[],[],[],[],lb,ub,[],options)
    J1 = NMS(solution2,filtx11, force, time,Samp_freq, phi_0);
    Enhancement = (J-J1)/J*100
% J=zeros(1,5);
% J1=100*ones(1,5);
% Best_sol=zeros(1,6);
% %Best_val=zeros(1,6);
% for i=2:1:3
%     x0=(lb+ub)/ceil((10*rand(1)));
%     J(i) = NMS(x0,filtx11, force, time,Samp_freq, phi_0);
%     [solution2,objectiveValue2] = fmincon(objfun,x0,[],[],[],[],lb,ub,[],options)
%     J1(i) = NMS(solution2,filtx11, force, time,Samp_freq, phi_0);
%     if J1(i)<J1(i-1)
%         Best_val=J1(i);
%         ind=i;
%         for j=1:1:length(solution2)
%            Best_sol(j)=solution2(j);
%         end
%         
%     end
%     Enhancement = (J-J1)/J*100
% end



