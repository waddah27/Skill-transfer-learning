function emg1_env_norm=DataProc(iemg,fsample,lsample,Fs)
T1=1; 
fc_smooth=3; % cutoff frequency Butterworth filter 
%FS: sampling frequency 

%% ECG LP filter design
Fpass1 = 0.1;         % Passband Frequency
% Change the stopband frequency from 3 to 5
Fstop1 = 5;          % Stopband Frequency
Apass1 = 1;           % Passband Ripple (dB)
Astop1 = 80;          % Stopband Attenuation (dB)
match = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY1 method.
h1  = fdesign.lowpass(Fpass1, Fstop1, Apass1, Astop1, Fs);
Hd1 = design(h1, 'cheby1', 'MatchExactly', match);
%% EMG processing 
figure,plot(iemg)
title('RawEMG') 
emg1=iemg;
% raw EMG signal recorded with MIE processed with Simpon's commands 
plot(emg1)
title('EMG records')
xlabel('Time samples')
ylabel('EMG')
fc=[50:50:Fs/2-1]; 
b=notch_filter(fc,T1,Fs); 
emg1=filtfilt(b,1,emg1); % notch filtering
emg1_detr=detrend(emg1)/std(emg1);% detrend
%length(emg1_detr)

emg1_rect=abs(emg1_detr); %  rectified


[pxx,f]=pwelch(emg1_rect,[],[],[],Fs);
pwelch(emg1_rect,[],[],[],Fs);

%Finding the top significant peaks of the signal in spectral domain
NumPks=3; %Number of peaks to be found
[pks,lox]=findpeaks(pxx,'NPeaks',NumPks,'SortStr','descend');
hold on, plot(f(lox)/1000,10*log10(pks),'or'),hold off
components=sort([f(lox)])

% Design of LP cheby filter
N     = 10;   % Order
Fpass = components(NumPks);  % Passband Frequency
Apass = 1;    % Passband Ripple (dB)

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.lowpass('N,Fp,Ap', N, Fpass, Apass, Fs);
Hd = design(h, 'cheby1');
emg1_rect=filter(Hd,emg1_rect);
%length(emg1_rect)
emg1_rect=medfilt1(emg1_rect,20);%Remove spikes
%length(emg1_rect)
[b_smooth,a_smooth]=butter(2,fc_smooth/Fs*2); % Butterworth Filter
emg1_env=filtfilt(b_smooth,a_smooth,emg1_rect); %Extensors, linear envelope 
%length(emg1_env)  
emg1_env_norm=emg1_env(fsample:lsample)/max(emg1_env(fsample:lsample)); 
figure, plot(emg1_env_norm)
title('EMG processed')
xlabel('Time samples')
ylabel('EMG')

end
