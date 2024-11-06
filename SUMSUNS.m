% Parameters
fs = 1024;               
N = 2048;                
timeVec = (0:N-1) / fs;        
A1 = 10; f1 = 100; phi1 = 0;
A2 = 5; f2 = 200; phi2 = pi/6;
A3 = 2.5; f3 = 300; phi3 = pi/4;

% Sum of three sinusoids
signal = A1 * cos(2*pi*f1*timeVec + phi1) + ...
         A2 * cos(2*pi*f2*timeVec + phi2) + ...
         A3 * cos(2*pi*f3*timeVec + phi3);


disp('Generated signal:');
disp(signal);


f_cutoff1 = (f1 + f2) / 2;  % Midpoint between f1 and f2
[b1, a1] = butter(4, f_cutoff1 / (fs/2), 'low');

f_cutoff2 = [(f1 + f2) / 2, (f2 + f3) / 2] / (fs/2);
[b2, a2] = butter(4, f_cutoff2, 'bandpass');

f_cutoff3 = (f2 + f3) / 2;
[b3, a3] = butter(4, f_cutoff3 / (fs/2), 'high');

% Apply filters to the signal
filtered_signal1 = filtfilt(b1, a1, signal);
filtered_signal2 = filtfilt(b2, a2, signal);
filtered_signal3 = filtfilt(b3, a3, signal);

% periodograms
figure;
subplot(2, 2, 1);
periodogram(signal, [], [], fs);
title('Periodogram of Original Signal');

subplot(2, 2, 2);
periodogram(filtered_signal1, [], [], fs);
title('Periodogram of Filtered Signal 1 (Low-pass)');

subplot(2, 2, 3);
periodogram(filtered_signal2, [], [], fs);
title('Periodogram of Filtered Signal 2 (Band-pass)');

subplot(2, 2, 4);
periodogram(filtered_signal3, [], [], fs);
title('Periodogram of Filtered Signal 3 (High-pass)');


