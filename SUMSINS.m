% Parameters
fs = 1024;               
N = 2048;               
timeVec = (0:N-1)/fs;          
A1 = 10; f1 = 100; phi1 = 0;
A2 = 5; f2 = 200; phi2 = pi/6;
A3 = 2.5; f3 = 300; phi3 = pi/4;

% sum of three sinusoids
sigVec = A1 * cos(2*pi*f1*timeVec + phi1) + ...
         A2 * cos(2*pi*f2*timeVec + phi2) + ...
         A3 * cos(2*pi*f3*timeVec + phi3);

% Design filters
% Filter for Signal 1 (Low-pass filter)
f_cutoff1 = (f1 + f2) / 2; % Midpoint between f1 and f2
filter1 = fir1(128, f_cutoff1/(fs/2), 'low');

% Filter for Signal 2 (Band-pass filter)
f_cutoff2 = [(f1 + f2) / 2, (f2 + f3) / 2] / (fs/2);
filter2 = fir1(128, f_cutoff2, 'bandpass');

% Filter for Signal 3 (High-pass filter)
f_cutoff3 = (f2 + f3) / 2;
filter3 = fir1(128, f_cutoff3/(fs/2), 'high');

% Apply filters to the signal
filtered_signal1 = filter(filter1, 1, sigVec);
filtered_signal2 = filter(filter2, 1, sigVec);
filtered_signal3 = filter(filter3, 1, sigVec);

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

% Optional: Plot time domain signals if desired
figure;
subplot(4, 1, 1);
plot(timeVec, sigVec);
title('Original Signal in Time Domain');


subplot(4, 1, 2);
plot(timeVec, filtered_signal1);
title('Filtered Signal 1 in Time Domain');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(4, 1, 3);
plot(timeVec, filtered_signal2);
title('Filtered Signal 2 in Time Domain');


subplot(4, 1, 4);
plot(timeVec, filtered_signal3);
title('Filtered Signal 3 in Time Domain');
