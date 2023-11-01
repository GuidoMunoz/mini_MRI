% Define the pulse parameters
T = 0.2;          % Pulse length in milliseconds
original_TE = 8; % Original echo time in milliseconds
total_duration = original_TE/2 + T/2; % Total duration of the signal
desired_samples = 32000; % Desired number of samples
frequency = 4.6e6; % Frequency in Hz (4.6 MHz)

% Calculate the phase between both pulses
phase_difference_degrees = 90; % Desired phase difference in degrees
phase_difference_radians = deg2rad(phase_difference_degrees);

% Calculate the time shift required to achieve the desired phase difference
time_shift_seconds = (phase_difference_radians * 1e-3) / (2 * pi * frequency); % Convert milliseconds to seconds

% Calculate the adjusted TE based on the time shift
new_TE = original_TE + 2 * time_shift_seconds;

% Calculate the time vector with the desired number of samples
t_desired = linspace(0, total_duration, desired_samples);

% Create the envelope for the first pulse (90-degree pulse)
amplitude_pulse1 = 0.5; % Amplitude of the first pulse
pulse1_center = T / 2; % Center time of the first pulse
pulse1_envelope = amplitude_pulse1 * sinc((t_desired - pulse1_center) / (T/10));
pulse1_envelope(t_desired < 0 | t_desired > T) = 0; % Zero out values outside pulse duration

% Calculate the adjusted center time of the second pulse based on TE and phase
pulse2_center = (new_TE/2); % Adjusted center time of the second pulse
pulse2_envelope = amplitude_pulse2 * sinc((t_desired - pulse2_center) / (T/10));
pulse2_envelope(t_desired < (new_TE/2 - T/2) | t_desired > (new_TE/2 + T/2)) = 0; % Zero out values outside pulse duration

% Combine both pulse envelopes into a single envelope vector
combined_envelope = pulse1_envelope + pulse2_envelope;

% Generate the sine wave with the specified frequency and phase
sine_wave = sin(2 * pi * frequency * t_desired);

% Calculate the phase of the combined signal
combined_signal = combined_envelope .* sine_wave;
phase_signal = angle(sum(combined_signal));

% Calculate the number of samples to shift
samples_to_shift = round(time_shift_seconds / (t_desired(2) - t_desired(1)));

% Shift the signal horizontally by the calculated number of samples
shifted_signal = circshift(combined_signal, samples_to_shift);

% Plot the signal
figure;
subplot(2, 1, 1);
plot(t_desired, abs(combined_envelope), 'b', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Amplitude');
title('RF Pulse Envelopes (Sinc Pulses)');
grid on;

subplot(2, 1, 2);
plot(t_desired, real(shifted_signal), 'r', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Amplitude');
title('Signal after Phase Shift');
grid on;

% Display the new TE value and the calculated phase of the combined signal
fprintf('New TE value: %.2f ms\n', new_TE);
fprintf('Calculated phase of the combined signal: %.2f degrees\n', rad2deg(phase_signal));

file_name = sprintf('spinecho_TE%.2fms_freqs%.0fkhz.csv', new_TE, 1 / (t_desired(2) - t_desired(1)));

% Write the combined_envelope vector to a CSV file
writematrix(combined_envelope, file_name);
