function [k_x, k_y] = demod_miniMRI_fun(S, t, fc)
%Dado un nombre de archivo y una frecuencia central, devuelve el espacio k.

%M = readtable(filename);

%S = table2array(M(:,2))-mean(table2array(M(:,2)));
%t = table2array(M(:,1)) - 10.11 * 10^(-3);

S = S -mean(S);
t = t - 10.11 * 10^(-3);

T = t(2) - t(1);
% Define filter specifications
fs = 1/T;  % Sampling frequency (adjust as needed)
Fpass1 = 425000;  % First passband frequency (Hz)
Fpass2 = 429492;  % Second passband frequency (Hz)
Fstop1 = 426025;  % First stopband frequency (Hz)
Fstop2 = 428467;  % Second stopband frequency (Hz)
Apass = 1;  % Passband ripple (dB)
Astop = 40; % Stopband attenuation in dB
%d = designfilt('bandstopiir', 'FilterOrder', 10, ...
%    'StopbandFrequency1', Fstop1, 'StopbandFrequency2', Fstop2, ...
%    'PassbandFrequency1', Fpass1, 'PassbandFrequency2', Fpass2, ...
%    'StopbandAttenuation', Astop, 'PassbandRipple', Apass, ...
%    'SampleRate', fs);
d = designfilt('bandstopiir', 'FilterOrder', 10, 'PassbandFrequency1', ...
               Fpass1, 'PassbandFrequency2', Fpass2, 'PassbandRipple', ...
               Apass, 'SampleRate', fs*1);

S_f = filter(d, S);

% NUT = 1

N = length(t);
U = 1/(N*T);
u = (0:N-1)*U;

% FT of signal
fftc = @(x) fft(fftshift(x));
Sk = fftc(S_f);
%Filter definition and carrier frequency (fc)
f = .4;

W_half_kHz = 13;
Fpass = W_half_kHz * 1e3;
Fstop = W_half_kHz * 1.01e3;

fs = 1/T;
filtertype = 'IIR';
Rp = 0.1;
Astop = 40;
LPF = dsp.LowpassFilter('SampleRate',fs,...
                             'FilterType',filtertype,...
                             'PassbandFrequency',Fpass,...
                             'StopbandFrequency',Fstop,...
                             'PassbandRipple',Rp,...
                             'StopbandAttenuation',Astop);
%Demodulated signal
ReS = LPF(cos(2*pi*fc*t) .* S_f);
ImS = LPF(sin(2*pi*fc*t) .* S_f);
S_lp = ReS-1i*ImS;

a = (-8192/2 + 1):1:(8192/2 - 1);
S_lp_c = conv(S_lp, sinc(a/5));
S_lp_c = S_lp_c(4096:16384-4097);
S_lp_ds = downsample(S_lp,64);

G = 25*1e-3;
gamma = 42576384;
TE = t(end-length(t)/2);
deltak = (T*gamma*G);

% Define the notch frequency and Q factor
notchFrequency = 429077;  % Frequency in Hz
Q = 100;  % Adjust this value to control the filter's bandwidth (higher Q for narrower filter)





k_spacevec = linspace(-(length(S_lp)/2)*deltak, (length(S_lp)/2)*deltak, length(S_lp));

k_y = fftshift(fft(S_lp_ds));
k_x = k_spacevec;
