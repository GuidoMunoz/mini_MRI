%Importing Signal
M = readtable('TE8ms_TR6s.csv');
%S = output;
%t = t - 10.01*1e-3;
%M = readtable('acq3_04_1.5V.csv');
S = table2array(M(:,2)) - mean(table2array(M(:,2)));
t = table2array(M(:,1));
t = t - t(1);
%M = acq6(:,1);
%t = acq6(:,2):
%a = 3;

%S = (S).*exp(-t*1e4);

%t = table2array(M(:,1));


% NUT = 1
T = t(2) - t(1);
N = length(t);
U = 1/(N*T);
u = (0:N-1)*U;
% FT of signal
fftc = @(x) fft(fftshift(x));
Sk = fftc(S);
%Filter definition and carrier frequency (fc)
f = .4;
fc = 397000;
W_half_kHz = 50;
Fpass = W_half_kHz * 1e3;
Fstop = W_half_kHz * 1.01e3;
%Acquired Signal
figure(1);
subplot(2,1,1)
plot(t*1e6, S); xlabel("Time [us]")

axis([t(1)*1e6 t(end)*1e6 -0.02 0.02])
%axis([7 250 0.035 0.075])
%axis([7 250 0.180 0.220])

title("Acquired Signal")
%FT of Acquired Signal
subplot(2,1,2)
plot(u*1e-6, abs(Sk)); xlabel("Frequency [MHz]")
hold on
%plot([fc fc]*1e-6, [0 10], 'r--')
%plot([fc+Fpass fc+Fpass]*1e-6, [0 10], 'b--')
%plot([fc-Fpass fc-Fpass]*1e-6, [0 10], 'b--')
%plot([fc+Fstop fc+Fstop]*1e-6, [0 10], 'g--')
%plot([fc-Fstop fc-Fstop]*1e-6, [0 10], 'g--')
hold off
axis([fc*1e-6-0.3 fc*1e-6+0.3 0 15])
%axis([fc*1e-6-0.3 fc*1e-6+0.3 0 5])
title("FT of Acquired Signal")
%IIR Filter design
fs = 1/T;
filtertype = 'IIR';
Rp = 0.1;
Astop = 80;
LPF = dsp.LowpassFilter('SampleRate',fs,...
                             'FilterType',filtertype,...
                             'PassbandFrequency',Fpass,...
                             'StopbandFrequency',Fstop,...
                             'PassbandRipple',Rp,...
                             'StopbandAttenuation',Astop);
%Demodulated signal
ReS = LPF(cos(2*pi*fc*t) .* S);
ImS = LPF(sin(2*pi*fc*t) .* S);
S_lp = ReS-1i*ImS;
%Plot of demodulated signal
figure(2)
subplot(2,1,1)
plot(t*1e6, real(S_lp))
hold on
plot(t*1e6, imag(S_lp))
plot(t*1e6, abs(S_lp), 'k')
plot(t*1e6, -abs(S_lp), 'k')
hold off
xlabel("Time [us]")
legend("ReS","ImS","|S|")
title("Demodulated Signal")
axis([t(1)*1e6 t(end)*1e6 -10e-3 10e-3])
%FT of demodulated signal
subplot(2,1,2)
plot((u-(N-1)*U/2)*1e-3, abs(fftshift(fftc(S_lp))))
axis([-2*Fstop*1e-3 2*Fstop*1e-3 0 20])
xlabel("Frequency [kHz]")
title("FT of Demodulated Signal")