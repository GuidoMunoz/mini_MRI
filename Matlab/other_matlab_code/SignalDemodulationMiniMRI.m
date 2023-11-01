%Importing Signal
close all

M = readtable('TE16ms_TR1s.csv');

%S = output(110:end);
%t = t_2(110:end) - 10.01e-3;

S = table2array(M(100:end,2))-mean(table2array(M(:,2)));
t = table2array(M(100:end,1)) - 10.11 * 10^(-3);

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
fc = 409000;
W_half_kHz = 100;
Fpass = W_half_kHz * 1e3;
Fstop = W_half_kHz * 1.01e3;
%Acquired Signal
figure(1);
subplot 311
plot(t*1e3, S*1e3); xlabel("Time [ms]"); ylabel("Voltage (mv)");
axis([t(1)*1e3 t(end-100)*1e3 -8 8])
title("Acquired Signal")

%FT of Acquired Signal
subplot 412
plot(u*1e-6, abs(Sk)); xlabel("Frequency [MHz]")
 hold on
 plot([fc fc]*1e-6, [0 10], 'r--')
 plot([fc+Fpass fc+Fpass]*1e-6, [0 10], 'b--')
 plot([fc-Fpass fc-Fpass]*1e-6, [0 10], 'b--')
 plot([fc+Fstop fc+Fstop]*1e-6, [0 10], 'g--')
 plot([fc-Fstop fc-Fstop]*1e-6, [0 10], 'g--')
 hold off
 axis([fc*1e-6-0.3 fc*1e-6+0.3 0 5])
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

%Moving sinc
a = (-8192/2 + 1):1:(8192/2);
S_lp = conv(S_lp, sinc(a/50));

S_lp = S_lp(4096:16384-4097);
S_lp = downsample(S_lp, 64);
%S_lp = movingSinc(S_lp);
%Plot of demodulated signal

G = 13*1e-3;
gamma = 42576384;
T = (t(65) - t(1));
%TE = t(end-length(t)/2);
deltak = (gamma*G*T);

% Moving average window size
windowSize = 64;

% Calculate the number of output values
outputSize = floor(numel(S_lp) / windowSize);

outputVector = zeros(outputSize, 1);

% Reshape the input vector into a matrix with windowSize columns
for i = 1:outputSize
    if i == 1
        outputVector(i) = mean(S_lp(1:(i*windowSize+windowSize/2)));
    elseif i == outputSize
        outputVector(i) = mean(S_lp((i*windowSize-windowSize/2):end));
    else
        outputVector(i) = mean(S_lp((i*windowSize-windowSize/2):(i*windowSize+windowSize/2)));
    end
end
% Compute the moving average

reduceFactor = 64;
% Calculate the number of output values
outputSize2 = floor(numel(S) / reduceFactor);

% Initialize the output vector
outputVector2 = zeros(outputSize2, 1);

% Copy every reduceFactor-th value from the input to the output vector
for i = 1:outputSize2
    outputVector2(i) = t((i-1)*reduceFactor + 1);
end

for i = 1
    
end

k_spacevec = linspace(-(length(S_lp)/2)*deltak, (length(S_lp)/2)*deltak, length(S_lp));
figure(2)
subplot 211
plot(k_spacevec, real(S_lp), 'b-')
hold on
plot(k_spacevec, imag(S_lp), 'r-')
plot(k_spacevec, abs(S_lp), 'k-')
%plot((outputVector2-TE)*(64/FOV), -abs(outputVector),'ko-')
hold off
xlabel("k-space units [cycles/m]")
ylabel("Amplitude")
legend("ReS","ImS","|S|")
title("k-space signal (Demodulated Signal)")
axis([-200 200 -0.1 0.1])
%FT of demodulated signal

FOV = 1/deltak;
centered_vec_u = (u-(N-1)*U/2);

space_vec = linspace(-FOV/2, FOV/2, length(S_lp));
subplot 212
plot(space_vec*100, abs(fftshift(fftc(S_lp))), 'o-')
hold on

%plot([0.2,0.2],[0,120],'r')
%5plot([0.5,0.5],[0,120],'r')
%plot([-0.2,-0.2],[0,120],'r')
%plot([-0.5,-0.5],[0,120],'r')
%plot([-0.2,0.2],[0,0],'r')
%plot([0.5,1],[0,0],'r')
%plot([-0.5,-1],[0,0],'r')
%plot([-0.5,-0.2],[120,120],'r')
%plot([0.5,0.2],[120,120],'r')
axis([-3 3 0 3])
xlabel("Space [cm]")
title("Acquired image (FT of Demodulated Signal)")

espacio_K_y = outputVector;
espacio_K_x = k_spacevec;

imagen_x = abs(fftshift(fftc(outputVector)));
imagen_y = space_vec;
%%
f = figure(3);

t_graph = 1:1:11.4e3;

rf_pulse = zeros(1,length(t_graph));
inicio = 0;
fin = 200;

inicio_2 = 4000;
fin_2 = 4200;
for i = 1:length(t_graph)
    if t_graph(i) > inicio
        if t_graph(i) < fin
            rf_pulse(i) = 90;
        end
    end
    if t_graph(i) > inicio_2
        if t_graph(i) < fin_2
            rf_pulse(i) = 180;
        end
    end
end

Gx = zeros(1,length(t_graph));
inicio_slope = 300;
inicio = 400;
fin = 2800;
fin_slope = 2900;
inicio_slope2 = 5530;
inicio2 = 5730;
fin2 = 10530;
fin_slope2 = 10730;

for i = 1:length(t_graph)
    if t_graph(i) < inicio_slope
        Gx(i) = 0;
    elseif t_graph(i) < inicio
        Gx(i) = 25/(inicio - inicio_slope)*t_graph(i) - 25*inicio_slope/(inicio - inicio_slope);
    elseif t_graph(i) < fin
        Gx(i) = 25;
    elseif t_graph(i) < fin_slope
        Gx(i) = -25/(fin_slope - fin)*t_graph(i) + 25*fin_slope/(fin_slope - fin);
    elseif t_graph(i) < inicio_slope2
        Gx(i) = 0;
    elseif t_graph(i) < inicio2
        Gx(i) = 25/(inicio2 - inicio_slope2)*t_graph(i) - 25*inicio_slope2/(inicio2 - inicio_slope2);
    elseif t_graph(i) < fin2
        Gx(i) = 25;
    elseif t_graph(i) < fin_slope2
        Gx(i) = -25/(fin_slope2 - fin2)*t_graph(i) + 25*fin_slope2/(fin_slope2 - fin2);
    else
        Gx(i) = 0;
    end
end

AD = zeros(1,length(t_graph));
begin = 7330;
stop = 8930;
for i = 1:length(t_graph)
    if t_graph(i) > begin
        if t_graph(i) < stop
            AD(i) = 1;
        end
    end
end

subplot 311
plot((t_graph'-100)/1000, rf_pulse','LineWidth',2,'Color','r')
title('RF Pulse')
xlabel('Time (ms)')
ylabel('Angle (Degree)')
xlim([-0.1 11.3])
ylim([0 180])
yticks([0 90 180])
subplot 312
plot((t_graph'-100)/1000, Gx','LineWidth',2,'Color','b')
title('X Gradient')
xlabel('Time (ms)')
ylabel('Gradient Strength (mT/m)')
xlim([-0.1 11.3])
subplot 313
plot((t_graph'-100)/1000, AD','LineWidth',2,'Color','g')
title('Analog-to-Digital Converter')
xlabel('Time (ms)')
xlim([-0.1 11.3])
f.Position = [100 100 700 700];