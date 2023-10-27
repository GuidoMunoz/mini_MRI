% Número de mediciones
num_mediciones = 18;

% Ángulos de rotación en grados
angulos = 00:10:170;
fc = 412700;

% Generar las mediciones (reemplaza esta parte con tus datos de MRI)
% Aquí estoy generando mediciones de prueba con ruido gaussiano
mediciones = zeros(128, num_mediciones);
for i = 1:num_mediciones
    filename = string("rects_"+angulos(i) + "grads");
    addpath("E:/average_10grads/" + filename)
    filename = string("rects_"+angulos(i) + "grads0%03d.csv");
    if angulos(i) == 0
        filename = "rects_00grads0%03d.csv";
    end
    [output, t] = averageSignals(filename);
    [a ,mediciones(:, i)] = demod_miniMRI_fun(output, t, fc);
end
%%

%Arregla el orden de las adquisiciones
mediciones_ordenadas = [mediciones(:,11:18), mediciones(:,1:10)];
%a = (-8192/2 + 1):1:(8192/2 - 1);
%S_lp = conv(S_lp, sinc(a/50));

%mediciones_ordenadas_conv = zeros(256, 18);
%for j = 1:18
%    mediciones_ordenadas_conv(:,j) = conv(mediciones_ordenadas(:,j),sinc(a/1));
%end

% Reconstrucción mediante Backprojection
imagen_reconstruida = zeros(256);  % Tamaño de la imagen reconstruida (ajusta si es necesario)

% Matriz de coordenadas de píxeles en la imagen reconstruida
%[X, Y] = meshgrid(1:size(imagen_reconstruida, 2), 1:size(imagen_reconstruida, 1));
%centro_x = size(imagen_reconstruida, 2) / 2;
%centro_y = size(imagen_reconstruida, 1) / 2;

    
angulo_rad = deg2rad(angulos(i));
imagen_reconstruida = iradon(abs(mediciones_ordenadas), angulos);
% Interpolación lineal para agregar la contribución de la proyección a cada píxel
% En lugar de interpolar, acumulamos directamente la contribución de cada proyección

% Normalizar la imagen reconstruida
imagen_reconstruida = imagen_reconstruida / num_mediciones;

% Mostrar la imagen reconstruida
figure(3);
imagesc(abs(imagen_reconstruida));
colormap(gray);
% Puedes ajustar los parámetros de visualización si es necesario:
% figure;
% imshow(imagen_reconstruida, [min_value max_value]);
% colormap(gray);
% colorbar;
% axis on;
% title('Imagen reconstruida');
%%
T = t(128) - t(1);
N = length(t);
U = 1/(N*T);
u = (0:N-1)*U;

G = 15*1e-3;
gamma = 42576384;
deltak = (gamma*G*T);

figure(2)
FOV = 1/deltak;
centered_vec_u = (u-(N-1)*U/2);

space_vec = linspace(-FOV/2, FOV/2, length(mediciones_ordenadas(:,15)));
plot(space_vec*100, abs(mediciones_ordenadas(:,4)), 'o-')
hold on

plot([0.2,0.2],[0,10],'r')
plot([0.5,0.5],[0,10],'r')
plot([-0.2,-0.2],[0,10],'r')
plot([-0.5,-0.5],[0,10],'r')
plot([-0.2,0.2],[0,0],'r')
plot([0.5,1],[0,0],'r')
plot([-0.5,-1],[0,0],'r')
plot([-0.5,-0.2],[10,10],'r')
plot([0.5,0.2],[10,10],'r')
%axis([-5 5 0 10])
ylim([0 0.1])
xlabel("Space [cm]")
title("Acquired image (FT of Demodulated Signal)")
%%
B = imgaussfilt(imagen_reconstruida, 1.5);
figure(4)
imagesc(B);
colormap(gray);