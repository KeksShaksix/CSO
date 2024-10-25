% Параметры сигнала
A = 1;         % Амплитуда
f = 50;        % Частота (Гц)
phi = pi/4;    % Начальная фаза (рад)
Fs = 1000;     % Частота дискретизации (Гц)
N = 10^5;      % Количество отсчетов

% Временные отсчеты
t = (0:N-1)/Fs; 

% Создание сигнала
x = A * cos(2 * pi * f * t + phi);

% Визуализация сигнала
figure;
plot(t, x);
title('Гармонический сигнал');
xlabel('Время (с)');
ylabel('Амплитуда');
grid on;

