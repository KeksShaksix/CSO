% spectrum_analysis_full.m
% Полный анализ спектра сигнала: ДПФ, оценка ширины спектра, дополнение, и измерение скорости расчетов.

% Определение нового дискретного сигнала
Fs = 1000; % Частота дискретизации
t = 0:1/Fs:1; % Временной вектор
input_signal = sin(2*pi*50*t) + 0.5*sin(2*pi*120*t) + 0.3*sin(2*pi*200*t); % Новый сигнал

% 1. Расчет ДПФ
N = length(input_signal); % Длина сигнала
X = fft(input_signal); % ДПФ сигнала
f = (0:N-1)*(Fs/N); % Частотный вектор

% Модуль и фаза ДПФ
magnitude = abs(X);
phase = angle(X);

% Построение графиков модуля и фазы
figure;

% График модуля
subplot(2, 1, 1);
stem(f, magnitude, 'filled');
title('Модуль спектра (ДПФ)');
xlabel('Частота (Гц)');
ylabel('Амплитуда');
xlim([0 Fs/2]); % Ограничение по оси X

% График фазы
subplot(2, 1, 2);
stem(f, phase, 'filled');
title('Фаза спектра (ДПФ)');
xlabel('Частота (Гц)');
ylabel('Фаза (рад)');
xlim([0 Fs/2]); % Ограничение по оси X

% 2. Оценка ширины спектра сигнала
energy = sum(input_signal.^2); % Энергия исходного сигнала
disp(['Энергия сигнала: ', num2str(energy)]);

X_copy = X; % Копия ДПФ
Nmax = 0; % Начальное значение
cumulative_energy = 0; % Начальная сумма энергии
percentage_threshold = 0.9; % Порог для энергии

while true
    Nmax = Nmax + 1; % Увеличение индекса
    if Nmax > length(X_copy) / 2 % Проверка на предел
        break;
    end
    
    harmonic_energy = abs(X_copy(Nmax))^2; % Энергия текущей гармоники
    cumulative_energy = cumulative_energy + harmonic_energy; % Обновление суммы энергии
    
    % Проверка условия на 90% энергии
    if cumulative_energy / energy >= percentage_threshold
        break;
    end
end

disp(['Минимальное количество гармоник для 90% энергии: ', num2str(Nmax)]);

% 3. Дополнение сигнала нулями
input_signal_extended = [input_signal, zeros(1, N)]; % Копия с добавлением нулей
N_extended = length(input_signal_extended); % Новая длина сигнала

% Вычисление ДПФ для дополненного сигнала
X_extended = fft(input_signal_extended); % ДПФ дополненного сигнала
f_extended = (0:N_extended-1)*(Fs/N_extended); % Новый частотный вектор

% Модуль и фаза ДПФ дополненного сигнала
magnitude_extended = abs(X_extended);
phase_extended = angle(X_extended);

% Построение графиков модуля и фазы для дополненного сигнала
figure;

% График модуля
subplot(2, 1, 1);
stem(f_extended, magnitude_extended, 'filled');
title('Модуль спектра (Дополненный ДПФ)');
xlabel('Частота (Гц)');
ylabel('Амплитуда');
xlim([0 Fs]); % Ограничение по оси X

% График фазы
subplot(2, 1, 2);
stem(f_extended, phase_extended, 'filled');
title('Фаза спектра (Дополненный ДПФ)');
xlabel('Частота (Гц)');
ylabel('Фаза (рад)');
xlim([0 Fs]); % Ограничение по оси X

% 4. Измерение скорости расчетов при вычислении ДПФ по прямой формуле
% Создание матрицы преобразования
T = zeros(N); % Квадратная матрица
for k = 1:N
    for n = 1:N
        T(k, n) = exp(-2*pi*1i*(k-1)*(n-1)/N); % Заполнение матрицы
    end
end

% Измерение времени вычислений
tic;
X_direct = T * input_signal.'; % Прямое умножение
time_direct = toc;

disp(['Время вычисления ДПФ по прямой формуле: ', num2str(time_direct), ' секунд']);

% 5. Измерение скорости расчетов при вычислении ДПФ с использованием быстрого алгоритма
tic;
X_fft = fft(input_signal); % Быстрое преобразование
time_fft = toc;

disp(['Время вычисления ДПФ с использованием fft: ', num2str(time_fft), ' секунд']);
