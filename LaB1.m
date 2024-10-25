
matlab_script = ""
% Параметры
U1 = 1;
U2 = 2;
U3 = 3;
U4 = 4;
Fd = 1000;% Частота дискретизации
T = 1 / Fd;% Интервал дискретизации
T1 = 0.01;  % Время T1 (в секундах)         
T2 = 0.02;           % Время T2 (в секундах)
N = round(T2 * Fd);  % Количество отсчетов N
x = randn(1, N);     % Пример сигнала x (здесь используется случайный сигнал)

% Убедимся, что x - это строка (если не так, то транспонируем)
if size(x, 1) > 1
    x = x';
end

% Вектор времени
t = 0:T:T2;  % Время от 0 до T2 с шагом T

% Формирование сигнала
x = zeros(size(t));  % Создаем вектор сигналов нулями
for k = 1:length(t)
    if t(k) <= T1
        x(k) = U1 + (U2 - U1) * (t(k) / T1);  % Линейное увеличение от U1 до U2
    else
        x(k) = U2 + (U4 - U2) * ((t(k) - T1) / (T2 - T1));  % Линейное увеличение от U2 до U4
    end
end

% Вектор времени для первого фрагмента
t1 = 0:T:T1;

% Расчет параметров a1 и b1
a1 = (U2 - U1) / T1;
b1 = U1;

% Значения отсчетов для первого линейного фрагмента
u_t1 = a1 * t1 + b1;

% Отображение значений отсчетов
disp(u_t1);

% Вектор времени для второго фрагмента
t2 = T1:T:T2;

% Расчет параметров a2 и b2
a2 = (U4 - U2) / (T2 - T1);
b2 = U2 - a2 * T1;

% Значения отсчетов для второго линейного фрагмента
u_t2 = a2 * t2 + b2;

% Отображение значений отсчетов
disp(u_t2);

% Вектор времени для всего сигнала
t = 0:T:T2;

% Расчет параметров a1, b1 для первого фрагмента
a1 = (U2 - U1) / T1;
b1 = U1;

% Расчет отсчетов для первого линейного фрагмента
u_t1 = a1 * t(t <= T1) + b1;

% Расчет параметров a2, b2 для второго фрагмента
a2 = (U4 - U2) / (T2 - T1);
b2 = U2 - a2 * T1;

% Расчет отсчетов для второго линейного фрагмента
u_t2 = a2 * t(t > T1) + b2;

% Соединяем оба фрагмента
u_t = [u_t1, u_t2];

% Построение графика с помощью функции plot
figure;
subplot(2,1,1);  % Разделение на два графика в одной фигуре
plot(t, u_t, '-o');
xlabel('Время (с)');
ylabel('Амплитуда');
title('Дискретный сигнал с использованием plot');
grid on;

% Построение графика с помощью функции stem
subplot(2,1,2);
stem(t, u_t);
xlabel('Время (с)');
ylabel('Амплитуда');
title('Дискретный сигнал с использованием stem');
grid on;

% Количество отсчетов N
N = round(T2 * Fd);  % N = T2 * Fd

% Формирование вектора-столбца номеров отсчетов k
k = (0:N-1)';  % Вектор-столбец от 0 до N-1


% Задаем количество элементов M для вектора частот
M = N;  % Можно выбрать от 500 до 1000
w = linspace(-pi, pi, M);  % Вектор частот от -pi до pi

% Вычисление попарных произведений столбца k и строки w
K_W_matrix = k * w;

% Умножение матрицы на -j для получения комплексной экспоненты
complex_exp_matrix = -1i * K_W_matrix;

% Вектор комплексных экспонент
complex_exponentials = exp(complex_exp_matrix);

% Умножение вектора отсчетов сигнала x на матрицу комплексных экспонент
X = sum(x.' .* complex_exponentials, 1);  % x должен быть вектором-столбцом

% Построение амплитудного и фазового спектров
f = (w / (2*pi)) * Fd;  % Линейная частота
amplitude_spectrum = abs(X);
phase_spectrum = angle(X);

figure;
subplot(2, 1, 1);
plot(f, amplitude_spectrum);
title('Amplitude Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

subplot(2, 1, 2);
plot(f, phase_spectrum);
title('Phase Spectrum');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
grid on;

% Формирование вектора моментов времени для расчета восстановленного сигнала
time_step = T / 10;  % Шаг времени в 10 раз меньше исходного
t_min = -5 * T;
t_max = (N - 1 + 5) * T;
t = t_min:time_step:t_max;

% Создание вектора для восстановленного сигнала
s_reconstructed = zeros(size(t));

% Цикл по k для восстановления сигнала
for k = 0:N-1
    s_kT = x(k+1);  % Значение сигнала в момент времени kT
    sinc_term = sinc((t - k*T) / T);  % Sinc-функция
    s_reconstructed = s_reconstructed + s_kT * sinc_term;  % Восстановленный сигнал
end
% Интервал дискретизации
T = 1 / Fd;

% Вектор времени для первого фрагмента
t1 = 0:T:T1;

% Расчет параметров a1 и b1
a1 = (U2 - U1) / T1;
b1 = U1;

% Значения отсчетов для первого линейного фрагмента
u_t1 = a1 * t1 + b1;

% Отображение значений отсчетов
disp(u_t1);
% Построение графика восстановленного сигнала
figure;
plot(t, s_reconstructed, 'LineWidth', 1.5);
title('Reconstructed Signal using sinc Interpolation');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

% Сравнение дискретного и восстановленного сигналов
figure;
plot(t, s_reconstructed, 'LineWidth', 1.5, 'DisplayName', 'Reconstructed Signal');
hold on;
stem((0:N-1)*T, x, 'r', 'LineWidth', 1.5, 'DisplayName', 'Original Discrete Signal');
title('Comparison of Original Discrete and Reconstructed Signals');
xlabel('Time (seconds)');
ylabel('Amplitude');
legend('show');
grid on;
hold off;


