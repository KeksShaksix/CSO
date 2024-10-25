% filter_analysis.m
% Анализ фильтра: получение выходного сигнала, построение характеристик,
% и разложение функции передачи на простые дроби.

% 1. Определение коэффициентов фильтра
b = [0.2, 0.2, 0.2, 0.2, 0.2]; % Числитель (коэффициенты FIR-фильтра)
a = 1; % Знаменатель (для FIR-фильтра a=1)

% Создание входного сигнала
Fs = 1000; % Частота дискретизации
t = 0:1/Fs:1; % Временной вектор
input_signal = sin(2*pi*50*t) + 0.5*randn(size(t)); % Синусоидальный сигнал с шумом

% 2. Применение фильтра и получение выходного сигнала
output_signal = filter(b, a, input_signal);

% 3. Построение графиков характеристик фильтра
figure;

% АЧХ (Magnitude Response)
subplot(3, 2, 1);
freqz(b, a);
title('АЧХ (Magnitude Response)');

% ФЧХ (Phase Response)
subplot(3, 2, 2);
[H, w] = freqz(b, a);
plot(w/pi, unwrap(angle(H)));
title('ФЧХ (Phase Response)');
xlabel('Нормализованная частота (\times \pi рад/об)');
ylabel('Фаза (рад)');

% Групповая задержка (Group Delay Response)
subplot(3, 2, 3);
[groupDelay, w] = grpdelay(b, a);
plot(w/pi, groupDelay);
title('Групповая задержка (Group Delay Response)');
xlabel('Нормализованная частота (\times \pi рад/об)');
ylabel('Групповая задержка (отсчеты)');

% Импульсная характеристика (Impulse Response)
subplot(3, 2, 4);
impulse_response = filter(b, a, [1 zeros(1, 99)]); % Импульсный сигнал
stem(impulse_response);
title('Импульсная характеристика (Impulse Response)');
xlabel('Отсчеты');
ylabel('Амплитуда');

% Полюс-Ноль график (Pole/Zero Plot)
subplot(3, 2, 5);
zplane(b, a);
title('Расположение нулей и полюсов (Pole/Zero Plot)');

% Установка общего заголовка
sgtitle('Характеристики фильтра');

% Настройка размера окна
set(gcf, 'Position', [100, 100, 800, 600]);

% 4. Разложение функции передачи на простые дроби (для IIR-фильтра)
% Пример IIR-фильтра
b_iir = [0.2];      % Коэффициенты числителя
a_iir = [1, -0.5];  % Коэффициенты знаменателя

% Разложение на простые дроби
[R, P, K] = residuez(b_iir, a_iir);

% Вывод результатов
disp('Коэффициенты простых дробей:');
disp('R (резиды):');
disp(R);
disp('P (полюса):');
disp(P);
disp('K (остатки):');
disp(K);

% Аналитическая формула для импульсной характеристики
syms z;
H = poly2sym(b_iir, z) / poly2sym(a_iir, z); % Функция передачи H(z)
H_simplified = simplify(H);

disp('Функция передачи H(z):');
disp(H_simplified);

% 5. Анализ прямой формы реализации дискретного фильтра
% Получаем выходной сигнал и максимальное значение
output_signal_direct = filter(b, a, input_signal);
max_direct_output = max(abs(output_signal_direct));
disp(['Максимальное по модулю значение выходного сигнала (прямая форма): ', num2str(max_direct_output)]);

% 6. Анализ канонической формы реализации дискретного фильтра
% Для канонической формы
a_canonical = [1, -0.5]; % Примерный знаменатель
output_signal_canonical = filter(b_iir, a_iir, input_signal); % Используем IIR-фильтр

% Визуализация выходного сигнала канонической формы
figure;
subplot(2, 1, 1);
stem(t, output_signal_canonical, 'filled');
title('Выходной сигнал в канонической форме');
xlabel('Время (с)');
ylabel('Амплитуда');
grid on;

% Максимальное значение выходного сигнала канонической формы
max_canonical_output = max(abs(output_signal_canonical));
disp(['Максимальное по модулю значение выходного сигнала (каноническая форма): ', num2str(max_canonical_output)]);


% Анализ транспонированной формы реализации дискретного фильтра
output_signal_transposed = zeros(size(input_signal));
state = zeros(length(b_iir)-1, 1); % Инициализация состояний (длина должна быть length(b_iir)-1)

for n = 1:length(input_signal)
    % Вычисление выходного сигнала
    output_signal_transposed(n) = b_iir(1) * input_signal(n) + state' * b_iir(2:end)'; % Здесь размеры должны совпадать
    
    % Обновление состояний
    if n > 1
        state = [input_signal(n); state(1:end-1)] + a_iir(2:end)' * output_signal_transposed(n); % Обновление состояний
    else
        state(1) = input_signal(n); % Первое состояние
    end
end

% Визуализация выходного сигнала транспонированной формы
subplot(2, 1, 2);
stem(t, output_signal_transposed, 'filled');
title('Выходной сигнал в транспонированной форме');
xlabel('Время (с)');
ylabel('Амплитуда');
grid on;

% Максимальное значение выходного сигнала транспонированной формы
max_transposed_output = max(abs(output_signal_transposed));
disp(['Максимальное по модулю значение выходного сигнала (транспонированная форма): ', num2str(max_transposed_output)]);

% Сохранение графиков в файл (по желанию)
saveas(gcf, 'FilterCharacteristics.png');