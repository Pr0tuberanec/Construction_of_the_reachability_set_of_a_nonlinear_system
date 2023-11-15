%%
%main_1 (для единичного запуска с конкретными параметрами)
clear, clc;

[T, alpha] = ParamsInput();

hold on;
[X, Y, X1, Y1, X2, T2] = reachset(alpha, T); 

%Отображение границы множества достижимости и кривой переключения
plot(X, Y, 'Color', 'b', 'LineWidth', 2);
plot(X1, Y1, '.', 'Color', 'r', 'LineWidth', 2); 

%Отображение особых точек
plot(X2(1, T2 == 1), X2(2, T2 == 1), '.', 'Color', 'g', 'MarkerSize', 10);
plot(X2(1, T2 == 2), X2(2, T2 == 2), '.', 'Color', 'black', 'MarkerSize', 10);

%Настройка выходного окна
left_x = min(X) - 1;
right_x = max(X) + 1;
left_y = min(Y) - 1;
right_y = max(Y) + 1;

xlim([left_x, right_x]);
ylim([left_y, right_y]);
xlabel('x_1');
ylabel('x_2');

hold off; 
%%
%main_2 (для вывода анимации в файл или на экран)
clear, clc;

%Строит границы мн-ва дост-ти для моментов времени linspace(t1, t2, N)
[t1, t2, N, alpha, filename] = ParamsCin();
reachsetdyn(alpha, t1, t2, N, filename);

%%
function [T, alpha] = ParamsInput()

    T = input('Введите момент времени T, в который необходимо построить мн-во дост-ти\n');
    alpha = input('Введите alpha - параметр системы\n');
    
end

%Построение границы множества достижимости
function [X, Y, X1, Y1, X2, T2] = reachset(alpha, T)

    X = [];
    Y = [];
    X1 = [];
    Y1 = [];
    X2 = [];
    T2 = [];

    %Строим часть границы мн-ва дост-ти, которая достигается на траекториях,
    %берущих своё начало из решения S_{+} до первого обнуления скорости
    [X, Y, X1, Y1] = reachsetPartBuilding("S_{+}", alpha, T, X, Y);
    %Аналогичная часть для S_{-}
    [X, Y, X11, Y11] = reachsetPartBuilding("S_{-}", alpha, T, X, Y);
   
    %Собираем результаты в кучу
    X1 = [rot90(rot90(X11)); X1];
    Y1 = [rot90(rot90(Y11)); Y1];
    X = [X; X(1)];
    Y = [Y; Y(1)];
    
    %Удаляем самопересечения
    [X, Y] = SelfIntersectionsRemove(X, Y);

    %Ищем особые точки, только те, которые потребуется вывести
    left = min(X) - alpha - 1;
    right = max(X) + alpha + 1;
    left = min(left, -right);
    right = max(right, -left);

    s0 = [left : 0.1 : right];
    eps = 0.001;
    f = @(x) x .* cos(x.^2) - alpha;
    for i = 1 : size(s0, 2)
        options = optimoptions('fsolve','Display','off');
        a = fsolve(f, s0(i), options);
        if abs(a * cos(a^2) - alpha) < eps
            X2(end + 1) = a;
        end
    end
    X2 = round(X2, 4);
    X2 = unique(X2);
    T2 = ones(1, size(X2, 2));

    %Свойство особых точек для данной системы
    tmp = -X2;
    X2 = [X2, tmp];
    T2 = [T2, ones(1, size(tmp, 2)) + 1];

    X2 = [X2; zeros(1, size(X2, 2))];
end

function [X, Y, X1, Y1] = reachsetPartBuilding(str, alpha, T, X, Y) 

    x0 = [0, 0];
    tspan = [0, T];
    %Можно сильнее детализировать границу множества достижимости, если уменьшить шаг сетки
    %tspan = [0 : 0.0025 : T];
    if str == "S_{+}"
        dxdt = @(t, x) [x(2);
                        -x(1).*cos(x(1).^2) - 2.*x(1).*x(2) + alpha;]; 
        k = 1;
    else
        dxdt = @(t, x) [x(2);
                        -x(1).*cos(x(1).^2) - 2.*x(1).*x(2) - alpha;];
        k = 0;
    end
    [t, x] = ode45(dxdt, tspan, x0, odeset('Events', @VEvent));  
    
    t_switch = t(2 : end);
    x_0 = x(2 : end, 1);
    y_0 = x(2 : end, 2);
    X1 = [0; x_0];
    Y1 = [0; y_0];

    while any(t_switch - T < 0)

        for i = 1 : size(t_switch, 1)
            if (mod(k, 2) == 1) && (t_switch(i) < T)
                [t_switch(i), x_0(i), y_0(i)] = S_Solve("S_{-}", t_switch(i), x_0(i), y_0(i), T, alpha);
            else if (t_switch(i) < T)
                    [t_switch(i), x_0(i), y_0(i)] = S_Solve("S_{+}", t_switch(i), x_0(i), y_0(i), T, alpha);
                end
            end
        end

        X1 = [X1; x_0(t_switch - T < 0)];
        Y1 = [Y1; y_0(t_switch - T < 0)];
        i = i + 1;
    end
    
    X = [X; x_0];
    Y = [Y; y_0];        
end

%Функция события: обнуление x_2
function [value, isterminal, direction] = VEvent(t, x)
    value = x(2);  
    isterminal = 1;   
    direction = 0;   
end

function [t_switch, x_1, y_1] = S_Solve(str, t_switch, x_0, y_0, T, alpha) 
                                                  
    x0 = [x_0, y_0, -1, 0];
    tspan = [t_switch, T];                    
    if str == "S_{+}" 
        
        dxdt = @(t, x) [x(2);
                        -x(1).*cos(x(1).^2) - 2.*x(1).*x(2) + alpha;
                        2.*x(2).*x(4) + cos(x(1).^2).*x(4) - 2.*x(4).*(x(1).^2).*sin(x(1).^2);
                        -x(3) + 2.*x(1).*x(4);];         
        [t, x] = ode45(dxdt, tspan, x0, odeset('Events', @PsiEvent));
        
        %Возможность отобразить построенные по данной системе траектории
        %plot(x(:, 1), x(:, 2), 'Color', 'green');
    else
        
        dxdt = @(t, x) [x(2);
                        -x(1).*cos(x(1).^2) - 2.*x(1).*x(2) - alpha;
                        2.*x(2).*x(4) + cos(x(1).^2).*x(4) - 2.*x(4).*(x(1).^2).*sin(x(1).^2);
                        -x(3) + 2.*x(1).*x(4);];         
        [t, x] = ode45(dxdt, tspan, x0, odeset('Events', @PsiEvent));
        
        %Возможность отобразить построенные по данной системе траектории 
        %plot(x(:, 1), x(:, 2), 'Color', 'black');
    end
    
    x_1 = x(end, 1);
    y_1 = x(end, 2);
    t_switch = t(end);   
end

%Функция события: обнуление psi_2
function [value, isterminal, direction] = PsiEvent(t, x)
    value = x(4);  
    isterminal = 1;   
    direction = 0;   
end

function [X, Y] = SelfIntersectionsRemove(X, Y)

    i = 1;
    while i <= size(X, 1) - 1
        
        %Определяем уравнение прямой, содержащей i-ый отрезок
        line_i = lineBuilding(i, X, Y);
        
        j = 1;
        while j <= size(X, 1) - 1      
            
            if (i ~= j) && (i ~= j - 1) && (i ~= j + 1)

                %Если концы j-го отрезка лежат по разные стороны построенной 
                %ранее прямой, то находим прямую, содержащую j-ый отрезок 
                res_1 = line_i(1)*X(j) + line_i(2) - Y(j);
                res_2 = line_i(1)*X(j+1) + line_i(2) - Y(j+1);

                if (res_1 * res_2 < 0)
                    
                    line_j = lineBuilding(j, X, Y);      
                    
                    %Если концы i-го отрезка лежат по разные стороны прямой, содержащей j-ый
                    %отрезок, то находим точку пересечения отрезков и обновляем границу мн-ва дост-ти
                    res_3 = line_j(1)*X(i) + line_j(2) - Y(i);
                    res_4 = line_j(1)*X(i+1) + line_j(2) - Y(i+1);
                    
                    if (res_3 * res_4 < 0)
                        
                        A = [line_i(1), -1; line_j(1), -1];
                        B = [-line_i(2); -line_j(2)];
                        E = linsolve(A, B);
                        
                        X(i + 1 : j) = [];
                        X = [X(1 : i); E(1); X(i+1 : end)];
                        Y(i + 1 : j) = [];
                        Y = [Y(1 : i); E(2); Y(i+1 : end)];  
                    end
                end
            end
            
            j = j + 1;
        end
        
        i = i + 1;
    end
end

function line_params = lineBuilding(i, X, Y) 
    A = [X(i), 1; X(i + 1), 1];
    B = [Y(i); Y(i+1)];
    line_params = linsolve(A, B);   
end

function [t1, t2, N, alpha, filename] = ParamsCin()

    t1 = input('Введите начальный момент времени t1, от которого будет рассматриваться динамика мн-ва дост-ти\n');
    t2 = input('Введите конечный момент времени t2, до которого будет рассматриваться динамика мн-ва дост-ти\n');
    N = input('Введите N - число отсчётов, в которых будет происходить анимация\n');
    alpha = input('Введите alpha - параметр системы\n');
    filename = input('Введите имя файла, в который вы хотите записать анимацию\n');
    
end

%Анимация изменения множества достижимости со временем в файл или на экран
function reachsetdyn(alpha, t1, t2, N, filename)

    t = linspace(t1, t2, N);
    
    hold on;
    xlabel('x_1');
    ylabel('x_2');
    legend('Location', 'southeast');
    str = strcat('\alpha = ', ' ', num2str(alpha));        
    title(str);   
    
    if filename
        for i = 1 : N 
            [X, Y] = reachset(alpha, t(i));
            str = strcat('t = ', num2str(t(i)));
            plot(X, Y, 'LineWidth', 2, 'DisplayName', str);
            frame(i) = getframe(gcf)           
        end

        v = VideoWriter(filename);
        v.FrameRate = 1/3;
        open(v);
        writeVideo(v,frame);
        close(v);
    else
        for i = 1 : N 
            [X, Y] = reachset(alpha, t(i));
            str = strcat('t = ', num2str(t(i)));
            plot(X, Y, 'LineWidth', 2, 'DisplayName', str);
            pause(3);
        end
    end
    hold off; 
end
