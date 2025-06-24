clc, clear

% 基础参数设置
A = 0;
B = 1.7 / (2 * pi);

r_in = @(theta) A + B * theta;
r_out = @(theta) A - B * theta;

L = 13.621245;
v = 1;


loc = zeros(224 * 2, 200);

% 龙头随时间迭代
T1 = linspace(-100, -1, 100);
theta1 = get_theta1(v, T1);
for i = 1:length(theta1)
    loc(1, i) = r_in(theta1(i)) * cos(theta1(i));
    loc(2, i) = r_in(theta1(i)) * sin(theta1(i));
end

T = floor(L / v);
T2 = linspace(0, T, T + 1);
for i = 1:length(T2)
    [result_x, result_y] = get_location_stage_s(v, T2(i));
    loc(1, 100 + i) = result_x;
    loc(2, 100 + i) = result_y;
end

T3 = linspace(T + 1, 100, 100 - T);
for i = 1:length(T3)
    theta2 = get_theta2(v, T3(i));
    loc(1, 101 + T + i) = r_out(theta2) * cos(theta2);
    loc(2, 101 + T + i) = r_out(theta2) * sin(theta2);
end

% 龙身随空间迭代
for t = 1:201
    loc = iterative_process(loc(1, t), loc(2, t), t, loc);
end

%% 求速度
delta_loc = zeros(224 * 2, 200);
speed = zeros(224*2, 200);

delta_t = 0.001;
delta_T1 = linspace(-100 + delta_t , -1 + delta_t, 100);
theta1 = get_theta1(v, delta_T1);
for i = 1:length(theta1)
    delta_loc(1, i) = r_in(theta1(i)) * cos(theta1(i));
    delta_loc(2, i) = r_in(theta1(i)) * sin(theta1(i));
end

T = 13;
delta_T2 = linspace(0 + delta_t, T + delta_t, T + 1);
for i = 1:length(delta_T2)
    [result_x, result_y] = get_location_stage_s(v, delta_T2(i));
    delta_loc(1, 100 + i) = result_x;
    delta_loc(2, 100 + i) = result_y;
end

delta_T3 = linspace(T + 1 + delta_t, 100 + delta_t, 100 - T);
for i = 1:length(delta_T3)
    theta2 = get_theta2(v, delta_T3(i));
    delta_loc(1, 101 + T + i) = r_out(theta2) * cos(theta2);
    delta_loc(2, 101 + T + i) = r_out(theta2) * sin(theta2);
end

for t = 1:201
    delta_loc = iterative_process(delta_loc(1, t), delta_loc(2, t), t, delta_loc);
end

[rows, cols] = size(loc);
for i = 1:rows
    for j = 1:cols
        speed(i, j) = (delta_loc(i, j) - loc(i, j)) / delta_t;
    end
end

[m, n] = size(speed);

speed_magnitudes = zeros(m/2, n);

% 计算速度大小
for col = 1:n
    vx = speed(1:2:end, col);  % 获取 x 方向速度 (奇数行)
    vy = speed(2:2:end, col);  % 获取 y 方向速度 (偶数行)
    
    % 速度大小 = sqrt(vx^2 + vy^2)
    speed_magnitudes(:, col) = sqrt(vx.^2 + vy.^2);
end

%% 不同时刻的龙头迭代
%% 龙头在 -100-0 时间内在第一段等距螺线位置函数
function x = get_theta1(v, t)
% 等螺线基本参数
b = 1.7 / (2 * pi);
theta0 = 16.631961107240077;
r = @(theta) b * theta;
dr_dtheta = @(theta) b;
%龙头在t时间内走过的弧长
for i = 1:length(t)
    L = -v .* t(i);
    % 定义被积函数
    integrand = @(theta) sqrt((dr_dtheta(theta)).^2 + (r(theta)).^2);
    integral_eq = @(x) integral(integrand, theta0, x) - L;
    % 使用 fsolve 求解 theta0
    theta0_initial_guess = 16*pi; 
    result = fzero(integral_eq, theta0_initial_guess);
    x(1, i) = result;
end
end
%% 龙头在 0-T 时间内在S形圆弧的直角坐标位置函数
function [result_x, result_y] = get_location_stage_s(v, t)
x0 = -2.711855863706671;
y0 = -3.591077522761064;
x1 = -0.760009116655550;
y1 = -1.305726426346273;
x2 = 1.735932490181134;
y2 = 2.448401974553573;
% l <= l_long + l_short
l = v * t;
R = 1.502708833894531;
theta = 3.021486841980670;
l_long = 2 * R * theta;
l_short = R * theta;
if l <= l_long
delta_theta = l / l_long * theta;
% 顺时针转动delta_theta
Location = [x1; y1] + [cos(delta_theta), -sin(-delta_theta); sin(-delta_theta), cos(delta_theta)] * [x0 - x1; y0 - y1];
result_x = Location(1);
result_y = Location(2);
else 
l_extra = l - l_long;
delta_theta = theta - l_extra / l_short * theta;
% 顺时针转动delta_theta
Location = [x2; y2] + [cos(delta_theta), -sin(-delta_theta); sin(-delta_theta), cos(delta_theta)] * [-x0 - x2; -y0 - y2];
result_x = Location(1);
result_y = Location(2);
end
end
%% 龙头在 T-100 时间内在第二段等距螺线位置函数
function theta0 = get_theta2(v, t)
% 等螺线基本参数
b = 1.7 / (2 * pi);
thetak = 16.631961107240077+pi; % 给定的终点角度
T = 13.621244906821248;
% 定义 r(theta) 螺旋线的半径
r = @(theta) b * (theta-pi);
% 定义 r(theta) 的导数 dr/dtheta
dr_dtheta = @(theta) b;
% 龙头在 t 时间内走过的弧长 L (负号表示从外向内运动)
L = v * (t-T);
% 定义被积函数，用于计算弧长
integrand = @(theta) sqrt((dr_dtheta(theta)).^2 + (r(theta)).^2);
% 定义积分方程，目标是求解 theta0
integral_eq = @(x) integral(integrand, thetak, x) - L;
% 使用 fzero 求解 theta0
theta0_initial_guess = 0; % 初始猜测点
result = fzero(integral_eq, theta0_initial_guess);
% 计算出的 theta0 
theta0 = result - pi;
end

%% 不同位置的龙身迭代
%%  逆螺线迭代函数1
function [d, e, exitflag] = spiral_spiral1(c, a, b, a_last, b_last)
    exitflag = 1;
    % 等距螺线基本参数
    A = 0;
    B = 1.7 / (2 * pi);

    eqn = @(var)[
        (var(1) * cos(var(2)) - a)^2 + (var(1) * sin(var(2)) - b)^2 - c^2;
        var(1) - A - (-B * var(2));
    ];
    
    % 初始猜测值
    guess = [-sqrt(a^2 + b^2), -sqrt(a^2 + b^2) / (-B) + 1];

    % 创建 fsolve 的选项
    options = optimoptions('fsolve', 'Display', 'off');  % 关闭显示

    % 调用 fsolve 进行求解
    result = fsolve(eqn, guess, options);
    
    % 结果
    d_temp = result(1) * cos(result(2));
    e_temp = result(1) * sin(result(2));

    % 检查是否需要调整解
    if [d_temp - a, e_temp - b] * [a_last - a, b_last - b].' < 0
        % "option A"
        d = d_temp;
        e = e_temp;
    else
        % "option B"
        % 更新猜测值，并再次求解
        % guess = [a + c, b + c];
        guess = [-sqrt(a^2 + b^2), -sqrt(a^2 + b^2) / (-B) - 1];
        result = fsolve(eqn, guess, options);
        d = result(1) * cos(result(2));
        e = result(1) * sin(result(2));
    end
    hold on;
    plot(d, e);

    % 如果解的模小于等于 4.5, 设置 exitflag 为 -1
    if sqrt(d^2 + e^2) <= 4.5
        sqrt(d^2 + e^2);
        exitflag = -1;
    end
end

function x = check_spiral(exitflag)
    if exitflag <= 0
        fprintf(" option 1 ")
        x = 0;
    else
        x = 1;
    end
end
%% 内部圆弧迭代函数
function [d,e,exitflag] = smallcircle_smallcircle(c,a,b,a_last,b_last)
    exitflag = 1;
    % 已知坐标
    x0 = -2.711855863706671;
    y0 = -3.591077522761064;
    x1 = -0.760009116655550;
    y1 = -1.305726426346273;
    x2 = 1.735932490181134;
    y2 = 2.448401974553573 ;
    r = 1.502708833894531;

    % 定义方程组，注意这里的 vars 是一个向量 [x, y]
    eq_fun = @(vars) [
        (a - vars(1))^2 + (b - vars(2))^2 - c^2; 
        (vars(1) - x2)^2 + (vars(2) - y2)^2 - r^2;
    ];

    % 初始猜测
    guess = [a - c, b - c];
    
    % 优化选项
    options = optimoptions('lsqnonlin', 'Display', 'off');  % 使用 lsqnonlin 的选项
    % 如果没有上下限，可以设为 []
    lb = [0.90395, -inf];
    ub = [inf, inf];
    % 使用 lsqnonlin 求解
    [result, ~, ~, exitflag] = lsqnonlin(eq_fun, guess, lb, ub, options);

    % 结果
    d_temp = result(1);
    e_temp = result(2);

    % size([d_temp - a, e_temp - b])
    % size([a_last - a, b_last - b])
    if [d_temp - a, e_temp - b] * [a_last - a, b_last - b].' < 0
        % "option A"
        % [d_temp - a, e_temp - b] * [a_last - a, b_last - b].'
        d = d_temp;
        e = e_temp;
    else
        % "option B"
        % [d_temp - a, e_temp - b] * [a_last - a, b_last - b].'
        guess = [a + c, b + c];
        % result = fsolve(eq_fun, guess, options);
        [result, ~, ~, exitflag] = lsqnonlin(eq_fun, guess, lb, ub, options);
        d = result(1);
        e = result(2);
        % [d - a, e - b] * [a_last - a, b_last - b].'

        if d == d_temp && e == e_temp
            exitflag = -1;
        end
    end
end

function [d,e,exitflag] = bigcircle_bigcircle(c,a,b,a_last,b_last)
    % 已知坐标
    exitflag = 1;
    x0 = -2.711855863706671;
    y0 = -3.591077522761064;
    x1 = -0.760009116655550;
    y1 = -1.305726426346273;
    x2 = 1.735932490181134;
    y2 = 2.448401974553573 ;
    r = 1.502708833894531 * 2;  % 大圆的半径

    % 定义方程组，注意这里的 vars 是一个向量 [x, y]
    eq_fun = @(vars) [
        (a - vars(1))^2 + (b - vars(2))^2 - c^2; 
        (vars(1) - x1)^2 + (vars(2) - y1)^2 - r^2;
    ];

    % 初始猜测
    guess = [a+c, b+c];

    % 优化选项
    options = optimoptions('lsqnonlin', 'Display', 'off');  % 使用 lsqnonlin 的选项
    % 如果没有上下限，可以设为 []
    lb = [-4, -inf];
    ub = [0.90395, inf];
    % 使用 lsqnonlin 求解
    [result, ~, ~, exitflag] = lsqnonlin(eq_fun, guess, lb, ub, options);

    % 结果
    d_temp = result(1);
    e_temp = result(2);

    % size([d_temp - a, e_temp - b])
    % size([a_last - a, b_last - b])
    if [d_temp - a, e_temp - b] * [a_last - a, b_last - b].' < 0
        d = d_temp;
        e = e_temp;
    else
        guess = [a-c, b-c];
        % result = fsolve(eq_fun, guess, options);
        [result, ~, ~, exitflag] = lsqnonlin(eq_fun, guess, lb, ub, options);
        d = result(1);
        e = result(2);
        if d == d_temp && e == e_temp
            exitflag = -1;
        end
    end
end

function [x, check_degree] = check_small(check_degree,c,d,e,exitflag)
x0 = -2.711855863706671;
y0 = -3.591077522761064;
x1 = -0.760009116655550;
y1 = -1.305726426346273;
x2 = 1.735932490181134;
y2 = 2.448401974553573 ;
r = 1.502708833894531;
theta1 = 3.021486841980670;

delta_degree = 2 * asin(c / (2 * r));
if check_degree == -114514
    check_degree = acos( ((d-x2)*(-x0-x2)+(e-y2)*(-y0-y2))/(r*r));
else
    check_degree = check_degree + delta_degree;
end

if exitflag<=0
    fprintf(" option 1 ");
    x=0;
    return;
end 

if   check_degree > theta1
    fprintf(" option 2 ");
    x=0;
    return;
end    
x=1;    
end

function [x, check_degree] = check_big(check_degree,c,d,e,exitflag)
x0 = -2.711855863706671;
y0 = -3.591077522761064;
x1 = -0.760009116655550;
y1 = -1.305726426346273;
x2 = 1.735932490181134;
y2 = 2.448401974553573 ;
r = 1.502708833894531 * 2;
theta1 = 3.021486841980670;

delta_degree = 2 * asin(c / (2 * r));
if check_degree == -114514
    check_degree = theta1 - acos( ((d-x1)*(x0-x1)+(e-y1)*(y0-y1))/(r*r));
else
    check_degree = check_degree + delta_degree;
end

if exitflag<=0
    fprintf(" option 1\n ");
    x=0;
    return;
end 

if  check_degree > theta1
    fprintf(" option 2 ");
    x=0;
    return;
end    
x=1;    
end

%% 螺线迭代函数2
function [d, e, exitflag] = spiral_spiral2(c, a, b, a_last, b_last)
    exitflag = 1;
    % 等距螺线基本参数
    A = 0;
    B = 1.7 / (2 * pi);

    eqn = @(var)[
        (var(1) * cos(var(2)) - a)^2 + (var(1) * sin(var(2)) - b)^2 - c^2;
        var(1) - A - (B * var(2));
    ];
    
    % 初始猜测值
    guess = [sqrt(a^2 + b^2), sqrt(a^2 + b^2) / (B) - 1];

    % 创建 fsolve 的选项
    options = optimoptions('fsolve', 'Display', 'off');  % 关闭显示

    % 调用 fsolve 进行求解
    result = fsolve(eqn, guess, options);
    % 结果
    d_temp = result(1) * cos(result(2));
    e_temp = result(1) * sin(result(2));

    % 检查是否需要调整解
    if [d_temp - a, e_temp - b] * [a_last - a, b_last - b].' < 0
        % "option A"
        d = d_temp;
        e = e_temp;
    else
        % "option B"
        % 更新猜测值，并再次求解

        guess = [sqrt(a^2 + b^2), sqrt(a^2 + b^2) / (B) + 1];
        result = fsolve(eqn, guess, options);
        d = result(1) * cos(result(2));
        e = result(1) * sin(result(2));
    end
    hold on;
    plot(d, e);

    % 如果解的模小于等于 4.5, 设置 exitflag 为 -1
    if sqrt(d^2 + e^2) <= 4.5
        sqrt(d^2 + e^2);
        exitflag = -1;
    end
end

%% 整体迭代过程
function loc = iterative_process(a, b, t, loc_in)
    max_iterations = 223;  % 设置最大迭代次数防止无限循环
    % 设定初始默认变量值
    a_last = a;
    b_last = b;
    check_degree = -114514;
    
    R0 = sqrt(a^2 + b^2);

    if t <= 103
        use_spiral1 = false;
        use_smallcircle = false;
        use_bigcircle = false;
        use_spiral2 = true;
    elseif t <= 113
        use_spiral1 = false;
        use_smallcircle = false;
        use_bigcircle = true;
        use_spiral2 = true;
    elseif t <= 118
        use_spiral1 = false;
        use_smallcircle = true;
        use_bigcircle = true;
        use_spiral2 = true;
    else
        use_spiral1 = true;
        use_smallcircle = true;
        use_bigcircle = true;
        use_spiral2 = true;
    end

    for iter = 1:max_iterations
        c = 1.65;
        if iter == 1
            c = 2.86;
        end
        
        if use_spiral1

            [d, e, exitflag] = spiral_spiral1(c, a, b, a_last, b_last);

            % 检查 smallcircle_smallcircle 结果
            x = check_spiral(exitflag);

            if x == 0
                [d, e, exitflag] = smallcircle_smallcircle(c, a, b, a_last, b_last);
                use_spiral1 = false;
            end
            
            % 打印当前迭代结果
            fprintf('Iteration %d (Spiral1): d = %.6f, e = %.6f, x = %d\n', iter, d, e, x);
            plot(d, e, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

            a_last = a;
            b_last = b;
            a = d;
            b = e;
        
        elseif use_smallcircle
            % 使用 smallcircle_smallcircle 函数
            [d, e, exitflag] = smallcircle_smallcircle(c, a, b, a_last, b_last); 
            % 检查 smallcircle_smallcircle 结果
            [x, check_degree] = check_small(check_degree, c, d, e, exitflag);

            if x == 0
                [d, e, exitflag] = bigcircle_bigcircle(c, a, b, a_last, b_last);
                check_degree = -114514;
                use_smallcircle = false;
            end
            
            % 打印当前迭代结果
            fprintf('Iteration %d (Smallcircle): d = %.6f, e = %.6f, x = %d\n', iter, d, e, x);
            plot(d, e, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

            
            a_last = a;
            b_last = b;
            a = d;
            b = e;
        elseif use_bigcircle
            % 使用 bigcircle_bigcircle 函数，不再进行 check
            [d, e, exitflag] = bigcircle_bigcircle(c, a, b, a_last, b_last);
            % 检查 bigcircle_bigcircle 结果
            [x, check_degree] = check_big(check_degree, c, d, e, exitflag);
            
            if x == 0
                [d, e, exitflag] = spiral_spiral2(c, a, b, a_last, b_last);
                check_degree = -114514;
                use_bigcircle = false;
            end

            % 打印 bigcircle 阶段的迭代结果
            fprintf('Iteration %d (Bigcircle): d = %.6f, e = %.6f, x = %d\n', iter, d, e, x);
            plot(d, e, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
            

            % 更新 a, b，继续 bigcircle_bigcircle 迭代
            a_last = a;
            b_last = b;
            a = d;
            b = e;
        elseif use_spiral2
            [d, e, exitflag] = spiral_spiral2(c, a, b, a_last, b_last);
            
            x = check_spiral(exitflag);

            if x == 0
                [d, e, exitflag] = smallcircle_smallcircle(c, a, b, a_last, b_last);
                % use_spiral2 = false;
            end
            
            % 打印当前迭代结果
            fprintf('Iteration %d (Spiral2): d = %.6f, e = %.6f, x = %d\n', iter, d, e, x);
            plot(d, e, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

            a_last = a;
            b_last = b;
            a = d;
            b = e;
        
        end
        loc_in(iter * 2 + 1, t) = d;
        loc_in(iter * 2 + 2, t) = e;
        loc = loc_in;
        % 可根据需求添加终止条件
    end
end

%% 坐标系转化函数(等距螺线上)

function [x, y] = theta2loc1(theta, a, b)

    % % 等距螺线基本参数
    % a = 16 * 0.55;
    % b = 0.55 / (2 * pi);

    % 定义 r(theta)
    r = @(theta) a + b * theta;

    x = r(theta) * cos(theta);
    y = r(theta) * sin(theta);

end

function theta = loc2theta1(loc, a, b)
    r = sqrt(loc(1)^2 + loc(2)^2);
    theta = (r - a) / (-b);
end

function theta = loc2theta2(loc, a, b)
    r = sqrt(loc(1)^2 + loc(2)^2);
    theta = (r - a) / (b);
end

function [r, theta] = loc2polar(x, y)
    % 计算 r
    r = sqrt(x^2 + y^2);
    
    % 计算 theta (弧度)
    theta = atan2(y, x);
end

%% 绘图函数
function draw(X0, Y0)
% 螺距设置
a = 1.7/(2*pi);  % 螺距
theta = linspace(0, 8 * pi, 1000);  % 角度范围

% 两条中心对称的等距螺旋线
r_in = a * theta;  % 顺时针盘入螺线
r_out = -a * theta;  % 逆时针盘出螺线（中心对称）

% 将极坐标转换为直角坐标
x_in = r_in .* cos(theta);
y_in = r_in .* sin(theta);

x_out = r_out .* cos(theta);
y_out = r_out .* sin(theta);

% 找到半径为4.5m的点
radius_target = 4.5;
[~, idx_in] = min(abs(r_in - radius_target));
[~, idx_out] = min(abs(r_out + radius_target));

% 对应半径为4.5m的点的坐标
x_in_4_5 = x_in(idx_in);
y_in_4_5 = y_in(idx_in);
x_out_4_5 = x_out(idx_out);
y_out_4_5 = y_out(idx_out);

% 圆弧参数
a1=3.2611;
a2= 3.0221;
 r= 1.5014;
 x1= -0.7639;
 x2= 1.7420;
 y1= -1.3018;
 y2= 2.4409;

% 创建大圆和小圆的参数
theta_circle = linspace(0, 2*pi, 1000);  % 完整圆的角度范围

x_circle1 = x1 +2 *r * cos(theta_circle);  % 大圆的x坐标
y_circle1 = y1 + 2*r * sin(theta_circle);  % 大圆的y坐标

x_circle2 = x2 + r * cos(theta_circle);  % 小圆的x坐标
y_circle2 = y2 + r * sin(theta_circle);  % 小圆的y坐标

% 绘制螺旋线
figure;
hold on;
plot(x_in, y_in, 'b', 'DisplayName', 'Clockwise Spiral (盘入螺线)');
plot(x_out, y_out, 'g', 'DisplayName', 'Counterclockwise Spiral (盘出螺线)');

% 标记半径为4.5m的点
scatter(x_in_4_5, y_in_4_5, 100, 'filled', 'r', 'DisplayName', 'Radius 4.5m Point (In Spiral)');
scatter(x_out_4_5, y_out_4_5, 100, 'filled', 'r', 'DisplayName', 'Radius 4.5m Point (Out Spiral)');

% 标注半径为4.5m的点
text(x_in_4_5, y_in_4_5, sprintf('(%.2f, %.2f)', x_in_4_5, y_in_4_5), 'Color', 'r', 'FontSize', 12);
text(x_out_4_5, y_out_4_5, sprintf('(%.2f, %.2f)', x_out_4_5, y_out_4_5), 'Color', 'r', 'FontSize', 12);

% 绘制完整的大圆和小圆
plot(x_circle1, y_circle1, 'r', 'LineWidth', 2, 'DisplayName', 'Large Circle');
plot(x_circle2, y_circle2, 'm', 'LineWidth', 2, 'DisplayName', 'Small Circle');

% 添加图例和标签
title('Center Symmetric Spirals and Circles with Radius 4.5m Points');
xlabel('x (m)');
ylabel('y (m)');
% legend('show');
axis equal;
grid on;


hold on;
line(xlim, [0 0], 'Color', 'k'); 
line([0 0], ylim, 'Color', 'k'); 
% 添加一个圆心在原点，半径为4.5米的圆
r1_circle = 4.5; 
theta_circle = linspace(0, 2*pi, 1000);
x1_circle = r1_circle * cos(theta_circle);
y1_circle = r1_circle * sin(theta_circle);
plot(x1_circle, y1_circle, 'k--', 'LineWidth', 1.5); % 绘制圆


% 坐标轴
hold on;
line(xlim, [0 0], 'Color', 'k');  % X轴
line([0 0], ylim, 'Color', 'k');  % Y轴

plot(X0, Y0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

iterative_process(1.65,X0,Y0);

hold off;
grid on;  % 显示网格
end