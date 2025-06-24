%% 基础作图
% 参数设置
r0 = 0;  % 初始半径偏移
p = 0.55;  % 螺旋增长的速度
theta_max = 32 * pi;  % 螺旋的终止角度

% 生成螺旋的角度和半径数据
theta = linspace(0, theta_max, 1000);  
r = r0 +  (p * theta) / (2 * pi);  

% 转换为笛卡尔坐标系
x = r .* cos(theta);
y = r .* sin(theta);

% 绘图
figure;
plot(x, y);
axis equal;  % 保持横纵坐标刻度一致
xlabel('x');
ylabel('y');
title('等距螺线');
grid on;  % 显示网格

% 去掉横纵坐标上的数字
xticks([]);
yticks([]);

% 坐标轴
hold on;
line(xlim, [0 0], 'Color', 'k');  % X轴
line([0 0], ylim, 'Color', 'k');  % Y轴
hold off;

%% 获得所有位置
% 设置时间的步长

i = 1;

T = [0:i:500];
location = get_location(T);

%% 获得所有速度
% 设置时间的步长

i = 1;

T = [0:i:500];
speed = get_speed(T);

% 问题一论文写作提交格式的位置表格

part_location = zeros(7*2, length(T));

for i = 1:length(T)
    part_location(1,i) = location(1,i);
    part_location(2,i) = location(2,i);
    part_location(3,i) = location(3,i);
    part_location(4,i) = location(4,i);
    part_location(5,i) = location(103,i);
    part_location(6,i) = location(104,i);
    part_location(7,i) = location(203,i);
    part_location(8,i) = location(204,i);
    part_location(9,i) = location(303,i);
    part_location(10,i) = location(304,i);
    part_location(11,i) = location(403,i);
    part_location(12,i) = location(404,i);
    part_location(13,i) = location(447,i);
    part_location(14,i) = location(448,i);
end

%% 问题一论文写作提交格式的速度表格

part_speed = zeros(7, length(T));

for i = 1:length(T)
    part_speed(1,i) = speed(1,i+1);
    part_speed(2,i) = speed(2,i+1);
    part_speed(3,i) = speed(52,i+1);
    part_speed(4,i) = speed(102,i+1);
    part_speed(5,i) = speed(152,i+1);
    part_speed(6,i) = speed(202,i+1);
    part_speed(7,i) = speed(224,i+1);

end

%% 问题二论文写作提交格式的位置和速度表格

T = [414];
location2 = get_location(T);

speed2 = get_speed(T);

result2 = zeros(224, 3);

result2(:,3) = speed2(:,2);

for i = 1:224
    result2(i,1) = location2(2*i-1, 1);
    result2(i,2) = location2(2*i, 1);

end
%% 获得所有速度的函数
% r(theta) * (\delta(theta)) / (\delta(t)) = (\delta(l)) / (\delta(t))
function result = get_speed(T)

    flag = 224;

    location_theta = get_location_theta(T);
    
    delta_t = 0.001;
    delta_T = T + delta_t;
    
    delta_location_theta = get_location_theta(delta_T);
    
    d_location_theta = location_theta - delta_location_theta;
    
    dtheta_dt = d_location_theta ./ delta_t;
    r_theta = get_r( (0.5 .* T + 0.5 .* delta_T) );
    
    speed = r_theta .* dtheta_dt;
    
    for i = 1:flag
        speed(i, 1) = i - 1;
    end

    result = speed;
end
%% 获得所有极径的函数

function result = get_r(T)

    a = 16 * 0.55;
    b = 0.55 / (2 * pi);

    location_theta = get_location_theta(T);
    result = a + b .* location_theta;

end

%% 获得所有位置(单参数)的函数
function result = get_location_theta(T)

    % 定义基本参数及函数
    a = 16 * 0.55;
    b = 0.55 / (2 * pi);
    
    % 龙头板登前把手和后把手的固定距离
    c0 = 3.41 - 2 * 0.275;
    % 其余板登前把手和后把手的固定距离
    c1 = 2.20 - 2 * 0.275;
    
    % 龙头的线速度
    v = 1;
    
    % 定义迭代次数
    flag = 224;  
    
    theta = zeros(flag, length(T)+1);
    
    for i = 1:flag
        theta(i, 1) = i - 1;
    end
    
    for i = 2:(length(T)+1)
        theta(1, i) = get_theta0(v, T(i-1));
        theta(2, i) = location_rec(c0, theta(1, i));
    
        for j = 2:(flag-1)
        theta(j+1, i) = location_rec(c1, theta(j, i));
        end
    
    end

    result = theta;
end

%% 获得所有位置(直角坐标系)的函数

function result = get_location(T)

    % 定义基本参数及函数
    a = 16 * 0.55;
    b = 0.55 / (2 * pi);
    
    % 龙头板登前把手和后把手的固定距离
    c0 = 3.41 - 2 * 0.275;
    % 其余板登前把手和后把手的固定距离
    c1 = 2.20 - 2 * 0.275;
    
    % 龙头的线速度
    v = 1;
    
    % 定义迭代次数
    flag = 224;  
    
    theta = zeros(flag, length(T)+1);
    
    for i = 1:flag
        theta(i, 1) = i - 1;
    end
    
    for i = 2:(length(T)+1)
        theta(1, i) = get_theta0(v, T(i-1));
        theta(2, i) = location_rec(c0, theta(1, i));
    
        for j = 2:(flag-1)
        theta(j+1, i) = location_rec(c1, theta(j, i));
        end
    
    end
    
    all_location = zeros(flag*2, length(T));
    
    for i = 1:length(T) 
        for j = 1:flag
        [all_location(2*j-1, i), all_location(2*j, i)]= change(theta(j, i+1));
        end
    end
    
    result = all_location;
end

%% 计算经过时间t龙头的位置函数

function theta0 = get_theta0(v, t)
    
    % 等距螺线基本参数
    a = 16 * 0.55;
    b = 0.55 / (2 * pi);

    % 定义 r(theta)
    r = @(theta) a + b * theta;

    % 它的导数 dr/dtheta
    dr_dtheta = @(theta) b;
    
    %龙头在t时间内走过的弧长
    L = v * t;

    % 定义被积函数
    integrand = @(theta) sqrt((dr_dtheta(theta)).^2 + (r(theta)).^2);
    
    % 定义积分方程：目标是求解 theta0，使得积分结果等于 v*t
    integral_eq = @(theta0) integral(integrand, theta0, 0) - L;
    
    % 使用 fsolve 求解 theta0
    theta0_initial_guess = 0; 
    result = fzero(integral_eq, theta0_initial_guess);
    
    theta0 = result;

end

%% 极坐标系转换成直角坐标系函数

function [x, y] = change(theta)

    % 等距螺线基本参数
    a = 16 * 0.55;
    b = 0.55 / (2 * pi);

    % 定义 r(theta)
    r = @(theta) a + b * theta;

    x = r(theta) * cos(theta);
    y = r(theta) * sin(theta);

end

%% 通过距离约束得到位置递推函数

function [theta_new] = location_rec(c, theta)

    % 等距螺线基本参数
    a = 16 * 0.55;
    b = 0.55 / (2 * pi);

    % 定义 r(theta)
    r = @(theta) a + b * theta;    

   % 计算已知 theta 对应的坐标
    [x_known, y_known] = change(theta);

    % 方程定义
    eq_fun = @(vars) [
        (vars(2) * cos(vars(1)) - x_known)^2 + ...
        (vars(2) * sin(vars(1)) - y_known)^2 - c^2; 
        vars(2) - (a + b * vars(1))
    ];

    % 初始猜测
    guess = [theta+0.1, r(theta)];

    % 使用 fsolve 求解方程
    options = optimoptions('fsolve', 'Display', 'off'); % 可选项
    result = fsolve(eq_fun, guess, options);
    
    % 输出结果
    theta_new = result(1);

end

