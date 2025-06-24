%% 基本图形

r0 = 0;         
p = 1.7;      
theta_max = 32 * pi; 
num_points = 1000;

theta = linspace(0, theta_max, num_points);

r1 = r0 + (p * theta) / (2 * pi);
r2 = r0 - (p * theta) / (2 * pi);

x1 = r1 .* cos(theta);
y1 = r1 .* sin(theta);

x2 = r2 .* cos(theta);
y2 = r2 .* sin(theta);

figure;
plot(x1, y1, 'b', 'LineWidth', 1);  
hold on;
plot(x2, y2, 'r', 'LineWidth', 1);  
axis equal;

xticks([]);
yticks([]);

xlabel('');
ylabel('');

legend off;

hold on;
line(xlim, [0 0], 'Color', 'k');  
line([0 0], ylim, 'Color', 'k'); 

% 添加一个圆心在原点，半径为4.5米的圆
r_circle = 4.5;  
theta_circle = linspace(0, 2*pi, num_points);
x_circle = r_circle * cos(theta_circle);
y_circle = r_circle * sin(theta_circle);
plot(x_circle, y_circle, 'k--', 'LineWidth', 1.5);  % 绘制圆

hold off;

%% 求解进入掉头区域坐标

% 基本参数定义
r0 = 0;         
p = 1.7;
r1_circle = 4.5;

r1 = @(theta) r0 + (p * theta) / (2 * pi);
r2 = @(theta) -r0 - (p * theta) / (2 * pi);

f1 = @(theta) r1(theta) - r1_circle;
f2 = @(theta) r2(theta) - r1_circle;

guess = 10;  

theta0 = fzero(f1, guess);
theta3 = fzero(f2, -guess);

x0 = r1(theta0) * cos(theta0);
y0 = r1(theta0) * sin(theta0);

x3 = r2(-theta3) * cos(-theta3);
y3 = r2(-theta3) * sin(-theta3);

fprintf('螺线与圆相切点的坐标: (%f, %f), %f\n', x0, y0, theta0);
fprintf('螺线与圆相切点的坐标: (%f, %f), %f\n', x3, y3, theta3);

% 在图上标注相切点
hold on;
plot(x0, y0, 'ko', 'MarkerSize', 10, 'LineWidth', 1);
plot(x3, y3, 'ko', 'MarkerSize', 10, 'LineWidth', 1);
hold off;


%% 主程序
% 定义初始猜测值
xmin = [-1, 1, -1, 1, 1.5, 3*pi/4, 3*pi/4]; % [x1, x2, y1, y2, R, theta1, theta2] 初始猜测

% 调用 fmincon 优化
options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'sqp');

% 最小化 2 * R * theta1 + R * theta2
[x_opt, fval] = fmincon(@objective, xmin, [], [], [], [], [], [], @constraints, options);

% 输出最优结果
disp('最优解:');
disp(x_opt);
disp('最小化目标函数的值:');
disp(fval);

%% 进一步作图

hold on;

x1 = x_opt(1);
x2 = x_opt(2);
y1 = x_opt(3);
y2 = x_opt(4);
R = x_opt(5);
theta1 = x_opt(6);
theta2 = x_opt(7);

num_points = 1000;
theta_circle = linspace(0, 2*pi, num_points);

r1_circle = 2 * R;  
x1_circle = x1 + r1_circle * cos(theta_circle);
y1_circle = y1 + r1_circle * sin(theta_circle);
plot(x1_circle, y1_circle, 'k--', 'LineWidth', 1.5); 

r2_circle = R;  
x2_circle = x2 + r2_circle * cos(theta_circle);
y2_circle = y2 + r2_circle * sin(theta_circle);
plot(x2_circle, y2_circle, 'k--', 'LineWidth', 1.5); 

plot([x0, x1], [y0, y1], 'k--', 'LineWidth', 1.5);
plot([x2, x1], [y2, y1], 'k--', 'LineWidth', 1.5);
plot([x3, x2], [y3, y2], 'k--', 'LineWidth', 1.5);

hold off;

%% 目标函数
function obj = objective(x)
    R = x(5);
    theta1 = x(6);
    theta2 = x(7);
    
    % 目标函数: 2 * R * theta1 + R * theta2
    obj = 2 * R * theta1 + R * theta2;
end

%% 约束方程
function [c, ceq] = constraints(x)
    x1 = x(1);
    x2 = x(2);
    y1 = x(3);
    y2 = x(4);
    R = x(5);
    theta1 = x(6);
    theta2 = x(7);
    x0 = -2.711855863706671;
    y0 = -3.591077522761064;
    theta0 = 16.631961107240077;
    p = 1.7;

    % k1(theta0) 函数定义
    k1 = @(theta) ( (p / (2 * pi)) * sin(theta) + (p * theta) / (2 * pi) * cos(theta)) / ...
                  ( (p / (2 * pi)) * cos(theta) - (p * theta) / (2 * pi) * sin(theta));

    % 约束1: (x1 - x0)^2 + (y1 - y0)^2 = (2 * R)^2
    ceq(1) = (x1 - x0)^2 + (y1 - y0)^2 - (2 * R)^2;
    
    % 约束2: (x2 + x0)^2 + (y2 + y0)^2 = R^2
    ceq(2) = (x2 + x0)^2 + (y2 + y0)^2 - R^2;
    
    % 约束3: [1, k1(theta0)] .* [x0 - x1, y0 - y1] = 0
    ceq(3) = [1, k1(theta0)] * [x0 - x1; y0 - y1];
   
    % 约束4: cos(theta1) = ([x0-x1, y0-y1] .* [x2-x1, y2-y1]) / (norm([x0-x1, y0-y1]) * norm([x2-x1, y2-y1]))
    ceq(4) = cos(theta1) - dot([x0 - x1, y0 - y1], [x2 - x1, y2 - y1]) / ...
           (norm([x0 - x1, y0 - y1]) * norm([x2 - x1, y2 - y1]));
    
    % 约束5: cos(theta2) = ([-x0-x2, -y0-y2] .* [x1-x2, y1-y2]) / (norm([-x0-x2, -y0-y2]) * norm([x1-x2, y1-y2]))
    ceq(5) = cos(theta2) - dot([-x0 - x2, -y0 - y2], [x1 - x2, y1 - y2]) / ...
           (norm([-x0 - x2, -y0 - y2]) * norm([x1 - x2, y1 - y2]));
    
    ceq(6) = (x1 - x2)^2 + (y1 - y2)^2 - (3 * R)^2;
    
    ceq(7) = [1, k1(theta0)] * [x0 + x2; y0 + y2];

    % 不等式约束为空
    c = [];
end
