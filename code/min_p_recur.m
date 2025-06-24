% 初始螺距参数
min_pitch = 0.43;  % 最小螺距
max_pitch = 0.46;  % 已知数据螺距为 0.55 cm
tolerance = 0.002; % 优化的精度
center_radius = 4.5; % 中心区域半径
half_width = 0.15;  % 设定平移的距离
extension = 0.275;  % 设定延长长度

% 优化螺距
while (max_pitch - min_pitch > tolerance)
    current_pitch = (min_pitch + max_pitch) / 2;

    % 生成轨迹数据
    data = get_trajectory(current_pitch);  % 调用 get_trajectory 生成轨迹
    collision_detected = false;

    % 获取点的数量
    num_points = size(data, 1) / 2;

    % 遍历时间步骤
    for t = 1:size(data, 2)
        x = data(1:2:end, t);  % 奇数行是 x 坐标
        y = data(2:2:end, t);  % 偶数行是 y 坐标
        
        % 清除旧图形
        clf;  
        
        % 绘制原始线段
        plot(x, y, '-o', 'LineWidth', 2, 'MarkerSize', 6);
        hold on;
        
        EX1 = [];
        EY1 = [];
        EX2 = [];
        EY2 = [];

        % 遍历每一对相邻点，计算线段并检测碰撞
        for i = num_points - 1:-1:1
            % 获取线段的端点
            x1 = x(i); y1 = y(i);
            x2 = x(i+1); y2 = y(i+1);

            % 计算线段方向
            dx = x2 - x1;
            dy = y2 - y1;
            length_original = sqrt(dx^2 + dy^2);  % 线段长度
            unit_vec = [dx, dy] / length_original;  % 方向单位向量

            % 计算垂直向量
            perp_vec = [-dy, dx];  % 垂直向量
            perp_unit_vec = perp_vec / norm(perp_vec);  % 垂直单位向量

            % 计算平移后的点
            shift_x = half_width * perp_unit_vec(1);
            shift_y = half_width * perp_unit_vec(2);

            % 平移后的线段
            shifted_x1 = x1 + shift_x; shifted_y1 = y1 + shift_y;
            shifted_x2 = x2 + shift_x; shifted_y2 = y2 + shift_y;

            % 延长平移后的线段
            extended_x1 = shifted_x1 - extension * unit_vec(1);
            extended_y1 = shifted_y1 - extension * unit_vec(2);
            extended_x2 = shifted_x2 + extension * unit_vec(1);
            extended_y2 = shifted_y2 + extension * unit_vec(2);

            EX1(end + 1) = extended_x1;
            EY1(end + 1) = extended_y1;
            EX2(end + 1) = extended_x2;
            EY2(end + 1) = extended_y2;

            % 绘制延长后的线段
            plot([extended_x1, extended_x2], [extended_y1, extended_y2], 'LineWidth', 2);

            % 对于首条线段m的检测，检查m与延长后的线段是否相交
            if i == 1
                m_x1 = [x1, x2] - shift_x;
                m_y1 = [y1, y2] - shift_y;
                extended_m_x1 = m_x1(1) - extension * unit_vec(1);
                extended_m_y1 = m_y1(1) - extension * unit_vec(2);
                extended_m_x2 = m_x1(2) + extension * unit_vec(1);
                extended_m_y2 = m_y1(2) + extension * unit_vec(2);

                % 绘制延长后的线段
                plot([extended_m_x1, extended_m_x2], [extended_m_y1, extended_m_y2], 'LineWidth', 2);

                % 只检测平移并延长的线段
                for j = 2:num_points-1
                    p1 = [EX1(j), EY1(j)];
                    p2 = [EX2(j), EY2(j)];

                    if check_intersection([extended_m_x1, extended_m_y1], [extended_m_x2, extended_m_y2], p1, p2)
                        disp(['碰撞发生在时间: t = ', num2str(t), ' 秒, 对应的螺距: ', num2str(current_pitch)]);
                        collision_detected = true;
                        break;
                    end
                end
            end
        end

        % 检查龙头是否到达中心
        if sqrt(x(1)^2 + y(1)^2) <= center_radius
            disp(['龙头到达中心，距离: ', num2str(sqrt(x(1)^2 + y(1)^2))]);
             disp([' 对应的螺距: ', num2str(current_pitch)]);
            break;
        end

        if collision_detected
            break;
        end
    end

    % 调整螺距
    if collision_detected
        min_pitch = current_pitch;
    else
        max_pitch = current_pitch;
    end
end

disp(['最小螺距为: ', num2str(current_pitch)]);

%% 轨迹生成函数 get_trajectory
function location = get_trajectory(current_pitch)
    % 参数设置
    T = 0:1:500;  % 时间范围
    location = get_location(T, current_pitch);  % 调用 get_location 生成轨迹
end

%% get_location 函数
function result = get_location(T, current_pitch)
    % 基础参数设置
    a = 8.8;
    b = current_pitch / (2 * pi);
    c0 = 3.41 - 2 * 0.275;  % 龙头板的前把手和后把手的固定距离
    c1 = 2.20 - 2 * 0.275;  % 其他板凳的前把手和后把手的固定距离
    v = 1;  % 线速度
    flag = 224;  % 龙身节数

    % 计算角度 theta
    theta = zeros(flag, length(T)+1);
    for i = 2:(length(T)+1)
        theta(1, i) = get_theta0(v, T(i-1), current_pitch);  % 龙头
        theta(2, i) = location_rec(c0, theta(1, i), current_pitch);  % 龙头板
        for j = 2:(flag-1)  % 其余板凳
            theta(j+1, i) = location_rec(c1, theta(j, i), current_pitch);
        end
    end

    % 计算每个板凳的位置
    all_location = zeros(flag*2, length(T));
    for i = 1:length(T)
        for j = 1:flag
            [all_location(2*j-1, i), all_location(2*j, i)] = change(theta(j, i+1), current_pitch);
        end
    end
    
    result = all_location;
end

%% get_theta0 函数
function theta0 = get_theta0(v, t, current_pitch)
    % 计算龙头在 t 时刻的角度 theta0
    a = 8.8;
    b = current_pitch / (2 * pi);
    r = @(theta) a + b * theta;
    dr_dtheta = b;
    L = v * t;  % 龙头走过的弧长
    integrand = @(theta) sqrt((dr_dtheta).^2 + (r(theta)).^2);  % 计算弧长的积分函数
    integral_eq = @(theta0) integral(integrand, theta0, 0) - L;
    theta0 = fzero(integral_eq, 0);  % 求解 theta0
end

%% 极坐标转换为直角坐标
function [x, y] = change(theta, current_pitch)
    a = 8.8;
    b = current_pitch / (2 * pi);
    r = a + b * theta;
    x = r * cos(theta);
    y = r * sin(theta);
end

%% 位置递推函数
function theta_new = location_rec(c, theta, current_pitch)
    % 根据当前点位置计算下一点的角度 theta_new
    a = 8.8;
    b = current_pitch / (2 * pi);
    r = @(theta) a + b * theta;    

    [x_known, y_known] = change(theta, current_pitch);
    eq_fun = @(vars) [
        (vars(2) * cos(vars(1)) - x_known)^2 + ...
        (vars(2) * sin(vars(1)) - y_known)^2 - c^2;
        vars(2) - (a + b * vars(1))
    ];

    result = fsolve(eq_fun, [theta+0.1, r(theta)], optimoptions('fsolve', 'Display', 'off'));
    theta_new = result(1);
end

%% 线段相交检测函数
function intersects = check_intersection(p1, p2, q1, q2)
    % 使用向量叉积法计算线段是否相交
    cross1 = (q1(1) - p1(1)) * (p2(2) - p1(2)) - (q1(2) - p1(2)) * (p2(1) - p1(1));
    cross2 = (q2(1) - p1(1)) * (p2(2) - p1(2)) - (q2(2) - p1(2)) * (p2(1) - p1(1));
    cross3 = (p1(1) - q1(1)) * (q2(2) - q1(2)) - (p1(2) - q1(2)) * (q2(1) - q1(1));
    cross4 = (p2(1) - q1(1)) * (q2(2) - q1(2)) - (p2(2) - q1(2)) * (q2(1) - q1(1));

    % 判断叉积符号是否相异，若相异则两线段相交
    intersects = (cross1 * cross2 < 0) && (cross3 * cross4 < 0);
end
