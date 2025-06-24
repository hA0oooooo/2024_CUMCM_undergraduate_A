% 读取 Excel 数据
filename = 'result4.xlsx';  % Excel 文件名
data = readmatrix(filename);  % 读取 Excel 文件为矩阵

% 设置平移的距离（half_width）和延长长度（extension）
half_width = 0.15;  % 设定平移的距离，可以根据实际调整
extension = 0.275;    % 设定延长长度

% 获取点的数量
global num_points;
num_points = size(data, 1) / 2;  % 假设是2N行，每个点有 x 和 y 坐标

% 创建动画
figure;



for t = 90:size(data, 2)  % 遍历每个时刻
    % 提取当前时刻的 x 和 y 坐标
    x = data(1:2:end, t);  % 奇数行是 x 坐标
    y = data(2:2:end, t);  % 偶数行是 y 坐标
    
    
    % 原始线段绘制
    plot(x, y, '-o', 'LineWidth', 2, 'MarkerSize', 6);
    hold on;
    text(0, 11, num_points, ['Frame: ', num2str(t)], 'FontSize', 12);
%%    
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
    
    % plot(X0, Y0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
%%

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
        
        % % 对于首条线段m的检测，检查m与延长后的线段是否相交
        % if i == 1
        %     m_x1 = [x1, x2] - shift_x;
        %     m_y1 = [y1, y2] - shift_y;
        %     extended_m_x1 = m_x1(1) - extension * unit_vec(1);
        %     extended_m_y1 = m_y1(1) - extension * unit_vec(2);
        %     extended_m_x2 = m_x1(2) + extension * unit_vec(1);
        %     extended_m_y2 = m_y1(2) + extension * unit_vec(2);
        %     % 绘制延长后的线段
        %     plot([extended_m_x1, extended_m_x2], [extended_m_y1, extended_m_y2], 'LineWidth', 2);
        % 
        %     % 只检测平移并延长的线段
        %     for j = 2:num_points-1
        %         p1 = [EX1(j), EY1(j)];
        %         p2 = [EX2(j), EY2(j)];
        % 
        %         if check_intersection([extended_m_x1, extended_m_y1], [extended_m_x2, extended_m_y2], p1, p2)
        %             disp(['碰撞发生在时刻: t = ', num2str(t)]);
        %             % disp(['此时位置为：', get_location(t)]);
        %             return;  % 停止动画
        %         end
        %     end
        % end
    end
    
    axis equal;  % 设置 x 和 y 轴比例相同
    axis([min(data(:)) max(data(:)) min(data(:)) max(data(:))]);  % 设置动态坐标轴范围
    grid on;


    % 控制动画速度
    pause(0.1);  % 暂停 0.1 秒
    clf;  % 清除旧图形
end


%% 定义交点检测函数，考虑延长后的线段
function intersects = check_intersection(p1, p2, q1, q2)
    % 使用向量叉积法计算线段是否相交
    cross1 = (q1(1) - p1(1)) * (p2(2) - p1(2)) - (q1(2) - p1(2)) * (p2(1) - p1(1));
    cross2 = (q2(1) - p1(1)) * (p2(2) - p1(2)) - (q2(2) - p1(2)) * (p2(1) - p1(1));
    cross3 = (p1(1) - q1(1)) * (q2(2) - q1(2)) - (p1(2) - q1(2)) * (q2(1) - q1(1));
    cross4 = (p2(1) - q1(1)) * (q2(2) - q1(2)) - (p2(2) - q1(2)) * (q2(1) - q1(1));
    
    % 判断叉积符号是否相异，若相异则两线段相交
    intersects = (cross1 * cross2 < 0) && (cross3 * cross4 < 0);
end

% %% 得到位置的函数
% function result = get_location(t)
%     global p;
% 
%     % 定义基本参数及函数
%     a = 16 * 0.55;
%     b = p / (2 * pi);
% 
%     % 龙头的线速度
%     v = 1;
% 
%     % 计算时间t时刻龙头的theta角度
%     theta0 = get_theta0(v, t);
% 
%     % 根据螺线公式，计算直角坐标位置
%     [x, y] = change(theta0);
% 
%     % 返回结果为一个2x1的向量，分别表示x和y坐标
%     result = [x; y];
% end
% 
% %% 计算经过时间t龙头的theta位置函数
% function theta0 = get_theta0(v, t)
%     global p;
%     % 等距螺线基本参数
%     a = 16 * 0.55;
%     b = p / (2 * pi);
% 
%     % 定义 r(theta)
%     r = @(theta) a + b * theta;
% 
%     % 它的导数 dr/dtheta
%     dr_dtheta = @(theta) b;
% 
%     % 龙头在t时间内走过的弧长
%     L = v * t;
% 
%     % 定义被积函数
%     integrand = @(theta) sqrt((dr_dtheta(theta)).^2 + (r(theta)).^2);
% 
%     % 定义积分方程：目标是求解 theta0，使得积分结果等于走过的弧长 L
%     integral_eq = @(theta0) integral(integrand, theta0, 0) - L;
% 
%     % 使用 fsolve 求解 theta0
%     theta0_initial_guess = 0; 
%     theta0 = fzero(integral_eq, theta0_initial_guess);
% end
% 
% %% 将theta转化为直角坐标系位置的函数
% function [x, y] = change(theta)
%     global p;
% 
%     % 定义螺线的基本参数
%     a = 16 * 0.55;
%     b = p / (2 * pi);
% 
%     % 计算螺线半径
%     r = a + b * theta;
% 
%     % 转换为直角坐标
%     x = r * cos(theta);
%     y = r * sin(theta);
% end
