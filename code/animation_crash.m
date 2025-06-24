% 读取 Excel 数据
filename = 'result1_500.xlsx';  % Excel 文件名
data = readmatrix(filename);  % 读取 Excel 文件为矩阵

% 设置平移的距离（half_width）和延长长度（extension）
half_width = 0.15;  % 设定平移的距离，可以根据实际调整
extension = 0.275;  % 设定延长长度

% 获取点的数量
num_points = size(data, 1) / 2;  % 假设是2N行，每个点有 x 和 y 坐标

% 创建动画
figure;

for t = 1:size(data, 2)  % 遍历每个时刻
    % 提取当前时刻的 x 和 y 坐标
    x = data(1:2:end, t);  % 奇数行是 x 坐标
    y = data(2:2:end, t);  % 偶数行是 y 坐标
    
    clf;  % 清除旧图形
    
    % 原始线段绘制
    plot(x, y, '-o', 'LineWidth', 2, 'MarkerSize', 6);
    hold on;
    
    % 添加显示时间的文本
    time_str = ['Time: ', num2str(t), ' s'];
    text(min(x), max(y), time_str, 'FontSize', 12, 'Color', 'r', 'FontWeight', 'bold');
    
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
                    disp(['碰撞发生在时刻: t = ', num2str(t)]);
                    return;  % 停止动画
                end
            end
        end
    end
    
     % 设置 x 和 y 轴比例相同
    axis([min(data(:)) max(data(:)) min(data(:)) max(data(:))]);  % 设置动态坐标轴范围
    grid on;
    
    % 控制动画速度
    pause(0.1);  % 暂停 0.1 秒
end

% 定义交点检测函数，考虑延长后的线段
function intersects = check_intersection(p1, p2, q1, q2)
    % 使用向量叉积法计算线段是否相交
    cross1 = (q1(1) - p1(1)) * (p2(2) - p1(2)) - (q1(2) - p1(2)) * (p2(1) - p1(1));
    cross2 = (q2(1) - p1(1)) * (p2(2) - p1(2)) - (q2(2) - p1(2)) * (p2(1) - p1(1));
    cross3 = (p1(1) - q1(1)) * (q2(2) - q1(2)) - (p1(2) - q1(2)) * (q2(1) - q1(1));
    cross4 = (p2(1) - q1(1)) * (q2(2) - q1(2)) - (p2(2) - q1(2)) * (q2(1) - q1(1));
    
    % 判断叉积符号是否相异，若相异则两线段相交
    intersects = (cross1 * cross2 < 0) && (cross3 * cross4 < 0);
end
