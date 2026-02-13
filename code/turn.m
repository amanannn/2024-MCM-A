function bench_dragon_turn_schematic
    % Setting up the figure
    figure('Color', 'w', 'Position', [100, 100, 800, 700]);
    hold on; axis equal; grid on;
    xlabel('x (m)'); ylabel('y (m)');
    title('问题4: 板凳龙S形调头路径示意图 (2:1圆弧)', 'FontSize', 14);

    %% 1. 参数设置
    p = 1.7;                % 螺距 (Pitch) = 1.7m
    b = p / (2 * pi);       % 螺线参数 r = b * theta
    R_zone = 4.5;           % 调头空间半径 (直径9m)
    
    % 绘制调头空间边界 (黄色区域示意)
    rectangle('Position', [-R_zone, -R_zone, 2*R_zone, 2*R_zone], ...
              'Curvature', [1, 1], 'EdgeColor', [0.9, 0.8, 0], ...
              'LineStyle', '--', 'LineWidth', 2);
    fill(R_zone*cos(0:0.01:2*pi), R_zone*sin(0:0.01:2*pi), ...
         [1, 1, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    %% 2. 求解S形曲线参数
    % 我们需要找到一个切点 theta_tan，在该点：
    % 1. 螺线切线与圆弧切线重合
    % 2. 圆弧1半径 r1 = 2 * r2
    % 3. 圆弧1和圆弧2在原点对称处相连 (基于中心对称假设)
    
    % 使用数值优化寻找最佳的 theta (螺线截断角度)
    initial_theta = 2 * pi; 
    options = optimset('Display', 'off', 'TolX', 1e-8);
    [theta_opt, fval] = fminsearch(@(t) solve_geometry(t, b), initial_theta, options);
    
    % 获取最优解下的几何参数
    [cost, r2, center1, P_in, start_angle1, end_angle1] = solve_geometry(theta_opt, b);
    r1 = 2 * r2;

    %% 3. 生成路径数据
    
    % --- A. 盘入螺线 (Spiral In) ---
    % 从外部画到切点 P_in
    theta_max_plot = 10 * pi; 
    theta_spiral_in = linspace(theta_opt, theta_max_plot, 500);
    x_spiral_in = b * theta_spiral_in .* cos(theta_spiral_in);
    y_spiral_in = b * theta_spiral_in .* sin(theta_spiral_in);
    
    % --- B. 盘出螺线 (Spiral Out) ---
    % 中心对称: (-x, -y)
    x_spiral_out = -x_spiral_in;
    y_spiral_out = -y_spiral_in;
    
    % --- C. S形曲线 - 第一段圆弧 (半径 2r) ---
    % 生成从切点到中间交点的圆弧
    angles1 = linspace(start_angle1, end_angle1, 100);
    x_arc1 = center1(1) + r1 * cos(angles1);
    y_arc1 = center1(2) + r1 * sin(angles1);
    
    % --- D. S形曲线 - 第二段圆弧 (半径 r) ---
    % 中心对称，圆心为 -center1，半径 r2
    center2 = -center1;
    % 对应的角度范围也需要旋转180度 (pi)
    angles2 = linspace(start_angle1 + pi, end_angle1 + pi, 100);
    x_arc2 = center2(1) + r2 * cos(angles2);
    y_arc2 = center2(2) + r2 * sin(angles2);

    %% 4. 绘图
    
    % 绘制螺线
    p1 = plot(x_spiral_in, y_spiral_in, 'b-', 'LineWidth', 1.5); % 盘入
    p4 = plot(x_spiral_out, y_spiral_out, 'r-', 'LineWidth', 1.5); % 盘出
    
    % 绘制S形圆弧
    p2 = plot(x_arc1, y_arc1, 'g-', 'LineWidth', 3); % 圆弧1 (大)
    p3 = plot(x_arc2, y_arc2, 'm-', 'LineWidth', 3); % 圆弧2 (小)
    
    % 标记关键点
    plot(P_in(1), P_in(2), 'ko', 'MarkerFaceColor', 'k'); % 切点1
    plot(-P_in(1), -P_in(2), 'ko', 'MarkerFaceColor', 'k'); % 切点2
    
    % 计算并标记中间连接点 (Arc1 和 Arc2 的切点)
    mid_point = (center1 * r2 + center2 * r1) / (r1 + r2);
    plot(mid_point(1), mid_point(2), 'ks', 'MarkerFaceColor', 'w'); 

    %% 5. 图例和标注
    legend([p1, p2, p3, p4], ...
           '盘入螺线 (Spiral In)', ...
           '第1段圆弧 (Radius = 2R)', ...
           '第2段圆弧 (Radius = R)', ...
           '盘出螺线 (Spiral Out)', ...
           'Location', 'best');
       
    text(P_in(1)+0.2, P_in(2), '切点 A', 'FontSize', 10);
    text(mid_point(1)+0.2, mid_point(2), '反向切点 B', 'FontSize', 10);
    
    % 标注半径关系
    text(center1(1), center1(2), 'O_1', 'HorizontalAlignment', 'center');
    plot([center1(1), P_in(1)], [center1(2), P_in(2)], 'k--', 'LineWidth', 0.5);
    
    text(center2(1), center2(2), 'O_2', 'HorizontalAlignment', 'center');
    
    xlim([-5, 5]); ylim([-5, 5]);
end

%% 辅助函数：几何解算
function [cost, r2, center1, P_in, start_angle, end_angle] = solve_geometry(theta, b)
    % 1. 计算螺线上切点 P 的坐标
    r = b * theta;
    x = r * cos(theta);
    y = r * sin(theta);
    P_in = [x, y];
    
    % 2. 计算螺线在该点的切向量和法向量
    dx = b * (cos(theta) - theta * sin(theta));
    dy = b * (sin(theta) + theta * cos(theta));
    
    tangent_vec = [dx, dy];
    tangent_vec = tangent_vec / norm(tangent_vec);
    
    % 法向量 (指向螺线内侧)
    normal_vec = [-tangent_vec(2), tangent_vec(1)];
    
    % 3. 确定圆心 O1
    % 几何约束推导：
    % 设小圆弧半径为 r2，大圆弧半径为 r1 = 2*r2。
    % 整个图形关于原点中心对称。
    % 联立方程求解 r2
    P_dot_N = dot(P_in, normal_vec);
    P_sq = dot(P_in, P_in);
    
    % 一元二次方程系数: 1.75*r2^2 + 4*(P·N)*r2 + |P|^2 = 0
    roots_r = roots([1.75, 4 * P_dot_N, P_sq]);
    
    % 取正实数解
    valid_r = roots_r(roots_r > 0 & imag(roots_r) == 0);
    
    if isempty(valid_r)
        cost = 1e6; 
        r2 = 1; center1 = [0,0]; start_angle=0; end_angle=0;
        return;
    end
    
    r2 = min(valid_r); 
    r1 = 2 * r2;
    
    center1 = P_in + normal_vec * r1;
    cost = 0; 
    
    % 计算绘图用的角度
    vec_start = P_in - center1;
    start_angle = atan2(vec_start(2), vec_start(1));
    
    center2 = -center1;
    vec_end = center2 - center1;
    end_angle = atan2(vec_end(2), vec_end(1));
    
    % 角度调整
    if end_angle < start_angle
       % end_angle = end_angle + 2*pi;
    end
end