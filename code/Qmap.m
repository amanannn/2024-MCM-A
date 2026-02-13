%% 板凳龙“盘入”示意图（龙头在螺线内侧，龙尾在外侧）
% 采用等距螺线 r = b*theta (螺距0.55m)，按相邻把手孔间距布点。
clear; clc; close all;

%% 参数设定
pitch = 0.55;    % 螺距
b = pitch/(2*pi); % 螺线系数
W = 0.30;        % 板宽
a = 0.275;       % 孔中心到板头距离
L_head = 3.41;   % 龙头板长
L_body = 2.20;   % 其他板长
N_bench = 223;   % 板凳数
N_pts = N_bench + 1; % 把手点数

d = [L_head - 2*a; repmat(L_body - 2*a, N_bench-1, 1)]; % 相邻把手孔间距

%% 设置龙头内侧角度
theta_head = 2*pi; 

%% 等距螺线弧长计算
arcS = @(th) (b/2) .* ( th.*sqrt(th.^2 + 1) + asinh(th) );

%% 放置把手点
theta = zeros(N_pts,1);
theta(1) = theta_head;

for k = 1:N_bench
    ds = d(k);
    target = arcS(theta(k)) + ds;
    dth0 = ds/(b*sqrt(theta(k)^2+1));
    th_guess = theta(k) + dth0;
    f = @(th) arcS(th) - target;
    theta(k+1) = fzero(f, th_guess);
end

r = b*theta;
P = [r.*cos(theta), r.*sin(theta)];

%% 绘制背景螺线
figure('Color','w'); hold on; axis equal; grid on;
title('板凳龙盘入示意图');

theta_bg = linspace(0, max(theta)*1.02, 6000);
r_bg = b*theta_bg;
plot(r_bg.*cos(theta_bg), r_bg.*sin(theta_bg), '-', 'LineWidth', 1);

%% 绘制板凳
for i = 1:N_bench
    pF = P(i,:); pR = P(i+1,:);
    v = pR - pF; nv = norm(v);
    if nv < 1e-10, continue; end
    u = v / nv; n = [-u(2), u(1)];
    eF = pF - a*u; eR = pR + a*u;
    c1 = eF + (W/2)*n; c2 = eF - (W/2)*n; c3 = eR - (W/2)*n; c4 = eR + (W/2)*n;
    patch([c1(1),c2(1),c3(1),c4(1)], [c1(2),c2(2),c3(2),c4(2)], ...
          0.8, 'FaceAlpha', 0.25, 'EdgeColor', 'k', 'LineWidth', 0.6);
end

%% 标记与标注
plot(P(:,1), P(:,2), '.', 'MarkerSize', 12);
text(P(1,1), P(1,2), "  龙头前把手", 'FontSize', 10);
text(P(end,1), P(end,2), "  龙尾后把手", 'FontSize', 10);
xlabel('x / m'); ylabel('y / m');

% 计算背景螺线坐标（确保维度一致）
x_bg = r_bg .* cos(theta_bg);
y_bg = r_bg .* sin(theta_bg);

% 自动调整视野
pad = 0.6;
xmin = min([P(:,1); x_bg(:)]) - pad;
xmax = max([P(:,1); x_bg(:)]) + pad;
ymin = min([P(:,2); y_bg(:)]) - pad;
ymax = max([P(:,2); y_bg(:)]) + pad;
xlim([xmin, xmax]); ylim([ymin, ymax]);

exportgraphics(gcf, 'bench_dragon_coil_head_inside.png', 'Resolution', 300);