clear; clc; close all;

%AADMM

% 初始化参数
k_max = 4;  % 最大迭代次数
x_t = zeros(k_max, 1);    % tilde{x}^k
y_t = zeros(k_max, 1);    % tilde{y}^k
lambda_t = zeros(k_max, 1); % tilde{lambda}^k
x = zeros(k_max, 1);          % x^k
y = zeros(k_max, 1);          % y^k
lambda = zeros(k_max, 1);     % lambda^k
t = zeros(k_max, 1);          % t_k

%nx
nx = zeros(k_max, 1);          % x^k
ny = zeros(k_max, 1);          % y^k
nlambda = zeros(k_max, 1);

% a=55;
% b=64;
% c=74;

a=60;
b=-14;
c=83;

% % 生成-5到5之间的随机整数
% a = randi([-100, 100]);
% b = randi([-100, 100]);
% c = randi([-100, 100]);

% 设置初值
x(1) = a;
y(1) = b;
lambda(1) = c;
t(1) = 1;

x_t0 = x(1);      % tilde{x}^0 = x^1
y_t0 = y(1);      % tilde{y}^0 = y^1
lambda_t0 = lambda(1); % tilde{lambda}^0 = lambda^1

%nx
nx(1)=x(1);
ny(1)=y(1);
nlambda(1)=lambda(1);

% k=1
k=1;
x_t(k) = -2*y(k) + lambda(k) + 1;
lambda_t(k) = lambda(k) - (1/2*x_t(k) + y(k) - 3/2);
y_t(k) = 1/2 * lambda_t(k);

t(k+1) = (1 + sqrt(1 + 4*t(k)^2)) / 2;

x(k+1) = x_t(k);
y(k+1) = (t(k+1)-1)/(2*t(k+1)) * y_t(k) - ...
    (t(k)-1)/(2*t(k+1)) * y_t0 + ...
    t(k)/(2*t(k+1)) * (lambda_t(k) - lambda(k) + y(k));
lambda(k+1) = lambda_t(k) + ...
    (t(k)-1)/t(k+1) * (lambda_t(k) - lambda_t0) + ...
    t(k)/t(k+1) * (lambda_t(k) - lambda(k) + y(k) - y_t(k)) - ...
    y(k+1);

for k = 2:k_max
    
    % 预测
    x_t(k) = -2*y(k) + lambda(k) + 1;
    lambda_t(k) = lambda(k) - (1/2*x_t(k) + y(k) - 3/2);
    y_t(k) = 1/2 * lambda_t(k);
    
    % 更新t_k+1
    t(k+1) = (1 + sqrt(1 + 4*t(k)^2)) / 2;
    
    % 校正
    x(k+1) = x_t(k);
    y(k+1) = (t(k+1)-1)/(2*t(k+1)) * y_t(k) - ...
        (t(k)-1)/(2*t(k+1)) * y_t(k-1) + ...
        t(k)/(2*t(k+1)) * (lambda_t(k) - lambda(k) + y(k));
    lambda(k+1) = lambda_t(k) + ...
        (t(k)-1)/t(k+1) * (lambda_t(k) - lambda_t(k-1)) + ...
        t(k)/t(k+1) * (lambda_t(k) - lambda(k) + y(k) - y_t(k)) - ...
        y(k+1);
end

% nx
for k=1:k_max
    nx(k+1)=-2*ny(k)+nlambda(k)+1;
    nlambda(k+1)=nlambda(k)-(1/2*nx(k+1)+ny(k)-3/2);
    ny(k+1)=1/2*nlambda(k+1);
end

tx=[x_t0;x_t];
ty=[y_t0;y_t];
tlambda=[lambda_t0;lambda_t];

% 绘制三维图像
figure;
set(gcf, 'Position', [100, 100, 1200, 800]);

% 三维散点图
subplot(2,3,[1,2,4,5]);

% 添加绿色点(1,1,2)
scatter3(1, 1, 2, 80, 'filled', 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
hold on;

% 添加AADMM
scatter3(tx, ty, tlambda, 80, 'filled', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
%添加无加速点
scatter3(nx, ny, nlambda, 80, 'filled', 'MarkerFaceColor', 'b', 'HandleVisibility', 'off');

scatter3(x(2), y(2), lambda(2), 80, 'filled', 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');


% 添加蓝色向量：从蓝色标注点1指向蓝色标注点2
if length(nx) >= 2
    quiver3(nx(2), ny(2), nlambda(2), ...
            nx(3)-nx(2), ny(3)-ny(2), nlambda(3)-nlambda(2), ...
            0, 'LineWidth', 2, 'Color', 'b', 'MaxHeadSize', 0.2, ...
            'HandleVisibility', 'off');
end

% 添加红色箭头1：从蓝色标注点2指向绿色点
if length(nx) >= 2
    quiver3(nx(2), ny(2), nlambda(2), ...
            x(2)-nx(2), y(2)-ny(2), lambda(2)-nlambda(2), ...
            0, 'LineWidth', 2, 'Color', 'r', 'MaxHeadSize', 0.3, ...
            'HandleVisibility', 'off');
end


% 添加红色箭头2：从绿色点指向红色标注点3
if length(tx) >= 3
    quiver3(x(2), y(2), lambda(2), ...
            tx(3)-x(2), ty(3)-y(2), tlambda(3)-lambda(2), ...
            0, 'LineWidth', 2, 'Color', 'r', 'MaxHeadSize', 0.2, ...
            'HandleVisibility', 'off');
end

% 添加图例并设置LaTeX解释器
leg = legend('Location', 'best');
set(leg, 'Interpreter', 'latex');

% 为红色点添加标注 \tilde{x}^{0} ...
for k = 1:length(tx)
    if k == 1
        % 第一个红色点标注为O
        text(tx(k), ty(k), tlambda(k)+0.1, '$O$', ...
            'Interpreter', 'latex', 'FontSize', 20, 'Color', 'red', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    elseif k == 2
        % 第二个红色点标注为A
        text(tx(k), ty(k), tlambda(k)+0.1, '$A$', ...
            'Interpreter', 'latex', 'FontSize', 20, 'Color', 'red', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    elseif k == 3
        % 第三个红色点标注为D
        text(tx(k), ty(k), tlambda(k)+2, '$F_{2}$', ...
            'Interpreter', 'latex', 'FontSize', 20, 'Color', 'red', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    else
        text(tx(k), ty(k), tlambda(k)+2, sprintf('$F_{%d}$', k-1), ...
            'Interpreter', 'latex', 'FontSize', 20, 'Color', 'red', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

% 为蓝色点添加标注
for k = 1:length(nx)
    if k == 1
        % 第一个蓝色点标注为O
        text(nx(k), ny(k), nlambda(k)+0.1, '$O$', ...
            'Interpreter', 'latex', 'FontSize', 20, 'Color', 'blue', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    elseif k == 2
        % 第二个蓝色点标注为A（左侧）
        text(nx(k), ny(k), nlambda(k)+3 , '$A$', ...
            'Interpreter', 'latex', 'FontSize', 20, 'Color', 'blue', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    elseif k == 3
        % 第三个蓝色点标注为B（左侧）
        text(nx(k), ny(k), nlambda(k)+2, '$P_{2}$', ...
            'Interpreter', 'latex', 'FontSize', 20, 'Color', 'blue', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    else
        % 其他蓝色点保持原来的数字标注（左侧）
        text(nx(k), ny(k), nlambda(k)+2, sprintf('$P_{%d}$', k-1), ...
            'Interpreter', 'latex', 'FontSize', 20, 'Color', 'blue', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    end
end

% 为黑色点添加标注B
text(x(2), y(2), lambda(2)+0.1, '$B$', ...
    'Interpreter', 'latex', 'FontSize', 12, 'Color', 'black', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

% 添加绿色点的图例标注（统一大小）
plot3(NaN, NaN, NaN, 'o', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'g', ...      % 填充绿色
    'MarkerEdgeColor', 'g', ...      % 边缘也设为绿色（关键修改）
    'DisplayName', '$(x^*,y^*,\lambda^*)$', 'LineStyle', 'none');% 添加红色点的图例标注（统一大小）
% 添加红色点的图例标注
plot3(NaN, NaN, NaN, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'Full prediction-correction');
% 添加蓝色点的图例标注
plot3(NaN, NaN, NaN, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', ...
    'DisplayName', 'Prediction-only', 'LineStyle', 'none');


% 添加标签
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 14);
grid on;


% 设置坐标轴刻度为4个整数
x_min = floor(min([tx; nx; 1]));
x_max = ceil(max([tx; nx; 1]));
y_min = floor(min([ty; ny; 1]));
y_max = ceil(max([ty; ny; 1]));
z_min = floor(min([tlambda; nlambda; 2]));
z_max = ceil(max([tlambda; nlambda; 2]));

% 生成4个整数刻度
x_ticks = unique(round(linspace(x_min, x_max, 4)));
y_ticks = unique(round(linspace(y_min, y_max, 4)));
z_ticks = unique(round(linspace(z_min, z_max, 4)));

% 确保至少有2个刻度
if length(x_ticks) < 2
    x_ticks = [x_min, x_max];
end
if length(y_ticks) < 2
    y_ticks = [y_min, y_max];
end
if length(z_ticks) < 2
    z_ticks = [z_min, z_max];
end

xticks(x_ticks);
yticks(y_ticks);
zticks(z_ticks);


% 添加图例并设置LaTeX解释器
leg = legend('Location', 'best');
set(leg, 'Interpreter', 'latex');

% 确保图形显示
drawnow;

% 方法1：调整图形尺寸后导出
set(gcf, 'Position', [100, 100, 2000, 1500]); % 增大窗口尺寸

% 方法2：使用print函数导出高质量EPS
print('-depsc', '-r300', '-painters', 'E1.eps');
