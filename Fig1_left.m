clear; clc; close all;

%AADMM

% 初始化参数
k_max = 5;  % 最大迭代次数
x_t = zeros(k_max, 1);    % tilde{x}^k
y_t = zeros(k_max, 1);    % tilde{y}^k
lambda_t = zeros(k_max, 1); % tilde{lambda}^k
x = zeros(k_max, 1);          % x^k
y = zeros(k_max, 1);          % y^k
lambda = zeros(k_max, 1);     % lambda^k
t = zeros(k_max, 1);          % t_k

% 初值
a=-95;
b=86;
c=46;

% 设置初值
x(1) = a;
y(1) = b;
lambda(1) = c;
t(1) = 1;

x_t0 = x(1);      % tilde{x}^0 = x^1
y_t0 = y(1);      % tilde{y}^0 = y^1  
lambda_t0 = lambda(1); % tilde{lambda}^0 = lambda^1

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

tx=[x_t0;x_t];
ty=[y_t0;y_t];
tlambda=[lambda_t0;lambda_t];

% ADMM

% 重新初始化ADMM变量
x_admm = zeros(k_max+1, 1);
y_admm = zeros(k_max+1, 1);
lambda_admm = zeros(k_max+1, 1);

% 设置初值
x_admm(1) = a;
y_admm(1) = b;
lambda_admm(1) = c;

for k = 1:k_max
    % x更新
    x_admm(k+1) = -2*y_admm(k) + lambda_admm(k) + 1;
    
    % y更新
    y_admm(k+1) = -1/6*x_admm(k+1) + 1/3*lambda_admm(k) + 1/2;
    
    % lambda更新
    lambda_admm(k+1) = lambda_admm(k) - (1/2*x_admm(k+1) + y_admm(k+1) - 3/2);
end

% 绘制三维图像
figure;
set(gcf, 'Position', [100, 100, 1200, 800]);

% 三维散点图
subplot(2,3,[1,2,4,5]);

% 添加绿色点(1,1,2)
scatter3(1, 1, 2, 80, 'filled', 'MarkerFaceColor', 'g', 'DisplayName', '$(x^*,y^*,\lambda^*)$');
hold on;

% AADMM - 红色
scatter3(tx, ty, tlambda, 80, 'filled', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
hold on;
plot3(tx, ty, tlambda, 'r--', 'LineWidth', 1.5,'DisplayName', 'AADMM');

% ADMM - 蓝色
scatter3(x_admm, y_admm, lambda_admm, 80, 'filled', 'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
plot3(x_admm, y_admm, lambda_admm, 'b-', 'LineWidth', 1.5,'DisplayName', 'ADMM');

% 为AADMM点添加数字标号（从0开始）
for i = 1:length(tx)
    if i == 5  % 序号4（因为从0开始，所以i=5对应序号4）
        text(tx(i)+3, ty(i), tlambda(i), num2str(i-1), ...
             'FontSize', 20, 'FontWeight', 'bold', 'Color', 'r', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    elseif i == 6  % 序号5（因为从0开始，所以i=6对应序号5）
        text(tx(i), ty(i)-1, tlambda(i), num2str(i-1), ...
             'FontSize', 20, 'FontWeight', 'bold', 'Color', 'r', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    else
        text(tx(i), ty(i), tlambda(i)+3, num2str(i-1), ...
             'FontSize', 20, 'FontWeight', 'bold', 'Color', 'r', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

% 为ADMM点添加数字标号（从0开始）
for i = 1:length(x_admm)
    text(x_admm(i), y_admm(i), lambda_admm(i)+3, num2str(i-1), ...
         'FontSize', 20, 'FontWeight', 'bold', 'Color', 'b', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% 添加标签
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 14);
grid on;

% 设置坐标轴刻度，减少数字密度并设置为整数
x_ticks = linspace(min([tx; x_admm; 1]), max([tx; x_admm; 1]), 4); % 4个刻度
y_ticks = linspace(min([ty; y_admm; 1]), max([ty; y_admm; 1]), 4); % 4个刻度
z_ticks = linspace(min([tlambda; lambda_admm; 2]), max([tlambda; lambda_admm; 2]), 4); % 4个刻度

% 将刻度值取整
x_ticks = round(x_ticks);
y_ticks = round(y_ticks);
z_ticks = round(z_ticks);

xticks(x_ticks);
yticks(y_ticks);
zticks(z_ticks);

% 添加图例并设置LaTeX解释器
leg = legend('Location', 'best');
set(leg, 'Interpreter', 'latex');

% 显示数值结果
fprintf('AADMM迭代结果:\n');
fprintf('初始点: x=%.4f, y=%.4f, lambda=%.4f\n', tx(1), ty(1),tlambda(1));
fprintf('最终点: x=%.4f, y=%.4f, lambda=%.4f\n', tx(end), ty(end),tlambda(end));

fprintf('\nADMM迭代结果:\n');
fprintf('初始点: x=%.4f, y=%.4f, lambda=%.4f\n', x_admm(1), y_admm(1), lambda_admm(1));
fprintf('最终点: x=%.4f, y=%.4f, lambda=%.4f\n', x_admm(end), y_admm(end), lambda_admm(end));

% 确保图形显示
drawnow;

% 方法1：调整图形尺寸后导出
set(gcf, 'Position', [100, 100, 2000, 1500]); % 增大窗口尺寸

% 方法2：使用print函数导出高质量EPS
% print('-depsc', '-r300', '-painters', 'compare.eps');
print('-djpeg', '-r150', 'output_image.jpg');




