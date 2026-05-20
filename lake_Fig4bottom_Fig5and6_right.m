close all;clc;clear all;

Ur=double(imread('lake.png'))/255;
ur=double(imread('noise_lake.png'))/255;

epsilon = 1e-4;

[rows, cols]=size(ur);
alpha = 0.1;
steps=0.2;
beta = 1;
t=1;
numIter = 5e4;
imgSz = rows*cols;
ne=1;

% AADMM
q =zeros(rows,cols);
Q=q;
tq=q;
Tq=q;

uu = ur;
tu=uu;

pp1 = zeros(rows, cols+1);
pp2 = zeros(rows+1, cols);

X1=[];
Y1=[];
AAresidual=[];  % residual abot AADMM  ||divp~+q~-f||_2

for i = 1:numIter
    
    %sub-step1: update tp^(k+1)=p^{k+1}
    for j=1:ne
        divp = pp1(:,2:cols+1)-pp1(:,1:cols)+pp2(2:rows+1,:)-pp2(1:rows,:);
        pt = divp +q- (ur+uu/beta);
        pp1(:,2:cols) = steps*(pt(:,2:cols) - pt(:,1:cols-1)) + pp1(:,2:cols);
        pp2(2:rows,:) = steps*(pt(2:rows,:) - pt(1:rows-1,:)) + pp2(2:rows,:);
    end
    
    % projection steps
    pp1 = min(pp1,alpha);
    pp1 = max(pp1,-alpha);
    pp2 = min(pp2,alpha);
    pp2 = max(pp2,-alpha);
    
    % update divp
    divp = pp1(:,2:cols+1)-pp1(:,1:cols)+pp2(2:rows+1,:)-pp2(1:rows,:);
    
    %sub-step2 update tu
    %keep tu
    Tu=tu;
    tu = uu + beta*(ur-divp-q);
    
    %sub-step update tq
    %keep tq
    Tq=tq;
    tq = tu;
    
    % residual ||divp~+q~-f||_2
    residual = norm(divp+tq-ur, 2);
    AAresidual(i) = residual;
    
    %sub-step4 update tk
    % keep tk
    T=t;
    t=(1+sqrt(1+4*t^2))/2;
    
    %sub-step5 update p^{k+1} see sub-step1
    
    %sub-step6 update q^{k+1}
    %keep q
    Q=q;
    q=(t-1)/(2*t)*tq-(T-1)/(2*t)*Tq+T/(2*t)*(tu-uu+(1-beta)*tq+beta*q);
    
    %sub-step7 update u^{k+1}
    % keep uu
    ku=uu;
    uu=(t-1)/(2*t)*(2*tu-beta*tq)+(T-1)/(2*t)*(beta*Tq-2*Tu)+...
        T/(2*t)*((beta-2)*uu+(4-beta)*tu+(beta^2-3*beta)*tq+(2*beta-beta^2)*Q);
    
    % ||uk+1-uk||
    Aerru =sqrt(sum(sum((ku-uu).*(ku-uu))));
    
    X1(i)=(i);
    Y1(i)=log(Aerru);
    
    if Aerru < epsilon
        break;
    end
    
end

% keep tu
u1=tu;

% average error of tu and Ur
errtuUr=sum(sum(abs(tu-Ur)))/imgSz;
msg1 = sprintf('number of AADMM = %u. \n', i);
disp(msg1);
msg2 = sprintf('average error of tu and Ur = %u. \n', errtuUr);
disp(msg2);
msg3 = sprintf('||uk+1-uk|| = %u. \n', Aerru);
disp(msg3);


% ADMM
X2=[];
Y2=[];
Aresidual=[];  % residual abot ADMM  ||divp^k+1+q^k+1-f||_2

q =zeros(rows,cols);
uu = ur;

pp1 =zeros(rows, cols+1);
pp2 =zeros(rows+1, cols);
divp = pp1(:,2:cols+1)-pp1(:,1:cols)+pp2(2:rows+1,:)-pp2(1:rows,:);

for i = 1:numIter
    
    %sub-step1: update p^(k+1)
    for j=1:ne
        divp = pp1(:,2:cols+1)-pp1(:,1:cols)+pp2(2:rows+1,:)-pp2(1:rows,:);
        pt = divp +q- (ur+uu/beta);
        pp1(:,2:cols) = steps*(pt(:,2:cols) - pt(:,1:cols-1)) + pp1(:,2:cols);
        pp2(2:rows,:) = steps*(pt(2:rows,:) - pt(1:rows-1,:)) + pp2(2:rows,:);
    end
    
    % projection steps
    pp1 = min(pp1,alpha);
    pp1 = max(pp1,-alpha);
    pp2 = min(pp2,alpha);
    pp2 = max(pp2,-alpha);
    
    % update divp
    Divp=divp;
    divp = pp1(:,2:cols+1)-pp1(:,1:cols)+pp2(2:rows+1,:)-pp2(1:rows,:);
    
    %sub-step2；update q^{k+1}
    %keep q
    Q=q;
    q = (uu -beta*(divp- ur)) / (1+beta);
    
    %sub-step3 update u^{k+1}
    Erru = beta*(ur-divp-q);
    uu = uu + Erru;
 
    % residual ||divp^k+1+q^k+1-f||_2
    residual = norm(divp+q-ur, 2);
    Aresidual(i) = residual;
    
    % ||uk+1-uk||
    erru =sqrt(sum(sum(Erru.*Erru)));
    
    
    X2(i)=(i);
    Y2(i)=log(erru);
    
    
    if erru < epsilon
        break;
    end
end

% keep tu
u2=uu;

% average error of uu and Ur
erruuUr=sum(sum(abs(uu-Ur)))/imgSz;
msg11 = sprintf('number of ADMM = %u. \n', i);
disp(msg11);
msg22 = sprintf('average error of uu and Ur = %u. \n', erruuUr);
disp(msg22);
msg33 = sprintf('||uk+1 - uk|| = %u. \n', erru);
disp(msg33);

figure('Position', [100, 100, 400, 300]);
plot(X2,Y2,'-b.','LineWidth', 0.1, 'MarkerSize', 1);
hold on;
plot(X1,Y1,'-r.','LineWidth', 0.1, 'MarkerSize', 1); %
xlabel('It.','FontName','Times');
ylabel('log(Error)','FontName','Times');
xlim([0 1800]);
xticks(0:500:1800);
% ylim([-20 12]);
% yticks(-20:10:12);
legend('ADMM','AADMM')

figure;
subplot(1,4,1);
imshow(u2,[]);
subplot(1,4,2);
imshow(ur,[]);
subplot(1,4,3);
imshow(Ur,[]);
subplot(1,4,4);
imshow(u1,[]);
saveas(gcf,'lake.eps','epsc')

figure;
iterations_AADMM = 1:length(AAresidual);
semilogy(iterations_AADMM, AAresidual, '-r.', 'LineWidth', 1.5);
xlabel('It.', 'FontName', 'Times');
ylabel('$\|\mathrm{div}\tilde{\mathbf{p}}^{k}+\tilde{q}^{k}-f\|$', 'FontName', 'Times', 'Interpreter', 'latex');
xlim([0 1200]);
xticks(0:500:1200);
% title('Residual of AADMM', 'FontName', 'Times');
grid on;

figure;
iterations_ADMM = 1:length(Aresidual);
semilogy(iterations_ADMM, Aresidual, '-b.', 'LineWidth', 1.5);
xlabel('It.', 'FontName', 'Times');
ylabel('$\Vert \mathrm{div}\mathbf{p}^{k+1}+q^{k+1}-f\Vert$', 'FontName', 'Times', 'Interpreter', 'latex');
xlim([0 1800]);
xticks(0:500:1800);
% title('Residual of ADMM', 'FontName', 'Times');
grid on;


