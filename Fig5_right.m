close all;clc;clear all;

Ur=double(imread('lake.png'))/255;
ur=double(imread('noise_lake.png'))/255;

epsilon = 1e-10;

[rows, cols]=size(ur);
alpha = 0.1;
steps=0.1;

t=1;
numIter = 5e3;
imgSz = rows*cols;

% AADMM beta=1
beta = 1;
q =zeros(rows,cols);
Q=q;
tq=q;
Tq=q;

uu = ur;
tu=uu;

pp1 = zeros(rows, cols+1);
pp2 = zeros(rows+1, cols);
divp = pp1(:,2:cols+1)-pp1(:,1:cols)+pp2(2:rows+1,:)-pp2(1:rows,:);

X1=[];
Y1=[];

for i = 1:numIter
    
    %sub-step1: update tp^(k+1)=p^{k+1}
    pt = divp +q- (ur+uu/beta);
    pp1(:,2:cols) = steps*(pt(:,2:cols) - pt(:,1:cols-1)) + pp1(:,2:cols);
    pp2(2:rows,:) = steps*(pt(2:rows,:) - pt(1:rows-1,:)) + pp2(2:rows,:);
    
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

    X1(i)=log(i);
    Y1(i)=log(Aerru);
 
    if Aerru < epsilon
        break;
    end
    
end

% average error of tu and Ur
errtuUr=sum(sum(abs(tu-Ur)))/imgSz;
msg1 = sprintf('number of AADMM (beta=1) = %u. \n', i);
disp(msg1);
msg2 = sprintf('average error of tu and Ur = %u. \n', errtuUr);
disp(msg2);
msg3 = sprintf('||uk+1-uk|| = %u. \n', Aerru);
disp(msg3);


% AADMM beta=0.5
beta = 0.5;
q =zeros(rows,cols);
Q=q;
tq=q;
Tq=q;

uu = ur;
tu=uu;

pp1 = zeros(rows, cols+1);
pp2 = zeros(rows+1, cols);
divp = pp1(:,2:cols+1)-pp1(:,1:cols)+pp2(2:rows+1,:)-pp2(1:rows,:);

X2=[];
Y2=[];

for i = 1:numIter
    
    %sub-step1: update tp^(k+1)=p^{k+1}
    pt = divp +q- (ur+uu/beta);
    pp1(:,2:cols) = steps*(pt(:,2:cols) - pt(:,1:cols-1)) + pp1(:,2:cols);
    pp2(2:rows,:) = steps*(pt(2:rows,:) - pt(1:rows-1,:)) + pp2(2:rows,:);
    
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

    X2(i)=log(i);
    Y2(i)=log(Aerru);
 
    if Aerru < epsilon
        break;
    end
    
end


% average error of tu and Ur
errtuUr=sum(sum(abs(tu-Ur)))/imgSz;

msg11 = sprintf('number of AADMM (beta=0.5) = %u. \n', i);
disp(msg11);
msg22 = sprintf('average error of tu and Ur = %u. \n', errtuUr);
disp(msg22);
msg33 = sprintf('||uk+1-uk|| = %u. \n', Aerru);
disp(msg33);


% AADMM beta=1.5
beta = 1.5;
q =zeros(rows,cols);
Q=q;
tq=q;
Tq=q;

uu = ur;
tu=uu;

pp1 = zeros(rows, cols+1);
pp2 = zeros(rows+1, cols);
divp = pp1(:,2:cols+1)-pp1(:,1:cols)+pp2(2:rows+1,:)-pp2(1:rows,:);

X3=[];
Y3=[];

for i = 1:numIter
    
    %sub-step1: update tp^(k+1)=p^{k+1}
    pt = divp +q- (ur+uu/beta);
    pp1(:,2:cols) = steps*(pt(:,2:cols) - pt(:,1:cols-1)) + pp1(:,2:cols);
    pp2(2:rows,:) = steps*(pt(2:rows,:) - pt(1:rows-1,:)) + pp2(2:rows,:);
    
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

    X3(i)=log(i);
    Y3(i)=log(Aerru);
 
    if Aerru < epsilon
        break;
    end
    
end


% average error of tu and Ur
errtuUr=sum(sum(abs(tu-Ur)))/imgSz;

msg111 = sprintf('number of AADMM (beta=1.5) = %u. \n', i);
disp(msg111);
msg222 = sprintf('average error of tu and Ur = %u. \n', errtuUr);
disp(msg222);
msg333 = sprintf('||uk+1-uk|| = %u. \n', Aerru);
disp(msg333);


figure
plot(X1,Y1,'-b.');
hold on;
plot(X2,Y2,'-r.'); 
hold on;
plot(X3,Y3,'-g.');
xlabel('log(It.)','FontName','Times');
ylabel('log(Error)','FontName','Times');
% xlim([0 500]);
% xticks(0:100:500);
% ylim([-20 12]);
% yticks(-20:10:12);
legend('\beta=1','\beta=0.5','\beta=1.5')

