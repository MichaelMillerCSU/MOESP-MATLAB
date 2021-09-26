%% 系统辨识总结一些问题：
% 1. 肯定能满足 计算系统辨识的时候 有三个assumption 需要进行设定:
%     1) rank(X0) = n
%     2) rank(Uo|k-1) = mk
%     3) rank(Uo|k-1;X0) = mk + n;
%     
%     怎么取数据使得上面三个假设满足？ 不满足怎么办？
%     
% 2. 没关系 在引入LQ Decomposition的时候， 需要用两个lemma （W = [Uo|k-1; Yo|k-1]）
%     1.rank(W) = rank(Uo|k-1; Yo|k-1) = km + n ;
%     2.W的线性组合可以表示任意输入输出长度为k的序列的组合值 （是由第一个lemma可以证明得到）
%     
%     如果数据不满足这些特征怎么办？ 我下面这个代码就不满足但是输出的结果还是好像差不多？
%     
% 3. SVD分解的时候，奇异值为什么降低这么快？ 感觉用阈值去在奇异值中选择系统阶数总是出问题，最后面还得自己来选？

%% Code
clc;
clear;
% load('newsample.mat');
% load('randData.mat')
load('data.mat')
k = 10;
Ts = 0.0001;
% u = [[0;4*ones(100000 - 1 ,1)] [zeros(10000,1);5000*ones(90000,1)] ]';
% y = [CA.Data' - max(CA.Data); T.Data' - min(T.Data)];
u = [4*ones(90000,1)'; 5000*ones(90000,1)'];
u = [zeros(2,10000)  u];
y = [CA.Data' - max(CA.Data); T.Data' - min(T.Data)];
y = [zeros(2,10000)  y];

% u = [CA0';Q'];
% y = [CA.Data' - 1.2; T.Data' - 439.78];

% u = [[0; 4 * ones(10000 - 1, 1)] [zeros(1000,1);5000 * ones(4000,1); zeros(1000,1);2500 * ones(4000,1)]]';
% y = [CA.Data' - max(CA.Data); T.Data' - min(T.Data)];
tt = 0 : Ts : size(u,2) * Ts - Ts;

m = size(u, 1);
p = size(y, 1);
% N = size(u, 2) - k + 1; %k order estimated
N = 40000 - k + 1; %k order estimated

hankel_u = ones(k * m, N);
hankel_y = ones(k * p, N);

for i = 1 : k
    hankel_u(2 * i - 1 : 2 * i - 1 + m - 1, :) = u(:, i : N + i - 1);
    hankel_y(2 * i - 1 : 2 * i - 1 + p - 1, :) = y(:, i : N + i - 1);
end
    
W = [hankel_u; hankel_y];
[Q,R] = qr(W');

L = triu(R)';
L11 = L(1 : k * m, 1 : k * m);
L12 = L(1 : k * m, k * m + 1 : k * m + k * p);
L21 = L(k * m + 1 : k * m + k * p, 1 : k * m);
L22 = L(k * m + 1 : k * m + k * p, k * m + 1 : k * m + k * p);

[U,S,V] = svd(L22);

% n=find(cumsum(S)>0.8*sum(S),1);

n = 2; %这里暂时没有什么好的方法

Ok = U(1 : k * m, 1 : n) * sqrt(S(1 : n, 1 : n));

C = Ok(1 : p , :);

A = pinv(Ok(1 : p * (k - 1), : )) * Ok(p + 1 : p * (k), :);


% % B and D idk how to do ??
% 
% U2 = U(:,n + 1 : size(U',1));
% 
% Z = U2'*L21* pinv(L11);  % Eq. 3.53
% XX = []; RR = [];
% for j = 1:k
%     XX = [XX; Z(:,m*(j-1)+1:m*j)];
%     Okj = Ok(1:p*(k-j),:);
%     Rj = [zeros(p*(j-1),p) zeros(p*(j-1),n);
%     eye(p) zeros(p,n); zeros(p*(k-j),p) Okj];
%     RR = [RR; U2'*Rj];
% end
% DB = pinv(RR)*XX;  % Eq. 3.57
% D = DB(1:p,:);
% B = DB(p+1:size(DB,1),:);

U2 = U(1 : k * p, n + 1 : end );
U2T = U2';

M = U2T * L21 * pinv(L11);

Next_to_U2 = U2T(: , p * 1 + 1 : end) * Ok(1 : p * (k - 1), :);

for i = 2 : k - 1
    Next_to_U2 = [Next_to_U2 ; U2T(: , p * i + 1 : end) * Ok(1 : p * (k - i), :)];
end

Next_to_U2 = [Next_to_U2 ; zeros(size(U2(1 : p , :)'))];

Li_Lk = U2T(: , 1 : p);

for i = 1 : k - 1
    Li_Lk = [Li_Lk ; U2T(: , p * i + 1 : p * (i + 1))];
end
Left = [Li_Lk, Next_to_U2];

Right = M(: , 1 : m);

for i = 1 : k - 1
    Right = [Right ; M(: , m * i + 1 : m * (i + 1))];
end

DB = pinv(Left) * Right;
D = DB(1 : p, 1 : m);
B = DB(p + 1 : end , 1 : m);

close
u = [[4*ones(100000 ,1)] [0*ones(10000,1) ; 5000*ones(90000,1)]]';
[Y,X] = dlsim(A,B,C,D,u);

figure (1) 
subplot(121) 
plot(tt(1, 1:100000), Y(1:100000, 1),'.')
hold on 
plot(tt(1, 1:100000), y(1,1:100000),'.')
hold off
ylabel('CA-CAS')
title('step inputs validation of CA - CAS')
subplot(122) 
plot(tt(1, 1:100000), Y(1:100000, 2),'.')
hold on
plot(tt(1, 1:100000), y(2,1:100000),'.')
hold off
ylabel('T-TS')
title('step inputs validation of T - TS')
% The step is in the heat rate input (u2) starting at 1 h with a magnitude of 5000 kJ/h.


u = [[4*ones(100000 ,1)] [0*ones(10000-100,1);1500*ones(100,1); 0*ones(90000,1)]]';
[Y,X] = dlsim(A,B,C,D,u);
figure (2) 
subplot(121) 
plot(tt(1, 1:100000), Y(1:100000, 1),'.')
ylabel('CA-CAS')
title('impulse inputs validation of CA - CAS')
subplot(122) 
plot(tt(1, 1:100000), Y(1:100000, 2),'.')
ylabel('T-TS')
title('impulse inputs validation of T - TS')
% To numerically simulate the impulse, a rectangular pulse of magnitude 1500 kJ/h in the heat rate input was applied for 36 s.


u = [[4*ones(100000 ,1)] [30000 * sin(8.72 * tt')]]';
[Y,X] = dlsim(A,B,C,D,u);
figure (3) 
load ('sinwave.mat')
subplot(221) 
plot(Y(1:100000, 1),'.')
ylabel('CA-CAS')
title('sinusoidal inputs validation of CA - CAS')
subplot(222) 
plot(CA.Data-1.2,'.');
ylabel('CA-CAS')
title('CSTR模型输出')
subplot(223) 
plot(Y(1:100000, 2),'.')
title('sinusoidal inputs validation of T - TS')
ylabel('T-TS')
subplot(224) 
plot(T.Data-438,'.');
ylabel('T-TS')
title('CSTR模型输出')
% The amplitude of the heat rate input sinusoid is 30,000 kJ/h with a frequency of 8.72 rad/h

A_T = (A - eye(2)) / 0.0001;
B_T = B / 0.0001;
C_T = C;
D_T = D;
disp("线性化的A：");
disp(C_T * A_T);
disp("线性化的B：");
disp(C_T * B_T + D_T );


%% trash

% AR = [-34.5 -0.473;1430 18.1];
% BR = [5.24 -8.09e-6;-11.6 4.57e-3];
% 
% % X0 = sqrt(S(1 : n, 1 : n))
% 
% x = zeros(2,10000);
% dt = 0.001;
%  
% 
% for i = 1:10000
%     dx = dt*(AR*x(:,i)+BR*u(:,i));
%     x(:,i+1) = x(:,i)+dx;
% end
% figure
% plot(x(1,:));
% figure;
% plot(x(2,:));

% Ck = eye(size(Ok, 1)) - Ok * inv(Ok' * Ok) * Ok';
% 
% phik = [eye(size(C)) zeros(size(C, 1), size(Ok(1 : p * (k - 1), :), 2))
%         zeros(size(Ok(1 : p * (k - 1), :), 1), size(C, 2)) Ok(p + 1 : p * (k), :)];
% 
% pinv(Ck) * Ck * phik
    
    
% % the model which is identified in paper(continuous time ) :
% A = [-34.5 -0.473;1430 18.1];
% B = [5.24 -8.09e-6;-11.6 4.57e-3];
% 
% A_T = [eye(2) + 0.001 * A];
% B_T = 0.001 * B;
% 
% disp(A)
% 
% disp(A_T)
% 
% 
% %% A    0.9655   -0.0005
%     1.4300    1.0181
% 
% 
% 
% %B    0.0052   -0.0000
% %    -0.0116    0.0000

