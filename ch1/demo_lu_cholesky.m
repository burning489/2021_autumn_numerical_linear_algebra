clear, clc, close all;

size_list = 10:10:500;
len = length(size_list);
flu_duration_list = zeros(1,len);
plu_duration_list = zeros(1,len);
cho_duration_list = zeros(1,len);

repeat = 10;

for i = 1:len
    n = size_list(i);
    A = rand(n,n);
    A = A + A' + diag(0.1*n*ones(1,n));
    b = rand(n,1);
    % full lu
    start = cputime;
    for k = 1:repeat
        func_lu_solver(A,b,"full");
    end
    flu_duration_list(i) = (cputime-start)/repeat;
    % partial lu
    start = cputime;
    for k = 1:repeat
        func_lu_solver(A,b,"partial");
    end
    plu_duration_list(i) = (cputime-start)/repeat;
    % cholesky
    start = cputime;
    for k = 1:repeat
        func_cholesky_solver(A,b);
    end
    cho_duration_list(i) = (cputime-start)/repeat;
end

figure();
plot(size_list,flu_duration_list,'-*','DisplayName','flu');
hold on;
plot(size_list,plu_duration_list,'-^','DisplayName','plu');
plot(size_list,cho_duration_list,'-o','DisplayName','cholesky');
grid on;
hold off;
legend;
xlabel('matrix size');
ylabel('duration (sec)');

saveas(gcf, "./efficiency.png");



