%bch编译码测试
clc
clear all
%% BCH（15,11,1）编码仿真
N = 15;
M = 11;
num_t = 100;        % 将长度M的原始码的段数
len = M * num_t;
coded = zeros(N*num_t,1);
coded_noise = zeros(N*num_t,1);
code_ori = ceil(2 * rand(len,1)) - 1;       %生成在0，1之间随机取值的原始矩阵

for num = 1:num_t
    code_tmp = zeros(N,1);
    regD = zeros(4,1);    % 每一组数据进来，寄存器清零
    for i = 1:N
       if i <= M
           code_tmp(i) = code_ori(i + (num - 1) * M);
           tmp0 = mod((code_ori(i + (num - 1) * M) + regD(4)),2);
           regD(4) = regD(3);
           regD(3) = regD(2);
           regD(2) = mod((regD(1) + tmp0),2);
           regD(1) = tmp0;
       else
           code_tmp(i) = regD(4);
           regD(4) = regD(3);
           regD(3) = regD(2);
           regD(2) = regD(1);
           regD(1) = 0;
       end
    end
    coded((num-1)*N+1:num*N) = code_tmp;  % 未加信道干扰前的码
    num_rand = ceil(N * rand);   % 随机选取一个位置进行电平翻转
    code_tmp(num_rand) = mod(code_tmp(num_rand) + 1,2);
    coded_noise((num-1)*N+1:num*N) = code_tmp;  % 加信道干扰后的码
end

%% BCH译码
code_dec = zeros(M*num_t,1);          %纠错解码后的序列
correctTable = [1,2,5,3,9,6,11,4,15,10,8,7,14,12,13];
correctTable = 16 - correctTable;
% correctTable = [12,13,1,14,5,2,7,15,11,6,4,3,10,8,9];
for num = 1:num_t
    code_tmp = coded_noise((num-1)*N+1:num*N);   % 逐段截取数据进行译码
    regD = zeros(4,1);             % 每一组数据进来，寄出去清零
    for i = 1:N
        tmp0 = regD(4);
        regD(4) = regD(3);
        regD(3) = regD(2);
        regD(2) = mod(regD(1)+tmp0,2);
        regD(1) = mod(code_tmp(i)+tmp0,2);
    end
    numReg = regD(4) * 8 + regD(3) * 4 + regD(2) * 2 + regD(1);
    if (numReg ~= 0)
        code_tmp(correctTable(numReg)) = mod(code_tmp(correctTable(numReg)) + 1,2);
    end
    code_dec((num-1)*M+1:num*M) = code_tmp(1:M);    %得到解码数据
end

code_compare = code_ori - code_dec;
s = sum(abs(code_compare))