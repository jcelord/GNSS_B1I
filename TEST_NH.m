
% 在i为1时取最大幅值20，有5个位置幅值为4，剩余14个点都为零
NHcode = [1 1 1 1 1 -1 1 1 -1 -1 1 -1 1 -1 1 1 -1 -1 -1 1];  %NH码序列
NHcode2 = [NHcode NHcode];
NHcode2 = NHcode2';
NHcal = zeros(20,1);
for i = 1:20
   NHcal(i) = NHcode * NHcode2(i:i+19); 
end
plot(NHcal);