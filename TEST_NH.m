
% ��iΪ1ʱȡ����ֵ20����5��λ�÷�ֵΪ4��ʣ��14���㶼Ϊ��
NHcode = [1 1 1 1 1 -1 1 1 -1 -1 1 -1 1 -1 1 1 -1 -1 -1 1];  %NH������
NHcode2 = [NHcode NHcode];
NHcode2 = NHcode2';
NHcal = zeros(20,1);
for i = 1:20
   NHcal(i) = NHcode * NHcode2(i:i+19); 
end
plot(NHcal);