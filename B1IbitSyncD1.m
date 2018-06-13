function [dataOut dataRes] = B1IbitSyncD1( dataIn )
%B1ISyncD1 ��B1I�źž����ٴ�����-1,1���н��б���λͬ��
%   input : dataIn    ����-1,1���У���������Ϊ1ms
%   output: dataOut   �����λͬ�����-1,1���У��������ڱ�Ϊ20ms
%           dataRes   �������������һ�飨20)������
%   zj     2018/6/13


%   �������е�37�ŷ���B1I�źŵ����ǣ����1-5ΪGEO���Ƕ�ӦD2��������
%   ��6-37������ΪMEO/IGSO���ǣ���ӦD1��������
%   D1��������������Ϊ50 bps,��������1 kbps�Ķ��α���

    if 1 == length(dataIn(:,1))                     % ��Ϊ�����У���ת��Ϊ������
        data0 = dataIn';
    else
        data0 = dataIn;
    end

    NHcode = [1 1 1 1 1 -1 1 1 -1 -1 1 -1 1 -1 1 1 -1 -1 -1 1];  %NH������
    NHflag = true;                                 % λͬ����־�����ҵ�λƫ����ʱ��������Ϊflase
    NHtd = 15;                                  % NH�����о���ֵ
    M = 9;                                      % M֮N�о��������ղ����о���
    N = 7;
    NHshift = 1;
    while NHflag
        count = 0;
        NHcal = NHcode * data0(NHshift:NHshift + 19);    % ����������������λѡȡ����Ϊ20��������NH���н������
    %     if 20 == abs(NHcal) 
        if abs(NHcal) > NHtd                        % ���ڿ��ܷ������룬��˽������С��20
          count = count + 1;
          for i = 1:M-1
              NHcal = NHcode * data0(NHshift + i*20:NHshift + i*20 + 19); 
              if abs(NHcal) > NHtd                        % ���ڿ��ܷ������룬��˽������С��20
                 count = count + 1;
              end
          end
          if count >= N                            % M�μ����г�����ֵ�Ĵ���������N,�о�ͨ��
              NHflag = false;
          end
        else
            NHshift = NHshift + 1;                  % ��ǰ���ǵ������ݱ�����ʼλ�ã���ƫ�Ƽ�һ
        end
        if NHshift > 19                             % ����19�������ʧ��
            error('û���ҵ�λͬ����ʼλ��');
        end
    end

    NHcount = floor((length(data0) - NHshift + 1)/20);
    dataOut  = zeros(NHcount,1);
    for i = 1:NHcount                              % ÿ20ms���н���
        dataOut(i) = ceil((NHcode * data0(NHshift+(i-1)*20:NHshift+i*20-1))/20);
    end
    dataRes = data0(NHshift + NHcount * 20 : end);
end

