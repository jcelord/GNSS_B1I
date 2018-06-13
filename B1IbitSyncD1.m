function [dataOut dataRes] = B1IbitSyncD1( dataIn )
%B1ISyncD1 对B1I信号经跟踪处理后的-1,1序列进行比特位同步
%   input : dataIn    输入-1,1序列，其码周期为1ms
%   output: dataOut   输出经位同步后的-1,1序列，其码周期变为20ms
%           dataRes   输入数据最后不足一组（20)的序列
%   zj     2018/6/13


%   对于现有的37颗发射B1I信号的卫星，编号1-5为GEO卫星对应D2导航电文
%   而6-37号卫星为MEO/IGSO卫星，对应D1导航电文
%   D1导航电文码速率为50 bps,并调制有1 kbps的二次编码

    if 1 == length(dataIn(:,1))                     % 若为行序列，则转置为列序列
        data0 = dataIn';
    else
        data0 = dataIn;
    end

    NHcode = [1 1 1 1 1 -1 1 1 -1 -1 1 -1 1 -1 1 1 -1 -1 -1 1];  %NH码序列
    NHflag = true;                                 % 位同步标志，当找到位偏移量时，将其置为flase
    NHtd = 15;                                  % NH解码判决阈值
    M = 9;                                      % M之N判决法（参照捕获判决）
    N = 7;
    NHshift = 1;
    while NHflag
        count = 0;
        NHcal = NHcode * data0(NHshift:NHshift + 19);    % 将从输入序列中移位选取长度为20的序列与NH序列进行相乘
    %     if 20 == abs(NHcal) 
        if abs(NHcal) > NHtd                        % 由于可能发生误码，相乘结果可能小于20
          count = count + 1;
          for i = 1:M-1
              NHcal = NHcode * data0(NHshift + i*20:NHshift + i*20 + 19); 
              if abs(NHcal) > NHtd                        % 由于可能发生误码，相乘结果可能小于20
                 count = count + 1;
              end
          end
          if count >= N                            % M次计算中超过阈值的次数不少于N,判决通过
              NHflag = false;
          end
        else
            NHshift = NHshift + 1;                  % 当前不是导航数据比特起始位置，码偏移加一
        end
        if NHshift > 19                             % 超过19则代表处理失败
            error('没有找到位同步起始位置');
        end
    end

    NHcount = floor((length(data0) - NHshift + 1)/20);
    dataOut  = zeros(NHcount,1);
    for i = 1:NHcount                              % 每20ms进行解码
        dataOut(i) = ceil((NHcode * data0(NHshift+(i-1)*20:NHshift+i*20-1))/20);
    end
    dataRes = data0(NHshift + NHcount * 20 : end);
end

