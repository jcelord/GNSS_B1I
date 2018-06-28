function simPara = initSimPara( )
%initSimPara 对仿真参数进行初始化
%   此处显示详细说明

    simPara.phi0_if = 0;        %[rad] Initial phase of RF signal;
    simPara.phi0_prn1 = 0;      %[rad] Initial phase of PRN signal;
    simPara.phi0_nh = 0;        %[rad] Initial phase of NH signal;
    simPara.phi0_bk = 0;        %[rad] Initial phase of BK signal;
    simPara.phi0_data = 0;      %[rad] Initial phase of data signal;
    
    simPara.doppler = 25;              %[Hz] Initial Doppler frequeny for RF-signal;
    simPara.codePhase = 3678;            %[sec] the delay in seconds rounded to a multiple of ts
    simPara.T = 3000e-3;                 %[sec] Signal length to be generated;
    simPara.PRN = 18;                    % PRN number [0;37];
    simPara.SNR = 0;                    % SNR  为零表示不加噪声，负数表示噪声信噪比

end

