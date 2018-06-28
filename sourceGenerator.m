function dataOut = sourceGenerator( settings,simPara )
%sourceGenerator 生成B1I仿真中频数据
%   input :     settings   系统参数，包括采样率，中频，等信息
%               simPara    仿真参数，包括PRN，多普勒，码相位等信息
%   output:     dataOut    仿真中频数据
    f_rf    = settings.f_rf;           %[Hz] BeiDou B1 nominal frequency;
    f_if    = settings.IF;           %[Hz] IF nominal frequency;
    f_prn   = settings.codeFreqBasis;          %[Hz] Nominal PRN-generator clock frequency;
    f_nh    = settings.f_nh;           %[Hz] Nominal Neiman-Huffman-generator clock frequency;
    f_data  = settings.f_data;         %[Hz] Nominal data rate of BeiDou D1 signal
    
    phi0_if   = simPara.phi0_if;      %[rad] Initial phase of RF signal;
    phi0_prn1 = simPara.phi0_prn1;      %[rad] Initial phase of PRN signal;
    phi0_nh   = simPara.phi0_nh;      %[rad] Initial phase of NH signal;
    phi0_bk   = simPara.phi0_bk;      %[rad] Initial phase of BK signal;
    phi0_data = simPara.phi0_data;      %[rad] Initial phase of data signal;
    
    fs = settings.samplingFreq;         %[Hz] Sampling frequency; 
    ts = 1/fs;                          %[sec]
    
    f_d = simPara.doppler;              %[Hz] Initial Doppler frequeny for RF-signal;
    code_delay = ts * simPara.codePhase;%[sec] the delay in seconds rounded to a multiple of ts
    T = simPara.T;                      %[sec] Signal length to be generated;
    prn_num = simPara.PRN;               % PRN number [0;37];
    noise = simPara.SNR;                 % SNR
    
    n_samples = floor(fs * T);
    prn_len  = settings.codeLength;        % [chips] PRN-code length in chips (bits);
    nh_len = 20;                           % Neumann-Hoffman code length in chips (bits);
    data_len = T * f_data;  %[bits] Navigation message length in bits;
    
    
    %% -----------Signal_generator--------------------------------------------------

t=(1:n_samples).*(1/fs);    % time samples for generating signal of T_elem length;

PRN1        = generateB1Icode(prn_num);                                % Generate PRN-code;
NH_original = [0 0 0 0 0  1 0 0  1  1 0  1 0  1 0 0  1  1  1  0];      % Generate Neiman-Huffman-code
NH          = [1 1 1 1 1 -1 1 1 -1 -1 1 -1 1 -1 1 1 -1 -1 -1  1];      % 0=>1, 1=>-1 
DATA        = (2 * round(rand(1, data_len))) - 1;                      % Temporary Stub for data bits after convolutional coder; (pseudorandom values of 1 and -1, of data_len values)

signal_I = [];

%Calculate phase for carrier;
phi_if = (phi0_if) + (2*pi*f_if*t) + (2*pi*f_d*t);                 % INITIAL PHASE + IF FREQ + DOPPLER FREQ
phi_if = mod(phi_if, (2*pi));                                      % Convert carrier phase to the range [0..2*pi];  

%Calculate phase for PRN1;
phi_prn1  = (phi0_prn1) + (f_prn*(t-code_delay));                  % (t-code_delay)
prn1_indx = (fix(mod(phi_prn1, prn_len)));					       

%Calculte phase for Neiman-Huffman;
phi_nh  = (phi0_nh) + (f_nh*(t-code_delay));                       % (t-code_delay) 
nh_indx = (fix( mod(phi_nh, nh_len) ));
 
%Calculate phase for data-message;
phi_data  = (phi0_data) + (f_data*(t-code_delay));                 % (t-code_delay)
data_indx = (fix(mod(phi_data, data_len)));

% rate_d = f_if ./ (f_if + f_d*t) ;                                     % DOPPLER Effect Factor
%   rate_d = (f_rf + f_d) /f_rf;
% df = 0;   % 多普勒变化率
% 
% 
% phi_if = (phi0_if) + (2*pi*f_if*t)+ 2*pi*f_d*t; % INITIAL PHASE + IF FREQ + DOPPLER FREQ (including Doppler frequency change)
% phi_if = mod(phi_if, (2*pi));                                       % Convert carrier phase to the range [0..2*pi];  
% 
% 
% %Calculate phase for PRN1;
% % phi_prn1  = (phi0_prn1) + (f_prn*(t-code_delay));                  % (t-code_delay)
% phi_prn1  = (phi0_prn1) + (f_prn*(t-code_delay).*rate_d);           % (t-code_delay) with DOPPLER Effect
% prn1_indx = (fix(mod(phi_prn1, prn_len)));                         % The function of "fix" is the same as "floor" 					       
% 
% %Calculte phase for Neiman-Huffman;
% % phi_nh  = (phi0_nh) + (f_nh*(t-code_delay));                       % (t-code_delay) 
% phi_nh  = (phi0_nh) + (f_nh*(t-code_delay).*rate_d);                % (t-code_delay) with DOPPLER Effect
% nh_indx = (fix( mod(phi_nh, nh_len) ));
%  
% %Calculate phase for data-message;
% % phi_data  = (phi0_data) + (f_data*(t-code_delay));                 % (t-code_delay)
% phi_data  = (phi0_data) + (f_data*(t-code_delay).*rate_d);          % (t-code_delay) with DOPPLER Effect
% data_indx = (fix(mod(phi_data, data_len)));


%Generate carrier;
carr_cos = cos(phi_if);      	% generate carrier (I)
carr_sin = sin(phi_if);         % generate carrier (Q)

%Generate PRN;
prn1 = PRN1(prn1_indx+1);    	% generate PRN for I-channel;

%Generate NH;
nh = NH(nh_indx+1);          	% generate NH-code;

%Generate DATA;
data = DATA(data_indx+1);  		% generate DATA;

s_I         = ((prn1 .* nh).* data);
%s_I         = ((prn1).* data);
signal_I = ((s_I ).*carr_cos);
signal_Q = ((s_I ).*carr_sin);
signal_test = signal_I + 1j * signal_Q;
%Let's add some noise:
%Noise generation is under development...
%noise = grand(1,length(signal_RSLT),"nor",0,noise_pwr);
if noise ~=0
    signal_test = awgn(signal_test,noise);
end
%signal_RSLT = signal_RSLT + noise;
%signal_RSLT = signal_RSLT / max(abs(signal_RSLT)); %Normalization after adding noise;

%Next step is to write signal to file;
%write_apend_complex_binary(signal_test, fileName, n_samples);
signal_save = zeros(1,2*n_samples);
signal_save(1:2:end-1) = real(signal_test);
signal_save(2:2:end)   = imag(signal_test);

fid = fopen(settings.fileName,'wb');
fwrite(fid,signal_save,settings.dataType);
fclose(fid);
dataOut = signal_test;
end

