function [trackResults, channel]= tracking(fid, channel, settings)
% Performs code and carrier tracking for all channels.
%
% [trackResults, channel] = tracking(fid, channel, settings)
%
%   Inputs:
%       fid                 - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all satellites to be tracked 
%                              (prepared by preRum.m from acquisition results).
%       settings         - receiver settings.
%   Outputs:
%       trackResults  - tracking results (structure array).
%                              Contains in-phase prompt outputs and absolute spreading code's starting positions, 
%                              together with other observation data from the tracking loops.
%                              All are saved every millisecond.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis, Dennis M. Akos
% Some ideas by Dennis M. Akos
%--------------------------------------------------------------------------
%This program is free software; 
%you can redistribute it and/or modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation;
%either version 2 of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%See the GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License along with this program;
% if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,USA.
%--------------------------------------------------------------------------

%CVS record: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

% 2018/6/8  zj  
% add FLL , the carrier loop switch between FLL and PLL
% Firstly, the FLL works, when frequence error is less than 1Hz, the PLL
% works

%% ================== Initialize result structure =========================
% Channel status
trackResults.status = '-';                                      % No tracked signal, or lost lock

% The absolute sample in the record of the C/A code start
trackResults.absoluteSample = zeros(1, settings.msToProcess);

% Freq of the C/A code
trackResults.codeFreq = inf(1, settings.msToProcess);

% Frequency of the tracked carrier wave
trackResults.carrFreq = inf(1, settings.msToProcess);

% Outputs from the correlators (In-phase)
trackResults.I_P = zeros(1, settings.msToProcess);
trackResults.I_E = zeros(1, settings.msToProcess);
trackResults.I_L = zeros(1, settings.msToProcess);

% Outputs from the correlators (Quadrature-phase)
trackResults.Q_E = zeros(1, settings.msToProcess);
trackResults.Q_P = zeros(1, settings.msToProcess);
trackResults.Q_L  = zeros(1, settings.msToProcess);

% Loop discriminators
trackResults.dllDiscr = inf(1, settings.msToProcess);
trackResults.dllDiscrFilt  = inf(1, settings.msToProcess);
trackResults.pllDiscr = inf(1, settings.msToProcess);
trackResults.pllDiscrFilt = inf(1, settings.msToProcess);

% Copy initial settings for all channels
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% =============== Initialize tracking variables ===========================

codePeriods = settings.msToProcess;                   % for GPS one C/A code is one ms

% ----------------------------------------- DLL variables ----------------------------
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

% Summation interval
PDIcode = 0.001;

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, settings.dllDampingRatio, 1.0);

% ----------------------------------------- PLL variables ----------------------------
% Summation interval
PDIpll = 0.001;

% Calculate filter coefficient values
[tau1pll, tau2pll] = calcLoopCoef(settings.pllNoiseBandwidth, settings.pllDampingRatio, 0.25);

% ----------------------------------------- FLL variables ----------------------------
% Summation interval
PDIfll = 0.001;
IP_old = 1;
QP_old = 1;
% Calculate filter coefficient values
[tau1fll, tau2fll] = calcLoopCoef(settings.fllNoiseBandwidth, settings.fllDampingRatio, 0.25);

hwb = waitbar(0,'Tracking...');

%% ================= Start processing channels ==========================
for channelNr = 1:settings.numberOfChannels
    
    % Only process if PRN is non zero (acquisition was successful)
    if (channel(channelNr).PRN ~= 0)
        % Save additional information - each channel's tracked PRN
        trackResults(channelNr).PRN = channel(channelNr).PRN;
        
        switch_FLL_PLL = 0;   % the FLL_PLL switch, 0 -- FLL; 1 -- PLL
        count_FLL_PLL = 0;
        
        % Move the starting point of processing.
        % Can be used to start thesignal processing at any point in the data record (e.g. for long records).
        % In addition skip through that data file to start at the appropriate sample (corresponding to code phase).
        % Assumes sample type is schar (or 1 byte per sample) 
        %fseek(fid, settings.skipNumberOfBytes + channel(channelNr).codePhase-1, 'bof');
        fseek(fid, 2*(settings.skipNumberOfBytes + channel(channelNr).codePhase-1), 'bof');
        
        % Get a vector with the C/A code sampled 1x/chip
        b1iCode = codegen_B1I(channel(channelNr).PRN);
        
        % Then make it possible to do early and late versions
        b1iCode = [b1iCode(settings.codeLength) b1iCode b1iCode(1)];

        % -------------------------- Perform various initializations ------------------------------
        % define initial code frequency basis of NCO
        codeFreq = settings.codeFreqBasis;
        
        % define residual code phase (in chips)
        remCodePhase  = 0.0;
        
        % define carrier frequency which is used over whole tracking period
        carrFreq = channel(channelNr).acquiredFreq;
        carrFreqBasis = channel(channelNr).acquiredFreq;
        carrFreq_old = carrFreq;
        
        % define residual carrier phase
        remCarrPhase  = 0.0;

        % code tracking loop parameters
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;

        % carrier/Costas loop parameters
        oldCarrNCO_fll   = 0.0;
        oldCarrError_fll = 0.0;
        oldCarrError_pll = 0.0;
        carrError_pll = 0.0;
        oldCarrNco_pll = 0.0;
        carrNco_pll = 0.0;

        %=== Process the number of specified code periods =================
        for loopCnt =  1:codePeriods
            
%% ===================== GUI update =================================
            % The GUI is updated every 50ms.
            %  This way Matlab GUI is still responsive enough.
            % At the same time Matlab is not occupied all the time with GUI task.
            if (rem(loopCnt, 50) == 0)
                try
                    waitbar(loopCnt/codePeriods, hwb, ...
                                ['Tracking: Ch ', int2str(channelNr), ' of ', int2str(settings.numberOfChannels), '; PRN#', ...
                                int2str(channel(channelNr).PRN), '; Completed ',int2str(loopCnt), ' of ', int2str(codePeriods), ' msec']);           
                catch
                    % The progress bar was closed.
                    %  It is used as a signal to stop, "cancel" processing. Exit.
                    disp('Progress bar closed, exiting...');
                    return
                end
            end

%% ================== Read next block of data =========================           
            % Find the size of a "block" or code period in whole samples
            % Update the phasestep based on code freq (variable) and sampling frequency (fixed)
            codePhaseStep = codeFreq / settings.samplingFreq;
            blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);
            
            % Read in the appropriate number of samples to process this interation 
            [rawSignal, samplesRead] = fread(fid, blksize, settings.dataType);
%             [rawSignal_t, samplesRead] = fread(fid, 2 * blksize, settings.dataType);
%             rawSignal = rawSignal_t(1:2:end-1) + 1j * rawSignal_t(2:2:end);
            rawSignal = rawSignal';                            %transpose vector
            
            % If did not read in enough samples, then could be out of  data - better exit 
            if (samplesRead ~= blksize)
%             if (samplesRead ~=  2 * blksize)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
                % fclose(fid);
                return
            end

%% =========== Set up all the code phase tracking information =================== 
            % Define index into early code vector
            tcode = (remCodePhase-earlyLateSpc) : codePhaseStep : ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
            tcode2 = ceil(tcode) + 1;
            earlyCode = b1iCode(tcode2);
            
            % Define index into late code vector
            tcode = (remCodePhase+earlyLateSpc) : codePhaseStep : ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
            tcode2 = ceil(tcode) + 1;
            lateCode = b1iCode(tcode2);
            
            % Define index into prompt code vector
            tcode = remCodePhase : codePhaseStep : ((blksize-1)*codePhaseStep+remCodePhase);
            tcode2 = ceil(tcode) + 1;
            promptCode = b1iCode(tcode2);
            remCodePhase = (tcode(blksize) + codePhaseStep) - settings.codeLength;

%% ======Generate the carrier frequency to mix the signal to baseband ==================
            time = (0:blksize) ./ settings.samplingFreq;
            
            % Get the argument to sin/cos functions
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
            remCarrPhase = rem(trigarg(blksize+1), (2 * pi));
            
            % Finally compute the signal to mix the collected data to bandband
            carrCos = cos(trigarg(1:blksize));
            carrSin = sin(trigarg(1:blksize));

%% ============Generate the six standard accumulated values ======================
            % First mix to baseband
            iBasebandSignal = carrSin .* rawSignal;
            qBasebandSignal = carrCos .* rawSignal;

            % Now get early, late, and prompt values for each
            I_E = sum(earlyCode .* iBasebandSignal);
            Q_E = sum(earlyCode .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode .* iBasebandSignal);
            Q_L = sum(lateCode .* qBasebandSignal);
            
%% ============ Find PLL/FLL error and update carrier NCO ==========================
            if (switch_FLL_PLL == 0)
                % ---------------FLL---------------------
                real_Q = IP_old * Q_P - QP_old * I_P;
                real_I = IP_old * I_P + QP_old * Q_P;
                %carrError_fll = atan(real_Q / real_I) / (2.0 * pi);
                carrError_fll = atan2(real_Q,real_I) / (2.0*pi);
                IP_old = I_P;
                QP_old = Q_P;
            
                % Implement carrier loop filter and generate NCO command
                carrNco_fll = oldCarrNCO_fll + (tau2fll/tau1fll) * (carrError_fll - oldCarrError_fll) + carrError_fll * (PDIfll/tau1fll);
                oldCarrNCO_fll = carrNco_fll;
                adj_flag = abs(carrNco_fll - oldCarrNCO_fll);
                oldCarrError_fll = carrError_fll;
                                
                if adj_flag < 10    %看相邻两次跟踪到的多普勒频率之差是否小于1Hz
                    count_FLL_PLL = count_FLL_PLL + 1;
                else
                    count_FLL_PLL = 0;
                end

                if count_FLL_PLL >= 2    %看是否有连续两次跟踪到的多普勒频率之差小于1Hz，若有，则认为频率跟踪已很稳定而精确，可以转入PLL
                    switch_FLL_PLL = 1;
                end
            end
                % ______________PLL____________________
                % Implement carrier loop discriminator (phase detector)
                carrError_pll = atan(Q_P / I_P) / (2.0 * pi);

                % Implement carrier loop filter and generate NCO command
                carrNco_pll = oldCarrNco_pll + (tau2pll/tau1pll) * (carrError_pll - oldCarrError_pll) + carrError_pll * (PDIpll/tau1pll);
                oldCarrNco_pll = carrNco_pll;
                oldCarrError_pll = carrError_pll;
            
%             oldCarrNco = carrNco;
%             oldCarrError = carrError;
                       
            % Modify carrier freq based on NCO command
%             carrFreq = carrFreqBasis + carrNco;
            carrFreq = carrFreqBasis + carrNco_fll + carrNco_pll;
            carrFreq_old = carrFreq;
            
            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;
            trackResults(channelNr).carrDoppler(loopCnt) = carrFreq - settings.IF;

%% ============ Find DLL error and update code NCO ===========================
            codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
                                (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));
            
            % Implement code loop filter and generate NCO command
            codeNco = oldCodeNco + (tau2code/tau1code) * (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco = codeNco;
            oldCodeError = codeError;
            
            % Modify code freq based on NCO command
            codeFreq = settings.codeFreqBasis - codeNco;
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;

%% ============ Record various measures to show in postprocessing ==================
            % Record sample number (based on 8bit samples)
            trackResults(channelNr).absoluteSample(loopCnt) = ftell(fid);

            trackResults(channelNr).dllDiscr(loopCnt) = codeError;
            trackResults(channelNr).dllDiscrFilt(loopCnt) = codeNco;
            trackResults(channelNr).pllDiscr(loopCnt) = carrError_pll;
            trackResults(channelNr).pllDiscrFilt(loopCnt) = carrNco_pll;

            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;
        end % for loopCnt

        % If we got so far, this means that the tracking was successful
        % Now we only copy status, but it can be update by a lock detector if implemented
        trackResults(channelNr).status  = channel(channelNr).status;        
        
    end % if a PRN is assigned
    
end % for channelNr 

% Close the waitbar
close(hwb)
