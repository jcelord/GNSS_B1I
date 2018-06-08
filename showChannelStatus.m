function showChannelStatus(acqResults,channel, settings)
% Prints the status of all channels in a table.
%
% showChannelStatus(channel, settings)
%
%   Inputs:
%       channel     - data for each channel. It is used to initialize and at the processing of the signal (tracking part).
%       settings     - receiver settings

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

%CVS record: showChannelStatus.m,v 1.4.2.8 2006/08/14 11:38:22 dpl Exp $

fprintf('\n*=======*=====*===========*==========*=============*===========*======*\n');
fprintf(  '|   Channel   |   PRN   |      Threshold       |        Doppler      |         Frequency       |      Code Offset   |      Status  |\n');
fprintf(  '*=======*=====*===========*==========*=============*===========*======*\n');

for channelNr = 1 : settings.numberOfChannels
    if (channel(channelNr).status ~= '-')
        fprintf('|       %2d        |    %3d    |       %f      |         %5.0f         |      %2.5e     |      %10.3f     |      %1s     |\n', ...
                channelNr, ...
                channel(channelNr).PRN, ...
                acqResults.peakMetric(channel(channelNr).PRN), ...
                channel(channelNr).acquiredFreq - settings.IF, ...
                channel(channelNr).acquiredFreq, ...
                channel(channelNr).codePhase, ...
                channel(channelNr).status);
    else
        fprintf('|       %2d        |    ---    |     ------------    |     ------------     |     ----------------     |     -----------     |    Off   |\n', channelNr);
    end
end

fprintf('*=======*=====*===========*==========*=============*===========*======*\n\n');

