
%==========================================================================
% Written by M.F.Pagnotta (June 2018).
%==========================================================================
% Simulation from Stokes and Purdon (PNAS, 2017) to test conditional
% Granger causality in the case of nonparametric methods
%--------------------------------------------------------------------------
% 3-node network with stationary interactions.
% The computation is performed usign the multitaper-based method.
%==========================================================================
% If you use nonparametricGGC_toolbox for a paper or talk please include
% the following references:
% 
% M.F. Pagnotta, M. Dhamala, G. Plomp, Benchmarking nonparametric Granger
% causality: robustness against downsampling and influence of spectral 
% decomposition parameters, NeuroImage. 183 (2018) 478–494. 
% https://doi.org/10.1016/j.neuroimage.2018.07.046
% 
% M. Dhamala, G. Rangarajan, M. Ding, Analyzing information flow in brain 
% networks with nonparametric Granger causality, NeuroImage. 41 (2008) 
% 354–362. https://doi.org/10.1016/j.neuroimage.2008.02.020.
% 
%_____________________
% Useful references:
% [1] Pagnotta et al., Data in Brief (submitted).
% [2] Pagnotta et al., NeuroImage, 2018.
% [3] Dhamala et al., NeuroImage, 2008.
% [4] Stokes & Purdon, PNAS, 2017.
% 
%--------------------------------------------------------------------------
% Extra:
% The function suplabel.m for plotting purpose has been made available by
% B.Barrowes at the following link:
% https://ch.mathworks.com/matlabcentral/fileexchange/7772-suplabel
% 
%__________________________________________________________________________
% [License]
%
% This file is part of nonparametricGGC_toolbox.
% Copyright (©) 2018, Mattia F. Pagnotta.
% 
% nonparametricGGC_toolbox is free software: you can redistribute it and/or
% modify it under the terms of the GNU General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or (at 
% your option) any later version.
% 
% nonparametricGGC_toolbox is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with nonparametricGGC_toolbox.  If not, see 
% <https://www.gnu.org/licenses/>.
%__________________________________________________________________________
%==========================================================================


% --- Model by Stokes & Purdon (2017, PNAS)-----
plag    = 3;
nNodes  = 3;
Fs      = 120;      % sampling frequency
L       = 500;      % length of each realization
N0      = 500;      % points to discard
nTrials = 1000;     % number of trials
sigma   = 1;

% --- nonparametric GGC: settings -----
doconditional = 1;              % 1 for conditional-GGC, 0 for pairwise-GGC
EndFreq       = Fs/2;           % Nyquist frequency
fRes          = Fs/L;           % frequency resolution (Hz)
freq = fRes:fRes:EndFreq;




%% Simulated model
% --- Interactions and frequencies -----
dt = 1/Fs;
r1 = 0.90;     f1o = 40;     f1 = f1o*dt;     cos_1 = cos(2*pi*f1);
r2 = 0.70;     f2o = 10;     f2 = f2o*dt;     cos_2 = cos(2*pi*f2);
r3 = 0.80;     f3o = 50;     f3 = f3o*dt;     cos_3 = cos(2*pi*f3);

% --- AR matrix -----
AR  = zeros(nNodes, nNodes, plag);

AR(:,:,1) = [2*r1*cos_1,             0,             0;
                 -0.356,    2*r2*cos_2,             0;
                      0,       -0.3098,    2*r3*cos_3;    ];

AR(:,:,2) = [     -r1^2,             0,             0;
                 0.7136,         -r2^2,             0;
                      0,           0.5,         -r3^2;   ];

AR(:,:,3) = [         0,             0,             0;
                 -0.356,             0,             0;
                      0,       -0.3098,             0;   ];

% --- Surrogate time series -----
w = randn(nNodes, N0+L, nTrials)*sigma;
y = w;
for nTr = 1:nTrials
    for i = plag+1:N0+L
        for po = 1:plag
            y(:,i,nTr) = y(:,i,nTr) + AR(:,:,po)*y(:,i-po,nTr);
        end
    end
end
clear nTr i
data = permute(y(:,N0+1:end,:), [2 1 3]);       %[L, nNodes, nTrials]
w    = permute(w(:,N0+1:end,:), [2 1 3]);




%% Nonparametric GGC (all trials)
X    = permute(data, [1 3 2]);                  %[L, nTrials, nNodes]
NW   = 4;
%..........................................................................
% Compute nonparametric GGC using multitaper-based approach (Dhamala et al., 2008)
[f, causality, cohTmp, icohTmp, powTmp] = compute_nonparGGC_multitaper(X, Fs, fRes, freq, doconditional, NW);
GGC = permute(causality, [3 2 1]);
if exist('cohTmp', 'var'),  COH  = permute(cohTmp,  [2 3 1]);     end       %[nNodes, nNodes, nFreqs]
if exist('icohTmp', 'var'), iCOH = permute(icohTmp, [2 3 1]);     end       %[nNodes, nNodes, nFreqs]
if exist('powTmp', 'var'),  POW  = permute(powTmp,  [2 1]);       end       %[nNodes, nFreqs]




%% Nonparametric GC (separate estimate for each single trial)
flgIntervals = 'percent';                   % option for confidence intervals
%-------------------------------------
if strcmp(flgIntervals, 'CI')
    alpha   = 0.05;
    nboot   = 10000;
    bootfun = @(x)(median(x));
elseif strcmp(flgIntervals, 'percent')
    percVal = [5 95];
end

GGC_trials = zeros([nNodes, nNodes, length(freq), nTrials]);
for indTrial = 1:size(data,3)
    X_1trial  = permute(data(:,:,indTrial), [1 3 2]);	%[L, 1, nNodes]
    [f, causality_1trial] = compute_nonparGGC_multitaper(X_1trial, Fs, fRes, freq, doconditional, NW);
    GGC_trials(:,:,:,indTrial) = permute(causality_1trial, [3 2 1]);
end
GGC_median = median(GGC_trials, 4);
if strcmp(flgIntervals, 'CI')
    % 95 percent CIs
    [CI_GGC, bootstat] = bootci(nboot, {bootfun, permute(GGC_trials, [4 1 2 3])}, 'alpha',alpha);
    CI_GGC = permute(squeeze(CI_GGC), [2 3 4 1]);
elseif strcmp(flgIntervals, 'percent')
    % Percentiles: 5% and 95%
    GGC_percentiles = prctile(GGC_trials, percVal, 4);
end





%% PLOTS
flg_plotSigleTr = 1;        % flag: plot single-trial GGC estimates
% ---
flg_plotPOW     = 1;
widthLINE_big   = 1.6;
widthLINE_small = 1.2;
fontName       = 'Times';
fontSizeSupLab = 15;  %title
fontSizeLabel  = 12;  %axis-labels
Ymax = [  0  0.4  0.4;
          6    0  0.4;
        0.4  0.8    0;];
% ---
RED = [255 0 0]/255; ORANGE = [255 165 0]/255; YELLOW = [255 255 0]/255;
GREEN = [0 210 0]/255; LIGHTBLUE = [0 191 255]/255; BLUE = [0 0 255]/255;
GREY = [128 128 128]/255; BLACK = [0 0 0]/255;

colorPOWback = [200 200 255]/255;
colorAREA = [210 210 210]/255;

for tmp_ind = 1:nNodes
    Sources{tmp_ind} = mat2str(tmp_ind);
end


% --- COMPREHENSIVE PLOT ------
TICKS_x = [0,10,40,50];
figure
pos = 0;
for i = 1:nNodes
    for j = 1:nNodes
        pos = pos + 1;
        if i~=j
            subplot(nNodes, nNodes, pos)
            if flg_plotSigleTr
                if strcmp(flgIntervals, 'CI')
                    y1 = squeeze(CI_GGC(i,j,:,1));
                    y2 = squeeze(CI_GGC(i,j,:,2));
                elseif strcmp(flgIntervals, 'percent')
                    y1 = squeeze(GGC_percentiles(i,j,:,1));
                    y2 = squeeze(GGC_percentiles(i,j,:,2));
                end
                y = [y1, (y2-y1)]; 
                ha = area(freq, y); hold on;
                set(ha(1), 'FaceColor','none')      % makes the bottom area invisible
                set(ha,    'LineStyle','none')
                set(ha(2), 'FaceColor', colorAREA);
                plot(freq, squeeze(GGC_median(i,j,:)), 'Color',GREY, 'LineWidth',widthLINE_small);	hold on;
            end
            plot(freq, squeeze(GGC(i,j,:)), 'Color',BLUE, 'LineWidth',widthLINE_big);	hold on;
            ylim([0,Ymax(i,j)]);
        elseif i==j
            if flg_plotPOW
                if exist('POW', 'var')
                    subplot(nNodes, nNodes, pos)
                    plot(freq, squeeze(POW(i,:)), 'Color',BLACK, 'LineWidth',widthLINE_big);
                    set(gca,'Color',colorPOWback)
                    title(sprintf('PSD_%d', i));
                end
            end
        end
        
        ax = gca;
        ax.XLim  = [0, EndFreq];
        ax.XTick = TICKS_x;
        if i == nNodes, ax.XLabel.String = ['from ' Sources{j}];  end
        if j == 1,      ax.YLabel.String = ['to ' Sources{i}];    end
        ax.XLabel.FontSize  = fontSizeLabel;
        ax.YLabel.FontSize  = fontSizeLabel;
        ax.FontName         = fontName;
    end
end
clear i j
[ax_X,h1] = suplabel('frequency (Hz)');
ax_X.XLabel.FontSize = fontSizeSupLab;
ax_X.XLabel.FontWeight = 'bold';
ax_X.FontName        = fontName;
[ax_Y,h2] = suplabel('GGC','y');
ax_Y.YLabel.FontSize = fontSizeSupLab;
ax_Y.YLabel.FontWeight = 'bold';
ax_Y.FontName        = fontName;


