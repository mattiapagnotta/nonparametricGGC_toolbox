
%==========================================================================
% Written by M.F.Pagnotta (June 2018).
%==========================================================================
% Simulation framework to evaluate:
% + "COMMON REFERENCE PROBLEM"
% 
% The effects are shown for both the regular definition of GGC and the
% tr-GGC (which is obtained using time reversal testing).
% The simulated model takes advantage of the relationship between AR(2)
% model coefficients and spectral density (Rodrigues and Andrade, 2015).
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
% [4] Trongnetrpunya et al., Front. Syst. Neurosci., 2016.
% [5] Bastos and Schoffelen, Front. Syst. Neurosci., 2016.
% [6] Rodrigues and Andrade, PeerJ, 2015.
% [7] Winkler et al., IEEE Trans. Signal Process., 2016.
% 
%--------------------------------------------------------------------------
% Extra:
% The function distinguishable_colors.m for plotting purpose has been made
% available by T.Holy at the following link:
% https://ch.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
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


% --- Common Reference settings -----
flg_froNorm = 1;                % Use Frobenius norm for the simulation of data
alphaCR     = 0:0.1:0.9;        % Control proportion between unipolar signals and common reference signal

flg_REFtype = 'FREQ';           % Define the type of common reference signal
%--options:--
% + 'WGN'  :    white Gaussian noise
% + 'SAME' :    using the signal from the third node (non-interacting), which has the same AR coefficients
% + 'FREQ' :    using the signal from the third node (non-interacting), which has a specific frequency component (f3o)
%------------

r1 = 0.8;      f1o = 40;
r2 = r1;       f2o = f1o;

if strcmp(flg_REFtype, 'SAME')
    r3  = r1;
    f3o = f1o;
elseif strcmp(flg_REFtype, 'FREQ')
    r3  = r1;
    f3o = 70;
%     f3o = 20;
end

% --- Model -----
flg_interact  = '1to2';         % Impose interaction
%--options:--
% + 'bidir' :   Bidirectional identical interactions
% + '1to2'  :   Unidirectional interaction 1->2
% + '2to1'  :   Unidirectional interaction 2->1
% + 'none'  :   NO interaction
% The third node is always non-interacting with the other two.
%------------


% --- Definition for time-reversal testing -----
flgTRGCdefinition = 'Diff';
%--options:--
% 'Conj':        conjustion-based tr-GGC
% 'Diff':        difference-based tr-GGC
% 'Net&Diff':    conjunction of net-GGC and difference-based tr-GGC
%------------


% Noise characteristics...
StandDEV = 1;

%--- MODEL ------
plag     = 2;
nNodes_0 = 3;
L        = 400;      % length of each realization
N0       = 1000;     % points to discard
nTrials  = 100;      % number of trials


% --- nonparametric GGC -----
doconditional = 1;              % 1 for conditional-GGC, 0 for pairwise-GGC
NW            = 4;              % select number of tapers: 2*NW-1
Fs            = 200;            % sampling frequency
EndFreq       = Fs/2;           % Nyquist frequency
fRes          = Fs/L;           % frequency resolution (Hz)


disp('-------------------------')
fprintf(' Option selected: %s \n', flg_REFtype );
disp('-------------------------')





%==========================================================================
%% Simulated MVAR-model
%--------------------------------------------------------------------------
% Variables initialization
AR  = zeros(nNodes_0, nNodes_0, plag);

dt = 1/Fs;
f1 = f1o*dt;    cos_1 = cos(2*pi*f1);
f2 = f2o*dt;    cos_2 = cos(2*pi*f2);
if or( strcmp(flg_REFtype,'SAME'), strcmp(flg_REFtype,'FREQ') )
f3 = f3o*dt;    cos_3 = cos(2*pi*f3);
end


% Univariate coefficients (diagonal)
AR(1,1,1) = 2*r1*cos_1;     AR(1,1,2) = -r1^2;
AR(2,2,1) = 2*r2*cos_2;     AR(2,2,2) = -r2^2;
if or( strcmp(flg_REFtype,'SAME'), strcmp(flg_REFtype,'FREQ') )
AR(3,3,1) = 2*r3*cos_3;     AR(3,3,2) = -r3^2;
end

% ---- Causal influences ---------
if strcmp(flg_interact,'bidir')
    % Bidirectional identical interactions
    AR(2,1,1) = -0.35;   AR(2,1,2) = +0.70;
    AR(1,2,1) = -0.35;   AR(1,2,2) = +0.70;
elseif strcmp(flg_interact,'1to2')
    % Unidirectional interaction 1->2
    AR(2,1,1) = -0.35;   AR(2,1,2) = +0.70;
elseif strcmp(flg_interact,'2to1')
    % Unidirectional interaction 2->1
    AR(1,2,1) = -0.35;   AR(1,2,2) = +0.70;
elseif strcmp(flg_interact,'none')
    % NO interaction
end





%==========================================================================
%% GGC ESTIMATION
%--------------------------------------------------------------------------
% Compute GGC on regular and time-reversed time series for each level of alphaCR
%--------------------------------------------------------------------------
disp('                        ')
disp('GGC estimation...       ')

% ----- Initialize variables ---------
freq = fRes:fRes:EndFreq;
GGC  = zeros(nNodes_0-1, nNodes_0-1, numel(freq), numel(alphaCR));      GGC_rev  = GGC;
COH  = zeros(nNodes_0-1, nNodes_0-1, numel(freq), numel(alphaCR));      COH_rev  = COH;
iCOH = zeros(nNodes_0-1, nNodes_0-1, numel(freq), numel(alphaCR));      iCOH_rev = iCOH;
POW  = zeros(nNodes_0-1, numel(freq), numel(alphaCR));                  POW_rev  = POW;


% ----- ESTIMATION... ---------
for kk = 1:numel(alphaCR)
    
    fprintf('\t alphaCR = %1.1f \n', alphaCR(kk));
    
    % --- Surrogate time series ----------------------
    w = StandDEV*randn(nNodes_0, N0+L, nTrials);
    y = w;
    for nTr = 1:nTrials
        for i = plag+1:N0+L
            for po = 1:plag
                y(:,i,nTr) = y(:,i,nTr) + AR(:,:,po)*y(:,i-po,nTr);
            end
        end
    end
    clear nTr i
    data = permute(y(:,N0+1:end,:), [2 1 3]);                               %[L, nNodes, nTrials]
    w    = permute(w(:,N0+1:end,:), [2 1 3]);
    
    
    % --- Common referencing ----------------------
    if flg_froNorm==0
        if strcmp(flg_REFtype, 'WGN')
            % Common reference as white Gaussian noise
            data(:,3,:) = [];                                               % Reduce data to the two nodes of interest (first step)
            nNodes = size(data,2);
            Cref = repmat( StandDEV*randn(L, 1, nTrials), [1 nNodes 1] );
            data = (1-alphaCR(kk))*data - alphaCR(kk)*Cref;
        elseif or( strcmp(flg_REFtype,'SAME'), strcmp(flg_REFtype,'FREQ') )
            % Consider channel 3 as common-reference
            data(:,1,:) = (1-alphaCR(kk))*data(:,1,:) - alphaCR(kk)*data(:,3,:);
            data(:,2,:) = (1-alphaCR(kk))*data(:,2,:) - alphaCR(kk)*data(:,3,:);
            data(:,3,:) = [];                                               % Reduce data to the two nodes of interest (last step)
            nNodes = size(data,2);
        end
        
    elseif flg_froNorm==1
        if strcmp(flg_REFtype, 'WGN')
            % Common reference as white Gaussian noise
            data(:,3,:) = [];                                               % Reduce data to the two nodes of interest (first step)
            nNodes = size(data,2);
            Cref = repmat( StandDEV*randn(L, 1, nTrials), [1 nNodes 1] );
            for nTr = 1:nTrials
                Fnorm_data = norm(data(:,:,nTr),'fro');
                Fnorm_Cref = norm(Cref(:,:,nTr),'fro');
                data(:,:,nTr) = (1-alphaCR(kk))*(data(:,:,nTr)./Fnorm_data) - alphaCR(kk)*(Cref(:,:,nTr)./Fnorm_Cref);
            end
            clear nTr
        elseif or( strcmp(flg_REFtype,'SAME'), strcmp(flg_REFtype,'FREQ') )
            % Consider channel 3 as common-reference
            for nTr = 1:nTrials
                Fnorm_data = norm(data(:,:,nTr),'fro');
                Fnorm_Cref = norm(data(:,3,nTr),'fro');
                data(:,1,nTr) = (1-alphaCR(kk))*(data(:,1,nTr)./Fnorm_data) - alphaCR(kk)*(data(:,3,nTr)./Fnorm_Cref);
                data(:,2,nTr) = (1-alphaCR(kk))*(data(:,2,nTr)./Fnorm_data) - alphaCR(kk)*(data(:,3,nTr)./Fnorm_Cref);
            end
            clear nTr
            data(:,3,:) = [];                                               % Reduce data to the two nodes of interest (last step)
            nNodes = size(data,2);
        end
        
    end
    
    
    % Time-reversed data
    data_rev = data(end:-1:1,:,:);
    
    
    % --- nonparametric GGC (Dhamala et al., 2008) ------
    % Traditional:
    X    = permute(data, [1 3 2]);                                          %[L, nTrials, nNodes]
    [~, causality, cohTmp, icohTmp, powTmp] = compute_nonparGGC_multitaper(X, Fs, fRes, freq, doconditional, NW);
    GGC(:,:,:,kk) = permute(causality, [3 2 1]);
    if exist('cohTmp', 'var'),  COH(:,:,:,kk)  = permute(cohTmp,  [2 3 1]);     end       %[nNodes, nNodes, nFreqs]
    if exist('icohTmp', 'var'), iCOH(:,:,:,kk) = permute(icohTmp, [2 3 1]);     end       %[nNodes, nNodes, nFreqs]
    if exist('powTmp', 'var'),  POW(:,:,kk)    = permute(powTmp,  [2 1]);       end       %[nNodes, nFreqs]
    
    % Time-reversed time series:
    X_rev    = permute(data_rev, [1 3 2]);                                  %[L, nTrials, nNodes]
    [~, causality_rev, cohTmp_rev, icohTmp_rev, powTmp_rev] = compute_nonparGGC_multitaper(X_rev, Fs, fRes, freq, doconditional, NW);
    GGC_rev(:,:,:,kk) = permute(causality_rev, [3 2 1]);
    if exist('cohTmp_rev', 'var'),  COH_rev(:,:,:,kk)  = permute(cohTmp_rev,  [2 3 1]);     end       %[nNodes, nNodes, nFreqs]
    if exist('icohTmp_rev', 'var'), iCOH_rev(:,:,:,kk) = permute(icohTmp_rev, [2 3 1]);     end       %[nNodes, nNodes, nFreqs]
    if exist('powTmp_rev', 'var'),  POW_rev(:,:,kk)    = permute(powTmp_rev,  [2 1]);       end       %[nNodes, nFreqs]
    
end
clear kk





%==========================================================================
%% TIME-REVERSAL TESTING
%--------------------------------------------------------------------------
% Apply time reversal testing according to the specific definition selected
% (flgTRGCdefinition)
%--------------------------------------------------------------------------
disp('                        ')
disp('Time reversal testing...')
TR_GGC   = zeros(size(GGC));

for kk = 1:numel(alphaCR)
    
    fprintf('\t alphaCR = %1.1f \n', alphaCR(kk));
    
    tmpGGC     = GGC(:,:,:,kk);
    tmpGGC_rev = GGC_rev(:,:,:,kk);
    
    % Compute NET Causality measures (for normal and time-reversed data)
    netGGC     = zeros(size(tmpGGC));
    netGGC_rev = zeros(size(tmpGGC_rev));
    for j = 1:1:nNodes
        for i = 1:1:nNodes
            netGGC(i,j,:)     = tmpGGC(i,j,:)     - tmpGGC(j,i,:);
            netGGC_rev(i,j,:) = tmpGGC_rev(i,j,:) - tmpGGC_rev(j,i,:);
        end
    end
    clear i j
    
    % Definitions of TRGC as in (Winkler et al., 2016):
    if strcmp(flgTRGCdefinition, 'Conj')
        Hlogic = and(  netGGC > 0 , netGGC_rev < 0  );
        tmpGGC = netGGC;                                % use netGGC as measure
        tmpGGC(~Hlogic) = 0;
    elseif strcmp(flgTRGCdefinition, 'Diff')
        Dmat   = netGGC - netGGC_rev;
        Hlogic = Dmat > 0;
        tmpGGC = Dmat;                                  % use Dmat as measure
        tmpGGC(~Hlogic) = 0;
    elseif strcmp(flgTRGCdefinition, 'Net&Diff')
        Dmat   = netGGC - netGGC_rev;
        Hlogic = and(  Dmat > 0 , netGGC > 0  );
        tmpGGC = Dmat;                                  % use Dmat as measure (or tmpGGC = GGC; % use forward-GC as measure)
        tmpGGC(~Hlogic) = 0;
    end
    
    TR_GGC(:,:,:,kk) = tmpGGC;      % Collect the time-reversed-GGC
    
    clear Hlogic tmpGGC tmpGGC_rev ...
        netGGC netGGC_rev
end
clear kk





%==========================================================================
%% PLOT RESULTS
%--------------------------------------------------------------------------
reg_MEASURE = GGC;              % traditional GGC
tr_MEASURE  = TR_GGC;           % time-reversed GGC

%-------
flg_plotPOW     = 1;
widthLINE       = 0.8;
fontName        = 'Times';
fontSize        = 6;                    %ticks
fontSizeLabel   = 9;                   %axis-labels
fontSizeLegend  = 8;                    %legend
TICKS_x         = [0:20:freq(end)];
%-------
SP_shift        = [ -0.08; -0.06; -0.04; ];         % shift of subplots (depending on the column)
VT_shift        = -0.03;                            % vertical shift for bottom subplots
limY_reg        = [0, 1.1*max(reg_MEASURE(:))];
limY_tr         = [0, 1.1*max(tr_MEASURE(:))];
%-------
COLORS          = distinguishable_colors(numel(alphaCR));
colorPOWback    = [220 220 220]/255;
%-------
for tmp_ind = 1:nNodes
    Sources{tmp_ind} = mat2str(tmp_ind);
end
clear tmp_ind
for kk = 1:numel(alphaCR)
    conditions{kk}   = mat2str(alphaCR(kk));
end
clear kk
%-------


% --- COMPREHENSIVE PLOT ------
dim1 = nNodes;
if flg_plotPOW, dim2 = nNodes+1;
else,           dim2 = nNodes;
end

figure

% -- GGC (from 1 to 2) ---
posCOL = 1;
h = subplot(dim1, dim2, posCOL);
pos_SP = get(h, 'Position');
pos_SP_new = pos_SP;
pos_SP_new(1) = pos_SP_new(1) + SP_shift(posCOL);
set(h, 'Position', pos_SP_new); clear pos_SP pos_SP_new
for kk = 1:numel(alphaCR)
    plot(freq, squeeze(reg_MEASURE(2,1,:,kk)),  'Color',COLORS(kk,:), 'LineWidth',widthLINE);  hold on;
end
clear kk
ax = gca;
ax.XLim  = [0, EndFreq];
ax.XTick = TICKS_x;
ax.YLim  = limY_reg;
ax.XLabel.String    = 'frequency (Hz)';
ax.YLabel.String    = 'GGC_2_1';
ax.FontSize         = fontSize;
ax.XLabel.FontSize  = fontSizeLabel;
ax.YLabel.FontSize  = fontSizeLabel;
ax.FontName         = fontName;


% -- tr-GGC (from 1 to 2) ---
posCOL = 2;
h = subplot(dim1, dim2, posCOL);
pos_SP = get(h, 'Position');
pos_SP_new = pos_SP;
pos_SP_new(1) = pos_SP_new(1) + SP_shift(posCOL);
set(h, 'Position', pos_SP_new); clear pos_SP pos_SP_new
for kk = 1:numel(alphaCR)
    plot(freq, squeeze(tr_MEASURE(2,1,:,kk)),  'Color',COLORS(kk,:), 'LineWidth',widthLINE);  hold on;
end
clear kk
ax = gca;
ax.XLim  = [0, EndFreq];
ax.XTick = TICKS_x;
ax.YLim  = limY_tr;
ax.XLabel.String    = 'frequency (Hz)';
ax.YLabel.String    = 'tr-GGC_2_1';
ax.FontSize         = fontSize;
ax.XLabel.FontSize  = fontSizeLabel;
ax.YLabel.FontSize  = fontSizeLabel;
ax.FontName         = fontName;



% -- GGC (from 2 to 1) ---
posCOL = 1;
h = subplot(dim1, dim2, posCOL+dim2);
pos_SP = get(h, 'Position');
pos_SP_new = pos_SP;
pos_SP_new(1) = pos_SP_new(1) + SP_shift(posCOL);
pos_SP_new(2) = pos_SP_new(2) + VT_shift;
set(h, 'Position', pos_SP_new); clear pos_SP pos_SP_new
for kk = 1:numel(alphaCR)
    plot(freq, squeeze(reg_MEASURE(1,2,:,kk)),  'Color',COLORS(kk,:), 'LineWidth',widthLINE);  hold on;
end
clear kk
ax = gca;
ax.XLim  = [0, EndFreq];
ax.XTick = TICKS_x;
ax.YLim  = limY_reg;
ax.XLabel.String    = 'frequency (Hz)';
ax.YLabel.String    = 'GGC_1_2';
ax.FontSize         = fontSize;
ax.XLabel.FontSize  = fontSizeLabel;
ax.YLabel.FontSize  = fontSizeLabel;
ax.FontName         = fontName;

% -- tr-GGC (from 2 to 1) ---
posCOL = 2;
h = subplot(dim1, dim2, posCOL+dim2);
pos_SP = get(h, 'Position');
pos_SP_new = pos_SP;
pos_SP_new(1) = pos_SP_new(1) + SP_shift(posCOL);
pos_SP_new(2) = pos_SP_new(2) + VT_shift;
set(h, 'Position', pos_SP_new); clear pos_SP pos_SP_new
for kk = 1:numel(alphaCR)
    plot(freq, squeeze(tr_MEASURE(1,2,:,kk)),  'Color',COLORS(kk,:), 'LineWidth',widthLINE);  hold on;
end
clear kk
ax = gca;
ax.XLim  = [0, EndFreq];
ax.XTick = TICKS_x;
ax.YLim  = limY_tr;
ax.XLabel.String    = 'frequency (Hz)';
ax.YLabel.String    = 'tr-GGC_1_2';
ax.FontSize         = fontSize;
ax.XLabel.FontSize  = fontSizeLabel;
ax.YLabel.FontSize  = fontSizeLabel;
ax.FontName         = fontName;
legend(conditions, ...
    'Fontsize',fontSizeLegend, ...
    'FontWeight','bold', ...
    'FontName',fontName, ...
    'EdgeColor',[1 1 1]);



% PSD
if flg_plotPOW
    if exist('POW', 'var')
        posCOL = 3;
        h = subplot(dim1, dim2, posCOL);
        pos_SP = get(h, 'Position');
        pos_SP_new = pos_SP;
        pos_SP_new(1) = pos_SP_new(1) + SP_shift(posCOL);
        set(h, 'Position', pos_SP_new); clear pos_SP pos_SP_new
        for kk = 1:numel(alphaCR)
            plot(freq, squeeze(POW(1,:,kk)), 'Color',COLORS(kk,:), 'LineWidth',widthLINE);       hold on
        end
        clear kk
        set(gca,'Color',colorPOWback)
        ax = gca;
        ax.XLim  = [0, EndFreq];
        ax.XTick = TICKS_x;
        ax.XLabel.String    = 'frequency (Hz)';
        ax.YLabel.String    = 'PSD_1';
        ax.FontSize         = fontSize;
        ax.XLabel.FontSize  = fontSizeLabel;
        ax.YLabel.FontSize  = fontSizeLabel;
        ax.FontName         = fontName;
        
        h = subplot(dim1, dim2, posCOL+dim2);
        pos_SP = get(h, 'Position');
        pos_SP_new = pos_SP;
        pos_SP_new(1) = pos_SP_new(1) + SP_shift(posCOL);
        pos_SP_new(2) = pos_SP_new(2) + VT_shift;
        set(h, 'Position', pos_SP_new); clear pos_SP pos_SP_new
        for kk = 1:numel(alphaCR)
            plot(freq, squeeze(POW(2,:,kk)), 'Color',COLORS(kk,:), 'LineWidth',widthLINE);       hold on
        end
        clear kk
        set(gca,'Color',colorPOWback)
        ax = gca;
        ax.XLim  = [0, EndFreq];
        ax.XTick = TICKS_x;
        ax.XLabel.String    = 'frequency (Hz)';
        ax.YLabel.String    = 'PSD_2';
        ax.FontSize         = fontSize;
        ax.XLabel.FontSize  = fontSizeLabel;
        ax.YLabel.FontSize  = fontSizeLabel;
        ax.FontName         = fontName;
    end
end


