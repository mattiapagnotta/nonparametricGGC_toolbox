
%==========================================================================
% Written by M.F.Pagnotta (June 2018).
%==========================================================================
% Simulation framework to evaluate:
% + SNR imbalance:              "SIGNAL-TO-NOISE RATIO IMBALANCE"
% + Independent noise:          "ADDITIVE NOISE"
% + Mixed noise:                "ADDITIVE NOISE"
% 
% The simulated model is adapted from Bastos & Schoffelen (2016) and
% extended to the multivariate case.
% The effects are shown for both the regular definition of GGC and the
% tr-GGC (which is obtained using time reversal testing).
% 
% Options:
% - possibility to vary number of nodes
% - option to select the define the influences in the network
% - option to selectively add noise to specific nodes only ('SNR' analysis)
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
% [4] Bastos & Schoffelen, Front. Syst. Neurosci., 2016.
% [5] Vinck et al., NeuroImage, 2015.
% [6] Nolte et al., Phys. Rev. Lett., 2008.
% [7] Haufe & Ewald, Brain Topogr., 2016.
% 
%--------------------------------------------------------------------------
% Extra:
% The function distinguishable_colors.m for plotting purpose has been made
% available by T.Holy at the following link:
% https://ch.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
% 
% The function mkpinknoise.m is part of the Berlin Brain Connectivity
% Benchmark (BBCB) simulion framework [7].
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


% --- Model settings -----
nNodes       = 3;                   % Select number of simulated nodes
sim_interact = [ 2,1;  3,2; ];      % Influences: 1->2 and 2->3
% sim_interact = [ 2,1;  3,1; ];      % Influences: 1->2 and 1->3


% --- Analysis settings -----
flg_froNorm  = 0;               % Use Frobenius norm for the simulation of data
alphaN       = 0:0.1:0.9;       % Control proportion between signals and noise sources

flg_Analysis = 'MIX';           % Define the analysis
%--options:--
% + 'SNR'  :    SNR imbalance between nodes
% + 'IND'  :    effects of independent (uncorrelated) noise
% + 'MIX'  :    effects of mixed (correlated) noise
%------------


if strcmp(flg_Analysis, 'SNR')
    addN_nodes = [2];                   % List of nodes where extra noise will be added
end
if strcmp(flg_Analysis, 'MIX')
    % Using the approach from Nolte et al. (2008)
    num_noiseS = nNodes;                % Number of noise sources: S = 2 with bivariate model in (Vinck et al., 2015)
    K = randn(nNodes, num_noiseS);      % (dimension: nodes-by-noise_sources)
    flg_Ncolor = 'white';               % 'white', 'pink', or 'white&pink'
end



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
fprintf(' Analysis: %s \n', flg_Analysis );
disp('-------------------------')





%==========================================================================
%% Simulated MVAR-model
%--------------------------------------------------------------------------
AR  = zeros(nNodes, nNodes, plag);

% Start with the model with no-interactions:
AR(:,:,1) = +0.5 * eye(nNodes);
AR(:,:,2) = -0.8 * eye(nNodes);

% Simulate interactions imposed (from j to i):
if ~isempty(sim_interact)
    for k = 1:size(sim_interact)
        ki = sim_interact(k,1);
        kj = sim_interact(k,2);
        AR(ki,kj,1) = +0.2;
        AR(ki,kj,2) = -0.1;
        clear ki kj
    end
    clear k
end





%==========================================================================
%% GGC ESTIMATION
%--------------------------------------------------------------------------
% Compute GGC on regular and time-reversed time series for each level of alphaN
%--------------------------------------------------------------------------
disp('                        ')
disp('GGC estimation...       ')

% ----- Initialize variables ---------
freq = fRes:fRes:EndFreq;
GGC  = zeros(nNodes, nNodes, numel(freq), numel(alphaN));      GGC_rev  = GGC;
COH  = zeros(nNodes, nNodes, numel(freq), numel(alphaN));      COH_rev  = COH;
iCOH = zeros(nNodes, nNodes, numel(freq), numel(alphaN));      iCOH_rev = iCOH;
POW  = zeros(nNodes, numel(freq), numel(alphaN));              POW_rev  = POW;


% ----- ESTIMATION... ---------
for kk = 1:numel(alphaN)
    
    fprintf('\t alphaN = %1.1f \n', alphaN(kk));
    
    % --- Surrogate time series ----------------------
    w = StandDEV*randn(nNodes, N0+L, nTrials);
    y = w;
    for nTr = 1:nTrials
        for i = plag+1:N0+L
            for po = 1:plag
                y(:,i,nTr) = y(:,i,nTr) + AR(:,:,po)*y(:,i-po,nTr);
            end
        end
    end
    clear nTr i
    data = permute(y(:,N0+1:end,:), [2 1 3]);                               % [L, nNodes, nTrials]
    w    = permute(w(:,N0+1:end,:), [2 1 3]);
    
    
    % --- Common referencing ----------------------
    if flg_froNorm==0
        if strcmp(flg_Analysis, 'SNR')
            if all(addN_nodes<=nNodes)
                % Add noise to specific nodes
                for nn = addN_nodes
                    N_nn         = StandDEV*randn(size(data(:,nn,:)));
                    data(:,nn,:) = (1-alphaN(kk))*data(:,nn,:) + alphaN(kk)*N_nn;
                end
                clear nn
            else
                error('Attempt to add noise to nonexistent node')
            end
        elseif strcmp(flg_Analysis, 'IND')
            % Independent (uncorrelated) noise processes
            indN = StandDEV*randn(size(data));
            data = (1-alphaN(kk))*data + alphaN(kk)*indN;
        elseif strcmp(flg_Analysis, 'MIX')
            % Mixed (correlated) noise processes (Nolte et al., 2008)
            if strcmp(flg_Ncolor,'white')
                z = StandDEV*randn(num_noiseS, L*nTrials);                  % WHITE noise
            elseif strcmp(flg_Ncolor,'pink')
                z_temp = mkpinknoise(L*nTrials, num_noiseS)';               % PINK noise (function by G.Nolte)
                z = StandDEV*(1/std(z_temp(:))).*z_temp;                    % ... control for the variance (basically normalize)
            elseif strcmp(flg_Ncolor,'white&pink')
                z_temp = mkpinknoise(L*nTrials, num_noiseS)';
                z = StandDEV*(1/std(z_temp(:))).*z_temp + StandDEV*randn(num_noiseS, L*nTrials) ;
            end
            E = K*z;
            mixN = reshape(E, [nNodes, L, nTrials]);
            mixN = permute(mixN, [2 1 3]);                                  % [L, nNodes, nTrials]
            data = (1-alphaN(kk))*data + alphaN(kk)*mixN;
        end
        
    elseif flg_froNorm==1
        if strcmp(flg_Analysis, 'SNR')
            % SNR imbalance between the two nodes (add noise only to node 1)
            N1          = StandDEV*randn(size(data(:,1,:)));
            for nTr = 1:nTrials
                Fnorm_data = norm(data(:,1,nTr), 'fro');
                Fnorm_N1   = norm(N1(:,:,nTr),   'fro');
                data(:,1,nTr) = (1-alphaN(kk))*(data(:,1,nTr)./Fnorm_data) + alphaN(kk)*(N1(:,:,nTr)./Fnorm_N1);
            end
            clear nTr
        elseif strcmp(flg_Analysis, 'IND')
            % Independent (uncorrelated) noise processes
            indN = StandDEV*randn(size(data));
            for nTr = 1:nTrials
                Fnorm_data = norm(data(:,:,nTr),'fro');
                Fnorm_N    = norm(indN(:,:,nTr),   'fro');
                data(:,:,nTr) = (1-alphaN(kk))*(data(:,:,nTr)./Fnorm_data) + alphaN(kk)*(indN(:,:,nTr)./Fnorm_N);
            end
            clear nTr
        elseif strcmp(flg_Analysis, 'MIX')
            % Mixed (correlated) noise processes (Nolte et al., 2008)
            if strcmp(flg_Ncolor,'white')
                z = StandDEV*randn(num_noiseS, L*nTrials);                  % WHITE noise
            elseif strcmp(flg_Ncolor,'pink')
                z_temp = mkpinknoise(L*nTrials, num_noiseS)';               % PINK noise (function by G.Nolte)
%                 z = StandDEV*(1/std(z_temp(:))).*z_temp;                    % ... control for the variance (basically normalize)
                z = StandDEV*z_temp;
            elseif strcmp(flg_Ncolor,'white&pink')
                z_temp = mkpinknoise(L*nTrials, num_noiseS)';
                z = StandDEV*(1/std(z_temp(:))).*z_temp + StandDEV*randn(num_noiseS, L*nTrials) ;
            end
            E = K*z;
            mixN = reshape(E, [nNodes, L, nTrials]);
            mixN = permute(mixN, [2 1 3]);                                  % [L, nNodes, nTrials]
            for nTr = 1:nTrials
                Fnorm_data = norm(data(:,:,nTr),'fro');
                Fnorm_N    = norm(mixN(:,:,nTr),   'fro');
                data(:,:,nTr) = (1-alphaN(kk))*(data(:,:,nTr)./Fnorm_data) + alphaN(kk)*(mixN(:,:,nTr)./Fnorm_N);
            end
            clear nTr
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

for kk = 1:numel(alphaN)
    
    fprintf('\t alphaN = %1.1f \n', alphaN(kk));
    
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
widthLINE       = 0.7;
fontName        = 'Times';
fontSize        = 6;                    %ticks
fontSizeLabel   = 9;                    %axis-labels
fontSizeTit     = 10;                   %title
fontSizeLegend  = 8;                    %legend
TICKS_x         = [0:20:freq(end)];
%-------
SP_shift        = [ -0.08; -0.06; -0.04; ];         % shift of subplots depending on the column - horizontal shift
VT_shift        = [ +0.02; -0.01; -0.04; ];         % shift of subplots depending on the row    - vertical shift
limY_reg        = [0, 1.1*max(reg_MEASURE(:))];
limY_tr         = [0, 1.1*max(tr_MEASURE(:))];
%-------
COLORS          = distinguishable_colors(numel(alphaN));
colorPOWback    = [220 220 220]/255;
%-------
for tmp_ind = 1:nNodes
    Sources{tmp_ind} = mat2str(tmp_ind);
end
clear tmp_ind
for kk = 1:numel(alphaN)
    conditions{kk}   = mat2str(alphaN(kk));
end
clear kk
%-------


%=====================================
% --- GGC (forward time) -------------
%=====================================
pos = 0;
figure
for ii = 1:nNodes
    for jj = 1:nNodes
        pos = pos + 1;
        % Causal estimates
        if ii~=jj
            h   = subplot(nNodes, nNodes, pos);
            pos_SP = get(h, 'Position');
            pos_SP_new = pos_SP;
            pos_SP_new(1) = pos_SP_new(1) + SP_shift(jj);
            pos_SP_new(2) = pos_SP_new(2) + VT_shift(ii);
            set(h, 'Position', pos_SP_new); clear pos_SP pos_SP_new
            for kk = 1:numel(alphaN)
                plot(freq, squeeze(reg_MEASURE(ii,jj,:,kk)),  'Color',COLORS(kk,:), 'LineWidth',widthLINE);  hold on;
            end
            clear kk
            ax = gca;
            ax.XLim  = [0, EndFreq];
            ax.XTick = TICKS_x;
            ax.YLim  = limY_reg;
            ax.XLabel.String    = 'frequency (Hz)';
            ax.YLabel.String    = sprintf('GGC_%d_%d',ii,jj);
            ax.FontSize         = fontSize;
            ax.XLabel.FontSize  = fontSizeLabel;
            ax.YLabel.FontSize  = fontSizeLabel;
            ax.FontName         = fontName;
            if pos==nNodes-1
                legend(conditions, ...
                    'Fontsize',fontSizeLegend, ...
                    'FontWeight','bold', ...
                    'FontName',fontName, ...
                    'EdgeColor',[1 1 1]);
            end
        % PSD estimates
        elseif ii==jj
            if flg_plotPOW
                if exist('POW', 'var')
                    h = subplot(nNodes, nNodes, pos);
                    pos_SP = get(h, 'Position');
                    pos_SP_new = pos_SP;
                    pos_SP_new(1) = pos_SP_new(1) + SP_shift(jj);
                    pos_SP_new(2) = pos_SP_new(2) + VT_shift(ii);
                    set(h, 'Position', pos_SP_new); clear pos_SP pos_SP_new
                    for kk = 1:numel(alphaN)
                        plot(freq, squeeze(POW(ii,:,kk)), 'Color',COLORS(kk,:), 'LineWidth',widthLINE);       hold on
                    end
                    clear kk
                    set(gca,'Color',colorPOWback)
                    ax = gca;
                    ax.XLim  = [0, EndFreq];
                    ax.XTick = TICKS_x;
                    ax.XLabel.String    = 'frequency (Hz)';
                    ax.YLabel.String    = sprintf('PSD_%d',ii);
                    ax.FontSize         = fontSize;
                    ax.XLabel.FontSize  = fontSizeLabel;
                    ax.YLabel.FontSize  = fontSizeLabel;
                    ax.FontName         = fontName;
                end
            end
        end
    end
end
clear ii jj


%=====================================
% --- tr-GGC -------------------------
%=====================================
pos = 0;
figure
for ii = 1:nNodes
    for jj = 1:nNodes
        pos = pos + 1;
        % Causal estimates
        if ii~=jj
            h   = subplot(nNodes, nNodes, pos);
            pos_SP = get(h, 'Position');
            pos_SP_new = pos_SP;
            pos_SP_new(1) = pos_SP_new(1) + SP_shift(jj);
            pos_SP_new(2) = pos_SP_new(2) + VT_shift(ii);
            set(h, 'Position', pos_SP_new); clear pos_SP pos_SP_new
            for kk = 1:numel(alphaN)
                plot(freq, squeeze(tr_MEASURE(ii,jj,:,kk)),  'Color',COLORS(kk,:), 'LineWidth',widthLINE);  hold on;
            end
            clear kk
            ax = gca;
            ax.XLim  = [0, EndFreq];
            ax.XTick = TICKS_x;
            ax.YLim  = limY_tr;
            ax.XLabel.String    = 'frequency (Hz)';
            ax.YLabel.String    = sprintf('tr-GGC_%d_%d',ii,jj);
            ax.FontSize         = fontSize;
            ax.XLabel.FontSize  = fontSizeLabel;
            ax.YLabel.FontSize  = fontSizeLabel;
            ax.FontName         = fontName;
            if pos==nNodes-1
                legend(conditions, ...
                    'Fontsize',fontSizeLegend, ...
                    'FontWeight','bold', ...
                    'FontName',fontName, ...
                    'EdgeColor',[1 1 1]);
            end
        % PSD estimates
        elseif ii==jj
            if flg_plotPOW
                if exist('POW_rev', 'var')
                    h = subplot(nNodes, nNodes, pos);
                    pos_SP = get(h, 'Position');
                    pos_SP_new = pos_SP;
                    pos_SP_new(1) = pos_SP_new(1) + SP_shift(jj);
                    pos_SP_new(2) = pos_SP_new(2) + VT_shift(ii);
                    set(h, 'Position', pos_SP_new); clear pos_SP pos_SP_new
                    for kk = 1:numel(alphaN)
                        plot(freq, squeeze(POW_rev(ii,:,kk)), 'Color',COLORS(kk,:), 'LineWidth',widthLINE);       hold on
                    end
                    clear kk
                    set(gca,'Color',colorPOWback)
                    ax = gca;
                    ax.XLim  = [0, EndFreq];
                    ax.XTick = TICKS_x;
                    ax.XLabel.String    = 'frequency (Hz)';
                    ax.YLabel.String    = sprintf('PSD_%d',ii);
                    ax.FontSize         = fontSize;
                    ax.XLabel.FontSize  = fontSizeLabel;
                    ax.YLabel.FontSize  = fontSizeLabel;
                    ax.FontName         = fontName;
                end
            end
        end
    end
end
clear ii jj


