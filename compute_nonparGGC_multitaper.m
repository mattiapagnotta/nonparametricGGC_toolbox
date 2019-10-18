
function [freq, causality, COH, iCOH, power, S] = compute_nonparGGC_multitaper(X, fs, fRes, freq, flgCond, NW)
%==========================================================================
% Written by M.Dhamala (August 2006).
% Modified by M.F.Pagnotta (June 2018).
%--------------------------------------------------------------------------
% INPUT
% - X:          bivariate or multivariate signals in the form of 3D matrix [nTime, nTrial, nChannel]
% - fs:         data sampling rate in Hz
% - fRes:       a desired frequency resolution (e.g. 1) - fRes is an optional input,  default value is fs/size(X,1)
% - freq:       frequencies at which causality, coherence and power are computed (vector)
% - flgCond:    if 1 do conditional Granger causality, while if 0 do pairwise causality
% - NW:         time-bandwidth parameter to control the number of tapers (given by 2*NW-1)
%--------------------------------------------------------------------------
% OUTPUT
% - freq:       frequencies at which causality, coherence and power are computed (vector)
% - causality:  Granger-Geweke causality between all pairs of signals (channels) [nFreqs, nChannel, nChannel]  
%               e.g., causality(:,2,1) means causality from 2 to 1, and 1 to 2 is causality(:,1,2)
%               self-causality are set to zero, i.e. causality(:,k,k) = 0 for all k
% - COH:        coherence between all pairs of signals [nFreqs, nChannel, nChannel]
% - iCOH:       imaginary-coherence between all pairs of signals [nFreqs, nChannel, nChannel]
% - power:      1-sided auto-spectra [nFreqs, nChannel]
% - S:          spectral density matrix [nChannel, nChannel, nFreqs]
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
% [1] Pagnotta et al., NeuroImage, 2018.
% [2] Dhamala et al., NeuroImage, 2008.
% [3] Dhamala et al., Phys. Rev. Lett., 2008.
% [4] Thomson, Proc. IEEE, 1982.
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


%--- Initialization --------------
[Nt, Ntr, Nc] = size(X);                                                    % Nt: timepoints, Ntr: trials, Nc: channels 

if nargin < 6,                  NW      = 4;        end                     % by default: # of tapers is 7
if nargin < 5,                  flgCond = 0;        end                     % by default: pairwise GGC
if nargin < 4,                  freq = 1:fRes:fs/2; end                     % vector of considered frequencies
if nargin < 3 || fRes > fs/Nt,  fRes = fs/Nt;       end                     % a lower frequency resolution is achieved by zero-padding

if fRes ~= abs(freq(2)-freq(1))
    error('Incompatibility between frequency resolution (fRes) and selected frequencies (freq)!'); %check compatibility
end



%--- Spectral quantities ----------
[S, f] = sig2mTspect_nv(X,fs,NW,fRes);  % Not vectorized, less memory-demanding, signals to multitapered auto&cross-spectra 

freq1  = [0 freq];                      % add zero freq
S      = S(:, :, 1:length(freq1));      % reduce S to frequency range of interest
causality = nan(length(freq1),Nc,Nc);

if nargout == 3
    spectra    = permute(S,[3 1 2]);
    COH        = S2coh(spectra);
    COH(1,:,:) = [];                    % omitting zero-frequency
elseif nargout > 3
    spectra     = permute(S,[3 1 2]);
    [COH, iCOH] = S2coh(spectra);
    COH(1,:,:)  = [];                   % omitting zero-frequency
    iCOH(1,:,:) = [];                   % omitting zero-frequency
    if nargout > 4
        for ichan = 1:Nc
            power(:,ichan) = spectra(:,ichan,ichan);    % one-sided power
        end
        power(1,:) = [];                %omitting zero-frequency
    end
end



%--- GGC computation --------------
% Bivariate case: pairwise causality of 2-channels
if Nc < 3
    [H, Z]    = wilson_sf(S, fs);                                           % Wilson's algorithm for spectral matrix factorization
    causality = hz2causality(H,S,Z,fs);                                     % Geweke's formula for Granger causality
    
% Multivariate case: more than 2 channels
elseif Nc >= 3
    % PAIRWISE causality (factorize two-channel spectra at a time)
    if ~flgCond
        for ii = 1:Nc-1
            for jj = ii+1:Nc
                S2       = S([ii jj],[ii jj],:);
                [H2, Z2] = wilson_sf(S2, fs);                               % Wilson's algorithm for spectral matrix factorization
                cs       = hz2causality(H2,S2,Z2,fs);                       % Geweke's formula for Granger causality
                causality(:,ii,jj) = cs(:,1,2);
                causality(:,jj,ii) = cs(:,2,1);
            end
            causality(:,ii,ii) = 0;                                         % self-causality is set to zero
        end
        causality(:,Nc,Nc) = 0;                                             % self-causality of the last channel is set to zero too
        
    % CONDITIONAL causality
    elseif flgCond
        [H0, SIGMA0] = wilson_sf(S, fs);                                    % FULL model (Wilson's algorithm)
        for y = 1:Nc
            y_omit = [1:y-1 y+1:Nc];
            [tmpHred0, tmpSIGMAred0] = wilson_sf(S(y_omit,y_omit,:), fs);   % REDUCED model (Wilson's algorithm)
            Hred0     = zeros(Nc,Nc,length(freq1)); 
            SIGMAred0 = zeros(Nc,Nc);
            Hred0(y_omit,y_omit,:)   = tmpHred0;
            SIGMAred0(y_omit,y_omit) = tmpSIGMAred0;
            for x = 1:Nc
                if x~=y
                    z = 1:Nc;
                    z([x y]) = [];                                          % indices of other variables (to condition out)
                    %-- FULL model ----------------
                    xyz = [x y z];
                    SIGMA = SIGMA0(xyz,xyz);
                    H     = H0(xyz,xyz,:);
                    b1  = 1;
                    b2  = b1+1;
                    b3  = (b2+1):length(xyz);
                    b   = [numel(b1) numel(b2) numel(b3)];
                    % normalization matrices
                    tmp1 = -SIGMA(b2,b1)/SIGMA(b1,b1);
                    tmp2 = -SIGMA(b3,b1)/SIGMA(b1,b1);
                    tmp3 = -(SIGMA(b3,b2)+tmp2*SIGMA(b1,b2))/(SIGMA(b2,b2)+tmp1*SIGMA(b1,b2));
                    p1_1 = [eye(b(1))    zeros(b(1),b(2))    zeros(b(1),b(3));
                            tmp1         eye(b(2))           zeros(b(2),b(3));
                            tmp2         zeros(b(3),b(2))    eye(b(3))];
                    p1_2 = [eye(b(1))            zeros(b(1),b(2))    zeros(b(1),b(3));
                            zeros(b(2),b(1))     eye(b(2))           zeros(b(2),b(3));
                            zeros(b(3),b(1))     tmp3                eye(b(3))];
                    P1   = p1_2*p1_1;
                    %-- REDUCED model --------------
                    xz   = [x z];
                    SIGMAred = SIGMAred0(xz,xz);
                    Hred     = Hred0(xz,xz,:);
                    bx1  = 1;
                    bx2  = (bx1+1):length(xz);
                    bx   = [numel(bx1) numel(bx2)];
                    P2   = [eye(bx(1))                              zeros(bx(1),bx(2));
                            -SIGMAred(bx2,bx1)/SIGMAred(bx1,bx1)    eye(bx(2))];
                    %-- conditional GC -------------
                    for kk = 1:size(S,3)
                        HH = H(:,:,kk)/P1;
                        B  = P2/Hred(:,:,kk);
                        BB = [B(bx1,bx1)        zeros(b(1),b(2))	B(bx1,bx2);
                              zeros(b(2),b(1))  eye(b(2))           zeros(b(2),b(3));
                              B(bx2,bx1)        zeros(b(3),b(2))    B(bx2,bx2)];
                        FF = BB*HH;
                        numer = abs(det(SIGMAred(bx1,bx1)));
                        denom = abs(det(FF(b1,b1)*SIGMA(b1,b1)*conj(FF(b1,b1))));
                        causality(kk,y,x) = log(numer./denom);
                    end
                    
                elseif x==y
                    causality(:,y,x) = 0;                                   % self-causality is set to zero
                    
                end
            end
        end
        
    end
end

causality(1,:,:) = [];                                                      % omitting zero-frequency causality
if nargout > 4
    S(:,:,1) = [];                                                          % omitting zero-frequency spectra (auto and cross)
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%        FUNCTIONS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
function [S, f] = sig2mTspect_nv(X, fs, nw, fRes)
%---------------------
% This function computes auto- & cross- spectra by using multitapers.
% 
% INPUT
% - X:          multichannel data [Nt, Ntr, Nc]  
% - fs:         data sampling rate in Hz
% - nw:         parameter to control the number of tapers (= 2*nw-1), good nw are 1.5, 2, 3, 4, 5, 6, or 7
% - fRes:       a desired frequency resolution (e.g. 1) - fRes is an optional input,  default value is fs/size(X,1)
% OUTPUT
% - S:          matrix of auto- and cross-spectra [Nc, Nc, nFreqs]
%---------------------

[Nt,Ntr,Nc] = size(X);        % Nt: timepoints, Ntr: trials, Nc: channels 

if nargin<4 || fRes>fs/Nt
       npad = 0;
       fRes = fs/Nt;  
end 
if (nargin==4) && (fRes<=fs/Nt) 
    npad = round((fs/fRes-Nt)/2);               % These many zeros will be padded on each side of the data
end
f = fs*(0:fix((Nt+2*npad)/2))/(Nt+2*npad);      % upto Nyquist-f

[tapers, v] = dpss(Nt+2*npad, nw, 2*nw-1);      % Return the 2*nw-1 most band-limited discrete prolate spheroidal sequences

S = zeros(Nc,Nc,Nt+2*npad);

for itrial = 1:Ntr
    for ii = 1:Nc
        Xft(:,:,ii) = mtfft(squeeze(X(:,itrial,ii)),tapers,fs,npad);
    end
    for ii = 1:Nc
        for jj = 1:Nc
            s(ii,jj,:) = squeeze(mean(Xft(:,:,ii).*conj(Xft(:,:,jj)),2)); %averaging over tapers 
        end
    end
    S = S + s; 
end
S = S/Ntr;                              % Averaging over trials
S = 2*S(:,:,1:fix(end/2)+1)/fs;         % factor 2 for make one-sided spectra
S(:,:,1) = S(:,:,1)/2;                  % dc-power doesn't double for one-sided case




%==========================================================================
function xf  = mtfft(data, tapers, fs, npad)
%---------------------
% This function multiply the time series from each trial by the preselected
% number of orthogonal tapers. The products are then Fourier-transformed.
% (Used inside the function sig2mTspect_nv.m)
% 
% INPUT
% - data:       time series from each trial
% - tapers:     preselected orthogonal tapers
% - fs:         sampling frequency
% - npad:       number of zeros for padding
% OUTPUT
% - xf:         Fourier-transformed of the products between each trial time series and the tapers
%---------------------

x0   = zeros(npad,size(data,2));
data = cat(1,x0,data);
data = cat(1,data,x0);
data = data(:,ones(1,size(tapers,2)));
data = data.*tapers;
xf   = fft(data,[],1);




%==========================================================================
function [Coh, iCoh] = S2coh(S)
%---------------------
% This function computes coherence and imaginary-coherence.
% 
% INPUT
% - S:          matrix of auto- and cross-spectra [nFreqs, Nc, Nc]
% OUTPUT
% - Coh:        coherence [nFreqs, Nc, Nc]
% - iCoh:       imaginary-coherence [nFreqs, Nc, Nc]
%---------------------

Nc = size(S,2);
for ii = 1: Nc
   for jj = 1: Nc
       Coh(:,ii,jj)  = abs(S(:,ii,jj)) ./ sqrt(S(:,ii,ii).*S(:,jj,jj)) ;            % Coherence (Coh)
       if nargout==2
           iCoh(:,ii,jj) = imag( S(:,ii,jj) ./ sqrt(S(:,ii,ii).*S(:,jj,jj)) );      % imaginary-Coherence (iCoh)
       end
   end
end




%==========================================================================
function causality = hz2causality(H, S, Z, fs)
%---------------------
% This function computes Granger-Geweke causality in the bivariate case.
% 
% INPUT
% - H:          transfer function
% - S:          matrix of auto- and cross-spectra [Nc, Nc, nFreqs]
% - Z:          noise covariance matrix
% - fs:         sampling frequency
% OUTPUT
% - causality:  Granger-Geweke causality
%---------------------

Nc = size(H,2);

for ii = 1:Nc
    for jj = 1:Nc
        if ii~=jj
            zc    = Z(jj,jj) - Z(ii,jj)^2/Z(ii,ii);
            numer = abs(S(ii,ii,:));
            denom = abs(S(ii,ii,:)-zc*abs(H(ii,jj,:)).^2/fs);
            causality(jj,ii,:) = log(numer./denom);
        end
    end
    causality(ii,ii,:) = 0;                 % self-causality set to zero
end
causality = permute(causality,[3 1 2]);     % [freq x channel from x channel to]


