
function [freq, causality, COH, iCOH, power, S] = compute_nonparGGC_wavelet_Paul(X, fs, freq, pad0, flgCond, m)
%==========================================================================
% Written by M.Dhamala and M.F.Pagnotta (June 2018).
%--------------------------------------------------------------------------
% INPUT
% - X:          bivariate or multivariate signals in the form of 3D matrix [nTime, nTrial, nChannel]
% - fs:         data sampling rate in Hz
% - freq:       frequencies at which causality, coherence and power are computed (vector)
% - pad0:       0 or 1 (0 to do not pad with zeros, 1 to do padding with zeros)
% - flgCond:    if 1 do conditional Granger causality, while if 0 do pairwise causality
% - m:          wavelet parameter in Paul wavelet prototype
%--------------------------------------------------------------------------
% OUTPUT
% - freq:       frequencies at which causality, coherence and power are computed (vector)
% - causality:  Granger-Geweke causality between all pairs of signals (channels) [nFreqs, nTime, nChannel, nChannel]  
%               e.g., causality(:,:,2,1) means causality from 2 to 1, and 1 to 2 is causality(:,1,2)
%               self-causality are set to zero, i.e. causality(:,k,k) = 0 for all k
% - COH:        coherence between all pairs of signals [nFreqs, nTime, nChannel, nChannel]
% - iCOH:       imaginary-coherence between all pairs of signals [nFreqs, nTime, nChannel, nChannel]
% - power:      1-sided auto-spectra [nFreqs, nTime, nChannel]
% - S:          spectral density matrix [nChannel, nChannel, nFreqs, nTime]
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
% [4] Torrence & Compo, Bull. Am. Meteorol. Soc., 1998.
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

if nargin < 3,  freq = 1:fs/2;  end
if nargin < 4,  pad0 = 1;       end
if nargin < 5,  flgCond = 1;    end                                         % by default: conditional GGC
if nargin < 6,  m = 4;          end                                         % parameter in Paul wavelet



%--- Spectral quantities ----------
if nargout == 2
    S = xwt_cpaul_nv(X,fs,freq,pad0,m);   % S: 4-D spectral matrix
end

if nargout == 3
    [S, COH] = xwt_cpaul_nv(X,fs,freq,pad0,m);
elseif nargout > 3
    [S, COH, iCOH] = xwt_cpaul_nv(X,fs,freq,pad0,m);
    if nargout > 4
        for ichan = 1:Nc
            power(:,:,ichan) = S(ichan,ichan,:,:);
        end
    end
end

freq1 = [0 freq];                   % add zero freq
S     = extrapolate_for_zfreq(S);   % extrapolating for zero-frequency

causality = nan(length(freq1), Nt, Nc, Nc);



%--- GGC computation --------------
% Bivariate case: pairwise causality of 2-channels
if Nc < 3
    for itime = 1:1:size(S,4)
        s         = squeeze(S(:,:,:,itime)); 
        [H, Z]    = wilson_sf(s,fs);                                        % Wilson's algorithm for spectral matrix factorization
        cs = hz2causality(H,s,Z,fs);                                        % Geweke's formula for Granger causality
        causality(:,itime,:,:) = permute(cs, [1 4 2 3]);
    end
    
% Multivariate case: more than 2 channels
elseif Nc >= 3
    % PAIRWISE causality (factorize two-channel spectra at a time)
    if ~flgCond
        for itime = 1:1:size(S,4)
            for ii = 1:Nc-1
                for jj = ii+1:Nc
                    S2 = squeeze(S([ii jj],[ii jj],:,itime));
                    [H2, Z2] = wilson_sf(S2,fs);                            % Wilson's algorithm for spectral matrix factorization
                    cs = hz2causality(H2,S2,Z2,fs);                         % Geweke's formula for Granger causality
                    causality(:,itime,ii,jj) = cs(:,1,2);
                    causality(:,itime,jj,ii) = cs(:,2,1);
                end
                causality(:,itime,ii,ii) = 0;                               % self-causality is set to zero
            end 
            causality(:,itime,Nc,Nc) = 0;                                   % self-causality of the last channel is set to zero too
        end
        
    % CONDITIONAL causality
    elseif flgCond
        for itime = 1:1:size(S,4)
            Sinst = squeeze(S(:,:,:,itime));
            [H0, SIGMA0] = wilson_sf(Sinst, fs);                                    % FULL model (Wilson's algorithm)
            for y = 1:Nc
                y_omit = [1:y-1 y+1:Nc];
                [tmpHred0, tmpSIGMAred0] = wilson_sf(Sinst(y_omit,y_omit,:), fs);   % REDUCED model (Wilson's algorithm)
                Hred0     = zeros(Nc,Nc,length(freq1)); 
                SIGMAred0 = zeros(Nc,Nc);
                Hred0(y_omit,y_omit,:)   = tmpHred0;
                SIGMAred0(y_omit,y_omit) = tmpSIGMAred0;
                for x = 1:Nc
                    if x~=y
                        z = 1:Nc;
                        z([x y]) = [];      % indices of other variables (to condition out)
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
                        for kk = 1:size(Sinst,3)
                            HH = H(:,:,kk)/P1;
                            B  = P2/Hred(:,:,kk);
                            BB = [B(bx1,bx1)        zeros(b(1),b(2))	B(bx1,bx2);
                                  zeros(b(2),b(1))  eye(b(2))           zeros(b(2),b(3));
                                  B(bx2,bx1)        zeros(b(3),b(2))    B(bx2,bx2)];
                            FF = BB*HH;
                            numer = abs(det(SIGMAred(bx1,bx1)));
                            denom = abs(det(FF(b1,b1)*SIGMA(b1,b1)*conj(FF(b1,b1))));
                            causality(kk,itime,y,x) = log(numer./denom);
                        end
                        
                    elseif x==y
                        causality(:,itime,y,x) = 0;       % self-causality is set to zero
                        
                    end
                end
            end
            
        end %Time instant
        
    end
    
end

causality(1,:,:,:) = [];                                                    % omitting zero-frequency causality
if nargout > 4
    S(:,:,1,:) = [];                                                        % omitting zero-frequency spectra (auto and cross)
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%        FUNCTIONS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
function [S, Coh, iCoh]  = xwt_cpaul_nv(X, fs, freq, pad0, m)
%---------------------
% This function computes auto- & cross- spectra by using Paul wavelet
% transform.
% 
% INPUT
% - X:          bivariate or multivariate signals in the form of 3D matrix [nTime, nTrial, nChannel]
% - fs:         data sampling rate in Hz
% - freq:       frequencies at which causality, coherence and power are computed (vector)
% - pad0:       0 or 1 (0 to do not pad with zeros, 1 to do padding with zeros)
% - m:          parameter in Paul wavelet
% OUTPUT
% - S:          matrix of auto- and cross-spectra [nChannel, nChannel, nFreqs, nTimes]
% - Coh:        coherence between all pairs of signals [nFreqs, nTime, nChannel, nChannel]
% - iCoh:       imaginary-coherence between all pairs of signals [nFreqs, nTime, nChannel, nChannel]
%---------------------

[Nt, Ntr, Nc] = size(X);        % Nt: timepoints, Ntr: trials, Nc: channels 

S = zeros(Nc,Nc,length(freq),Nt);

for itrial = 1:Ntr
    for ichan = 1:Nc   %.........now taking the wavelet transform....
        Wx(:,:,ichan) = cwt_cpaul(X(:,itrial,ichan),fs,freq,pad0,m);
    end
    %........ computing spectral matrix ........
    for ii = 1:Nc
        for jj = 1:Nc
            s(ii,jj,:,:) = Wx(:,:,ii).*conj(Wx(:,:,jj));
        end
    end
    S = S + s;
end
S = S/Ntr;      % dividing by the number of trials
clear X Wx s

%....... computing coherence .......
if nargout>1
    for ii = 1:Nc
        for jj = 1:Nc
            Coh(ii,jj,:,:) = abs(S(ii,jj,:,:)) ./ sqrt(S(ii,ii,:,:).*S(jj,jj,:,:)) ;                % Coherence (Coh)
            if nargout>2
                iCoh(ii,jj,:,:) = imag(  S(ii,jj,:,:) ./ sqrt(S(ii,ii,:,:).*S(jj,jj,:,:))  );       % imaginary-Coherence (iCoh)
            end
        end
    end
    Coh  = permute(Coh,  [3 4 1 2]);
    if nargout>2, iCoh = permute(iCoh, [3 4 1 2]); end
end



%==========================================================================
function  wave = cwt_cpaul(signal, fs, freq, pad, m)
%---------------------
% This function computes the complex Paul wavelet transform.
% (Used inside the function xwt_cpaul_nv.m)
% 
% INPUT
% - signal:     1-D timeseries
% - fs:         data sampling rate in Hz
% - freq:       frequencies at which causality, coherence and power are computed (vector)
% - pad:        0 or 1 (0 to do not pad with zeros, 1 to do padding with zeros)
% - m:          parameter in Paul wavelet
% OUTPUT
% - wave:       complex Paul wavelet transform
%---------------------
% Reference: Torrence & Compo, BAMS (1998), pp.65.
%---------------------

n1      = length(signal);               % minimum scale: 2/fs
x(1:n1) = signal; %- mean(signal);

%-----zero-padding to speed up and to reduce edge-effects
if (pad == 1)
    x  = [x, zeros(1,1*n1)];
    n2 = length(x);
    x  = [x, zeros(1,2^nextpow2(n2)-n2)];   % nextpower of 2 zeros padded
end
n  = length(x);
xf = fft(x);        % Fourier transform of the (padded) time series

%....construct wavenumber array used in transform 
k  = [1:fix(n/2)];
k  = k.*((2.*pi)/(n/fs));
k  = [0., k, -k(fix((n-1)/2):-1:1)];
fourier_factor = (4*pi)/(2*m + 1);      % Fourier-Factor (Paul wavelet)

%....scale and wave arrays
scale  = 1./(fourier_factor*freq);
Nscale = length(scale);

wave = zeros(Nscale,n);     % define the wavelet array
wave = wave + 1i*wave;      % make it complex

% Loop through all scales and compute transforms
for iscale = 1:Nscale
    daughter = wavelet_kernel(m,k,scale(iscale));	
    wave(iscale,:) = ifft(xf.*daughter);
end
wave = wave(:,1:n1);        % wavelet transform without the zero-padding




%==========================================================================
function daughter = wavelet_kernel(m, k, scale)
%---------------------
% This function generates the daughter wavelets.
% (Used inside the function cwt_cmorl.m)
% 
% INPUT
% - m:          parameter in Paul wavelet
% - k:          vector of Fourier frequencies 
% - scale:      wavelet scale
% OUTPUT
% - daughter:   vector of the wave function
%---------------------
% Reference: Torrence & Compo, BAMS (1998), pp.65.
%---------------------

n = length(k);

expnt    = -(scale.*k).*(k > 0.);
norm1    = sqrt(scale*k(2))*(2^m / sqrt(m*factorial(2*m-1)))*sqrt(n);
norm2    = (scale.*k).^m;

daughter = norm1*norm2.*exp(expnt);
daughter = daughter.*(k > 0.);                      % Heaviside step function




%==========================================================================
function S = extrapolate_for_zfreq(S)
%---------------------
% This function extrapolates 
% 
% INPUT
% - S:          matrix of auto- and cross-spectra [nChannel, nChannel, nFreqs, nTimes]
% OUTPUT
% - S:          matrix of auto- and cross-spectra [nChannel, nChannel, nFreqs+1, nTimes]
%               where the new frequncy vector is: freq1 = [0 freq]
%---------------------

[N1,~,~,Nt] = size(S);

for ii = 1:Nt
    for k1 = 1:N1
        for k2 = 1:N1
            Y  = log(squeeze(S(k1,k2,1:5,ii)));
            yi = interp1(1:5,Y,0:5,'linear','extrap');
            y(k1,k2,1,ii) = real(exp(yi(1)));
        end
    end
end
S = cat(3,y,S);




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


