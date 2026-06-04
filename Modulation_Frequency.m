%% Comparative Acoustic Analysis of a Car Engine Condition Using Audio Signal Processing Techniques
% Work developed by Muhammad Ali Saqib & Shahariar Ryehan
% Course: AUDIO PROCESSING AND SYNTHESIS

clear; clc; close all;

% ==================== USER INPUT ====================
% Specify full paths to your normal and faulty engine WAV files
wavList = {...
    'normal_engine.wav', ...   % change path if needed
    'faulty_engine.wav' ...    % change path if needed
};

% Configuration parameters
cfg.winShort = 1024;           % short window length for STFT
cfg.winLong  = 4096;           % long window length for STFT
cfg.overlap  = 0.90;           % overlap fraction (90%)
cfg.showCWT  = true;           % show continuous wavelet transform if available

% ==================== CHECK FILES EXIST ====================
for i = 1:length(wavList)
    if ~exist(wavList{i}, 'file')
        error('File not found: %s', wavList{i});
    end
end

% ==================== PROCESS EACH FILE ====================
featureTbl = table();

for k = 1:length(wavList)
    feat = analyseOneFile(wavList{k}, cfg);
    featureTbl = [featureTbl; struct2table(feat)];
end

% Display extracted features
disp(' ');
disp('=================== EXTRACTED FEATURES ===================');
disp(featureTbl);

% ==================== LOCAL FUNCTION ====================
function feat = analyseOneFile(wavFile, cfg)
    % Read audio file
    [x, fs] = audioread(wavFile);
    
    % Convert stereo to mono if needed
    if size(x, 2) > 1
        x = mean(x, 2);
    end
    
    % Remove DC offset (detrend)
    x = detrend(x, 'constant');
    
    % Duration
    dur = numel(x) / fs;
    fprintf('\n--- %s (%.1f s at %d Hz) ---\n', wavFile, dur, fs);
    
    % Define STFT function (uses built-in stft or fallback to spectrogram)
    mySTFT = @(sig, win, hop) stftOrSpect(sig, fs, win, hop);
    
    % Hop sizes based on overlap
    winS = cfg.winShort;
    hopS = round((1 - cfg.overlap) * winS);
    winL = cfg.winLong;
    hopL = round((1 - cfg.overlap) * winL);
    
    % Compute STFTs
    [Ss, Fs, Ts] = mySTFT(x, winS, hopS);
    [Sl, Fl, Tl] = mySTFT(x, winL, hopL);
    
    % Hilbert envelope
    env = abs(hilbert(x));
    [Se, Fm, Tm] = mySTFT(env, winS, hopS);
    
    % --- Plotting ---
    figure('Name', sprintf('STFT/Envelope – %s', wavFile), 'NumberTitle', 'off');
    tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % STFT short window
    nexttile;
    imagesc(Ts, Fs, mag2db(abs(Ss)));
    axis xy;
    title('STFT – Short (1024)');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar;
    caxis([-100 0]);
    
    % STFT long window
    nexttile;
    imagesc(Tl, Fl, mag2db(abs(Sl)));
    axis xy;
    title('STFT – Long (4096)');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar;
    caxis([-100 0]);
    
    % Envelope STFT
    nexttile;
    imagesc(Tm, Fm, mag2db(abs(Se)));
    axis xy;
    title('Envelope STFT (Modulation)');
    xlabel('Time (s)');
    ylabel('Modulation Frequency (Hz)');
    colorbar;
    caxis([-100 0]);
    
    % Spectral Kurtosis
    nexttile;
    [Sz, Fz] = spectrogram(env, hamming(winS, 'periodic'), ...
                           winS - hopS, winS, fs, 'yaxis');
    SK = kurtosis(abs(Sz), 0, 2);
    plot(SK, Fz, 'LineWidth', 1.5);
    set(gca, 'YDir', 'reverse');
    grid on;
    xlabel('Spectral Kurtosis');
    ylabel('Frequency (Hz)');
    title('DIY Spectral Kurtosis (impulsive fault indicator)');
    
    % Optional Continuous Wavelet Transform
    if cfg.showCWT && exist('cwt', 'file') == 2
        figure('Name', sprintf('CWT – %s', wavFile), 'NumberTitle', 'off');
        cwt(x, fs, 'Wavelet', 'amor');
        title(sprintf('CWT Scalogram – %s', wavFile));
    end
    
    % --- Feature extraction ---
    feat.File = string(wavFile);
    feat.Duration_s = dur;
    feat.RMS = rms(x);
    feat.CrestFactor = max(abs(x)) / feat.RMS;
    feat.Skewness = skewness(x);
    feat.Kurtosis = kurtosis(x);
    
    % Modulation frequency (peak in envelope STFT)
    [~, pkM] = max(abs(Se), [], 1);
    feat.ModFreq_Hz = median(Fm(pkM));
    
    % Carrier frequency (peak in long-window STFT)
    [~, pkC] = max(abs(Sl), [], 1);
    feat.CarrierFreq_Hz = median(Fl(pkC));
end

% ==================== HELPER FUNCTION ====================
function [S, F, T] = stftOrSpect(x, fs, win, hop)
    % Use built-in stft if available (MATLAB R2019b+), otherwise fallback to spectrogram
    if exist('stft', 'file') == 2
        [S, F, T] = stft(x, fs, 'Window', hamming(win, 'periodic'), ...
                         'OverlapLength', win - hop, 'FFTLength', win);
    else
        [S, F, T] = spectrogram(x, hamming(win, 'periodic'), win - hop, win, fs, 'yaxis');
    end
end
