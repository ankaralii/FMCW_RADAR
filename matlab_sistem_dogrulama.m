clear; clc;

%% FMCW Radar Parametreleri
c    = 3e8;          % ışık hızı (m/s)
f0   = 2.4e9;        % merkez frekans (Hz)
B    = 40e6;         % bant genişliği (Hz) — menzil çözünürlüğü = c/2B = 3.75 m
T    = 2e-3;         % chirp süresi (s)
fs   = 40e6;         % örnekleme hızı (Hz)
N    = T * fs;       % chirp başına örnek sayısı = 80000
lambda = c / f0;     % dalga boyu = 0.125 m

%% Türetilen parametreler
k    = B / T;        % frekans eğimi (Hz/s)
dR   = c / (2*B);   % menzil çözünürlüğü (m)
fprintf('Menzil çözünürlüğü: %.2f m\n', dR);
%% Chirp sinyali
t_vec = (0:N-1)' / fs;
chirp_tx = exp(1j * pi * k * t_vec.^2);
%% Hedef tanımla (birden fazla hedef eklenebilir)
hedefler = [
%  mesafe(m)  hız(m/s)   RCS(dBsm)
    10,        0,          10;
    20,        0,          10;
    35,        0,           5;
];

%% Alınan sinyali oluştur
rx_signal = zeros(N, 1);

for i = 1:size(hedefler, 1)
    R   = hedefler(i, 1);
    v   = hedefler(i, 2);
    rcs = hedefler(i, 3);
    
    % Gidiş-dönüş gecikmesi
    tau = 2 * R / c;
    
    % Gecikmiş chirp
    t_delayed = t_vec - tau;
    chirp_rx  = exp(1j * pi * k * t_delayed.^2);
    
    % Yol kaybı + RCS etkisi
    amplitude = 10^(rcs/20) / R^2;
    
    % Doppler fazı
    doppler = exp(1j * 2*pi * (2*v/lambda) * t_vec);
    
    rx_signal = rx_signal + amplitude * chirp_rx .* doppler;
end

% Gürültü ekle (SNR = 20 dB)
SNR_dB  = 20;
noise   = (10^(-SNR_dB/20)) * (randn(N,1) + 1j*randn(N,1)) / sqrt(2);
rx_signal = rx_signal + noise;
%% Dechirp
beat = rx_signal .* conj(chirp_tx);
beat = beat - mean(beat);  % DC gider

%% Range profile
NFFT = 4 * N;
win  = hann(N);
spec = abs(fft(beat .* win, NFFT));
spec = spec(1:NFFT/2);

f_ax = (0:NFFT/2-1) * fs / NFFT;
r_ax = f_ax * c * T / (2 * B);

figure;
plot(r_ax, 20*log10(spec + eps));
xlabel('Menzil (m)'); ylabel('Genlik (dB)');
title('FMCW Simülasyon - Range Profile');
xlim([0 60]); grid on;

% Hedef mesafelerini işaretle
for i = 1:size(hedefler,1)
    xline(hedefler(i,1), 'r--', sprintf('%.0f m', hedefler(i,1)));
end
%% Çoklu chirp — Range-Doppler
N_chirp = 64;        % CPI başına chirp sayısı
PRI     = T + 0.5e-3;  % chirp arası boşluk dahil tekrarlama süresi

rd_matrix = zeros(N, N_chirp);

for nc = 1:N_chirp
    rx_c = zeros(N, 1);
    t_c  = (0:N-1)' / fs + (nc-1) * PRI;  % mutlak zaman
    
    for i = 1:size(hedefler, 1)
        R   = hedefler(i,1);
        v   = hedefler(i,2);
        rcs = hedefler(i,3);
        tau = 2 * R / c;
        t_delayed = t_c - tau;
        chirp_rx  = exp(1j * pi * k * t_delayed.^2);
        amp = 10^(rcs/20) / R^2;
        dop = exp(1j * 2*pi * (2*v/lambda) * t_c);
        rx_c = rx_c + amp * chirp_rx .* dop;
    end
    
    noise = (10^(-SNR_dB/20)) * (randn(N,1)+1j*randn(N,1))/sqrt(2);
    rx_c  = rx_c + noise;
    
    beat_c = rx_c .* conj(exp(1j * pi * k * ((0:N-1)'/fs).^2));
    beat_c = beat_c - mean(beat_c);
    rd_matrix(:, nc) = beat_c;
end

%% Range-Doppler FFT
win2D   = hann(N) * hann(N_chirp)';
rd_fft  = fftshift(fft(fft(rd_matrix .* win2D, NFFT, 1), N_chirp, 2), 2);
rd_map  = 20*log10(abs(rd_fft(1:NFFT/2, :)) + eps);

%% Eksenler
v_max  = lambda / (4 * PRI);
v_ax   = linspace(-v_max, v_max, N_chirp);

figure;
imagesc(v_ax, r_ax, rd_map);
axis xy;
xlabel('Hız (m/s)'); ylabel('Menzil (m)');
title('Range-Doppler Haritası');
colorbar; clim([max(rd_map(:))-40, max(rd_map(:))]);
ylim([0 60]);