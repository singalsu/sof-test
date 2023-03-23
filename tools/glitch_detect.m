function glitch_detect(test_tone_f)

if nargin < 1
	test_tone_f = 0;
end

fpath = '/home/singalsu/Downloads/';
fn = 'bat.wav.3Vbz0k'; % total mess
%fn = 'bat.wav.zhJhMf'; % glitch in right channel
%fn = 'bat.wav.Wodg37'; % low quality
%fn = 'bat.wav.WR5of5'; % Low signal level but OK
%fn = 'bat.wav.05IgGp'; % good
%fn = 'bat.wav.3UiPoZ'; % good

param.t_ignore_from_start = 30e-3;
param.t_ignore_from_end = 0;
param.t_step = 1.0e-3;
param.n_fft = 1024;
param.max_snr_drop = 6;
param.min_snr = 50;

% Read audio
[x, fs, frames, channels] = read_audio(fullfile(fpath, fn));

fprintf(1, 'Ch, Pass, Frequency, SNR min, SNR avg, Signal avg, Noise avg, Noise max, Num glitch, 1st glitch\n');

for ch = 1 : channels
	% Sanity check
	if ~signal_found(x(:, ch))
		continue;
	end

	% STFT
	[stft, f, t, param] = compute_stft(x(:, ch), fs, param);

	% Plot
	idstr = sprintf('%s ch%d', fn, ch);
	plot_specgram(stft, f, t, idstr);

	% Check frequency, get STFT frequency index
	[signal_idx, tone_f] = find_test_tone(stft, param, f, t, test_tone_f);

	% Get levels and SNR estimates
	meas = get_tone_levels(stft, param, f, t, signal_idx);

	% Plot SNR
	plot_levels(meas, t, idstr);

	% Checks for levels data
	meas = check_tone_levels(stft, param, meas, t);

	fprintf(1, '%d, %2d,     %5.0f,  %6.1f,  %6.1f,     %6.1f,    %6.1f,    %6.1f,     %6d,  %8.3f\n', ...
		ch, meas.success, tone_f, meas.min_snr_db, meas.mean_snr_db, ...
		meas.mean_signal_db, meas.mean_noise_db, ...
		meas.max_noise_db, meas.num_glitches, meas.t_glitch);

	if meas.num_glitches
		figure;
		t_s = (0:(frames - 1)) / fs;
		plot(t_s, x(:,ch));
		t_start = meas.t_glitch - param.t_step;
		t_end = meas.t_glitch + param.t_step;
		ax = axis();
		axis([t_start t_end ax(3:4)]);
		grid on
	end

	if ~meas.success
		meas = check_glitch_periodicity(x(:,ch), param, meas);
	end

end

end

%
% Helper functions
%

function meas = check_glitch_periodicity(x, param, meas)

[b, a] = butter(2, 0.95, 'high');
afx = abs(filter(b, a, x));
[~, locs]= findpeaks(afx, param.fs/1e3, 'MinPeakHeight', 0.20*max(afx));

figure
hist(diff(locs))
grid on
ylabel('Glitch ccurrence');
xlabel('Glitch intelval (ms)');

end

function success = signal_found(x)

% All zeros or DC
if abs(min(x) - max(x)) < eps
	success = 0;
else
	success = 1;
end

end

function plot_levels(meas, t, idstr)

figure;
subplot(3, 1, 1);
plot(t, meas.snr_db);
grid on;
ylabel('SNR (dB)');
title(idstr);
subplot(3, 1, 2);
plot(t, meas.signal_db);
grid on;
ylabel('Signal (dBFS)');
subplot(3, 1, 3);
plot(t, meas.noise_db);
grid on;
ylabel('Noise (dBFS)');
xlabel('Time (s)');

end

function meas = get_tone_levels(stft, param, f, t, signal_idx)

signal_i1 = signal_idx - param.win_spread;
signal_i2 = signal_idx + param.win_spread;
if signal_i1 < 1
	error('Too low tone frequency, increase FFT length');
end

signal_db = zeros(param.n_stft, 1);
noise_db = zeros(param.n_stft, 1);
snr_db = zeros(param.n_stft, 1);
for i = 1 : param.n_stft
	% Integrate signal power
	p_signal = sum(stft(signal_i1 : signal_i2, i));

	% Integrate noise power, but replace DC and signal with
	% average noise level.
	noise = stft(:, i);
	noise_avg = mean(noise(signal_i2 : end));
	noise(1 : param.win_spread) = noise_avg;
	noise(signal_i1 : signal_i2) = noise_avg;
	p_noise = sum(noise);

	% Sign, noise, and "SNR" as dB
	signal_db(i) = 10*log10(p_signal);
	noise_db(i) = 10*log10(p_noise);
	snr_db(i) = signal_db(i) - noise_db(i);
end

meas.noise_db = noise_db - param.win_gain;
meas.signal_db = signal_db - param.win_gain;
meas.snr_db = signal_db - noise_db;

i1 = find(t > param.t_ignore_from_start, 1, 'first');
i2 = find(t <= t(end) - param.t_ignore_from_end, 1, 'last');

meas.mean_signal_db = mean(meas.signal_db(i1 :i2));
meas.mean_noise_db = mean(meas.noise_db(i1 :i2));
meas.mean_snr_db = mean(meas.snr_db(i1 :i2));
meas.max_noise_db = max(meas.noise_db(i1 :i2));
meas.min_snr_db = min(meas.snr_db(i1 :i2));

end

function meas = check_tone_levels(stft, param, meas, t)

meas.t_glitch = 0;
meas.num_glitches = 0;
meas.success = true;

% Find glitches from SNR curve drops
i1 = find(t > param.t_ignore_from_start, 1, 'first');
i2 = find(t <= t(end) - param.t_ignore_from_end, 1, 'last');
idx = find(meas.snr_db(i1:i2) < meas.mean_snr_db - param.max_snr_drop);

if ~isempty(idx)
	idx = idx + i1 - 1;
	didx = diff(idx);
	meas.num_glitches = 1 + length(find(didx > 2));
	start_idx = idx(1);
	cidx = find(didx(1:end) > 1, 1, 'first');
	if isempty(cidx)
		end_idx = idx(end);
	else
		end_idx = idx(cidx);
	end
	meas.t_glitch = param.t_step * mean([start_idx end_idx] - 1) + ...
		0.5 * param.n_fft / param.fs;
	meas.success = false;
end

if meas.min_snr_db < param.min_snr
	meas.success = false;
end

end

function [signal_idx, tone_f] = find_test_tone(stft, param, f, t, test_tone_f)

if test_tone_f > 0
	err_ms = (f - test_tone_f) .^ 2;
	signal_idx = find(err_ms == min(err_ms));
	tone_f = f(signal_idx);
	return
end

i1 = find(t > param.t_ignore_from_start, 1, 'first');
i2 = find(t <= t(end) - param.t_ignore_from_end, 1, 'last');

signal_idx_all = zeros(i2 - i1 + 1, 1);
for i = i1 : i2
	signal_idx_all(i - i1 + 1) = find(stft(:, i) == max(stft(:, i)),1, 'first');
end

signal_idx = round(mean(signal_idx_all));
tone_f = f(signal_idx);

end

function [x, fs, frames, channels] = read_audio(fn)

if 1
	[x, fs] = audioread(fn);
else
	fs = 48e3;
	a = 1;
	%a = 10^(-3 / 20);
	x(:,1) = multitone(fs, 997, a, 1);
	%x(:,2) = multitone(fs, 821, a, 1);
	x = round(2^15 * x) / 2^15;
end

sx = size(x);
frames = sx(1);
channels = sx(2);

% idx = round(rand(1,1)*frames); x(idx) = rand(1,1) * 2 - 1;

end

function [stft, f, t, param] = compute_stft(x, fs, param)

sx = size(x);
if sx(2) > 1
	error('One channel only');
end

frames = sx(1);
win = kaiser(param.n_fft, 20);
param.win_spread = 7;
param.win_gain = -13.0379;
param.fs = fs;

param.n_step = fs * param.t_step;
param.n_stft = floor((frames - param.n_fft) / param.n_step);
n_half_fft = param.n_fft / 2 + 1;
scale = 1 / param.n_fft;
f = linspace(0, fs/2, n_half_fft);
t = (0 : (param.n_stft - 1)) * param.t_step;
stft = zeros(n_half_fft, param.n_stft);

for i = 1 : param.n_stft
	i1 = (i - 1) * param.n_step + 1;
	i2 = i1 + param.n_fft - 1;
	s1 = fft(x(i1 : i2) .* win) * scale;
	s2 = s1(1 : n_half_fft);
	s = s2 .* conj(s2);
	stft(:, i) = s;
end

end


function plot_specgram(stft, f, t, tstr)

figure;
clims = [-160 0];
imagesc(1e3 * t, f, 10*log10(stft), clims)
set(gca, 'ydir', 'normal');
colorbar;
grid on;
title(tstr);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');

end
