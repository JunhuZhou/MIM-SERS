%% Raman Spectra Analysis Pipeline: Baseline Correction, Classification, and Visualization
% This script processes Raman spectra from an Excel file using the following steps:
% 1. Despiking
% 2. Savitzky-Golay smoothing
% 3. Baseline correction (airPLS)
% 4. Spectrum classification based on peak detection
% 5. Outlier filtering and statistical analysis
% 6. Visualization and heatmap plotting

clear; clc;

%% 1. Load data
[filename, filepath] = uigetfile('*.xlsx', 'Select the Raman Spectra File');
data = readmatrix(fullfile(filepath, filename));
raman_shift = data(:, 1);
intensities = data(:, 2:101);

%% 2. Remove overexposed spectra
OVEREXPOSED_THRESHOLD = 64535;
over_mask1 = max(intensities, [], 1) >= OVEREXPOSED_THRESHOLD;
intensities = intensities(:, ~over_mask1);
fprintf('Removed %d overexposed spectra (pre-correction).\n', sum(over_mask1));

%% 3. Set parameters
lambda = 1e7; order = 3; wep = 0.1; p = 0.05; itermax = 20;
sgolay_window = 11; sgolay_polynomial = 3;
despike_threshold = 5;
peak_region = 337:350;
ref_region = 60:73;

N = size(intensities, 2);
corrected_spectra = zeros(size(intensities));
despiked_spectra = zeros(size(intensities));
smoothed_spectra = zeros(size(intensities));
classification_result = zeros(1, N);
totalnumber = 0;

%% 4. Preprocess, correct, classify
for i = 1:N
    original = intensities(:, i);
    
    % Despiking
    medf = medfilt1(original, 3);
    spikes = abs(original - medf) > despike_threshold;
    original(spikes) = medf(spikes);
    despiked_spectra(:, i) = original;
    
    % Smoothing
    smoothed = sgolayfilt(original, sgolay_polynomial, sgolay_window);
    smoothed_spectra(:, i) = smoothed;
    
    % Baseline correction
    corr = airPLS(smoothed', lambda, order, wep, p, itermax)';
    corrected_spectra(:, i) = corr;
    
    % Classification
    Is = max(corr(peak_region));
    x = mean(corr(ref_region));
    sigma = std(corr(ref_region));
    if Is - x > 10 * sigma
        classification_result(i) = 1;
        totalnumber = totalnumber + 1;
    end
end

%% 5. Save interim data
base = fullfile(filepath, erase(filename, ".xlsx"));
writematrix([raman_shift despiked_spectra], [base '_Despiked_Spectra.xlsx']);
writematrix([raman_shift smoothed_spectra], [base '_Smoothed_Spectra.xlsx']);
writematrix([raman_shift corrected_spectra], [base '_Corrected_Spectra.xlsx']);
writematrix(classification_result', [base '_Classification.xlsx']);
fprintf('Saved despiked, smoothed, corrected, and classification data.\n');

%% 6. Heatmap of classification
classification_matrix = reshape(classification_result, [10, 10]);
figure('Color','w','Position',[100,100,800,600]);
h = heatmap(classification_matrix);
h.Title = 'Classification of Spectra ("0" or "1")';
h.XLabel = 'Spectrum Index (X)';
h.YLabel = 'Spectrum Index (Y)';
h.Colormap = [0.9 0.9 0.9; 0.2 0.6 1];
h.CellLabelColor = 'none';
h.FontSize = 18;
h.GridVisible = 'on';
h.ColorLimits = [0 1];

%% 7. Original & corrected spectral plots
figure('Color','w','Position',[100,100,800,900]);
subplot(3,1,1);
plot(raman_shift, intensities);
title('Original Raman Spectra');
xlabel('Raman Shift (cm^{-1})');
ylabel('Intensity');
box on; grid on;

subplot(3,1,2);
plot(raman_shift, corrected_spectra);
title('Baseline Corrected Spectra');
xlabel('Raman Shift (cm^{-1})');
ylabel('Intensity');
box on; grid on;

subplot(3,1,3);
text(0.5,0.5,sprintf('Total classified = %d', totalnumber),...
    'FontSize',12,'HorizontalAlignment','center');
title('Classification Summary');
axis off;

%% 8. Overexposed removal after correction
over_mask2 = max(corrected_spectra, [], 1) >= OVEREXPOSED_THRESHOLD - 1;
filtered = corrected_spectra(:, ~over_mask2);
fprintf('Removed %d overexposed spectra (post-correction).\n', sum(over_mask2));

%% 9. Statistical analysis & save
initial_mean = mean(filtered,2);
initial_std = std(filtered,0,2);
z_threshold = 1.96;
inlier = abs(filtered - initial_mean) <= z_threshold * initial_std;
filtered(~inlier) = NaN;
mean_spectrum = nanmean(filtered,2);
std_spectrum = nanstd(filtered,0,2);

writematrix([raman_shift mean_spectrum std_spectrum],...
    [base '_Mean_STD_Spectra.xlsx']);
fprintf('Saved mean and std spectra.\n');

%% 10. Plot mean ± SD
figure('Color','w','Position',[100,100,800,600]);
xfill = [raman_shift; flipud(raman_shift)];
yfill = [mean_spectrum + std_spectrum; flipud(mean_spectrum - std_spectrum)];
fill(xfill, yfill, [0.8 0.8 1], 'EdgeColor','none'); hold on;
plot(raman_shift, mean_spectrum, 'b','LineWidth',2);
title('Mean Corrected Spectrum with Standard Deviation');
xlabel('Raman Shift (cm^{-1})');
ylabel('Intensity (a.u.)');
legend('Mean ± 1 SD','Mean Spectrum','Location','northeast');
grid on; box on;

ax = gca;
ax.FontSize = 14;
ax.LineWidth = 1.5;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
set(gca,'MinorGridLineStyle',':','GridColor',[0.9 0.9 0.9],...
    'GridAlpha',0.5,'MinorGridAlpha',0.2);

%% Sub-function: airPLS baseline correction
function [Xc, Z] = airPLS(X, lambda, order, wep, p, itermax)
    if nargin < 6, itermax = 20; end
    if nargin < 5, p = 0.05; end
    if nargin < 4, wep = 0.1; end
    if nargin < 3, order = 2; end
    if nargin < 2, lambda = 1e7; end

    [m, n] = size(X);
    wi = [1:ceil(n*wep) floor(n-n*wep):n];
    D = diff(speye(n), order);
    DD = lambda * (D' * D);
    Z = zeros(m, n);

    for i = 1:m
        w = ones(n,1);
        x = X(i,:);
        for j = 1:itermax
            W = spdiags(w,0,n,n);
            C = chol(W + DD);
            z = (C\(C'\(w.*x')))';
            d = x - z;
            dssn = abs(sum(d(d<0)));
            if dssn < 0.001*sum(abs(x)), break; end
            w(d>=0) = 0;
            w(wi) = p;
            w(d<0) = j * exp(abs(d(d<0)) / dssn);
        end
        Z(i,:) = z;
    end
    Xc = X - Z;
end