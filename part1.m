%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSP Project 2021
% MATLAB Exercise
% Part 1
% Dvir Ben Asuli, Assaf Gadish
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_a_mone = [1 -0.4 -0.69 -0.964 -0.714];
x_b_mehane = [1 0.6 0.34 -0.034 -0.3807 -0.4212];

%%% Question 1
[Sxx_a, Sxx_b] = question1(x_a_mone, x_b_mehane);

%%% Question 2
question2(x_a_mone, x_b_mehane);

%%% Question 3
%question3(x_a_mone, x_b_mehane, Sxx_a, Sxx_b);

%%% Question 4
%question4(x_a_mone, x_b_mehane);



function [] = question4(x_a_mone, x_b_mehane)
    question4a(x_a_mone, x_b_mehane, 1000)
end

function [] = question4a(x_a_mone, x_b_mehane, M)
    K = 1024;
    N = 400;

    sxxa_array = zeros(M,K);
    sxxb_array = zeros(M,K);

    [sxxa, wa, ~, ~] = calc_real_spectrum(x_a_mone, [1], K);
    [sxxb, wb, ~, ~] = calc_real_spectrum([1], x_b_mehane, K);

    for m = (1:M)
        [xa, xb] = generate_xa_xb(x_a_mone, x_b_mehane);
        Xa_k = get_positive_fft(xa, K);
        Sxxa_pos_periodogram = 1 / N * (Xa_k .* conj(Xa_k));
        sxxa_array(m,:) = Sxxa_pos_periodogram;

        Xb_k = get_positive_fft(xb, K);
        Sxxb_pos_periodogram = 1 / N * (Xb_k .* conj(Xb_k));
        sxxb_array(m,:) = Sxxb_pos_periodogram;
    end

    [sxxa_hat, Ba, std_a, RMSE_a] = calc_empiric_values(sxxa_array, sxxa);
    [sxxb_hat, Bb, std_b, RMSE_b] = calc_empiric_values(sxxb_array, sxxb);

    plot_empiric_values(sxxa, Ba, std_a, RMSE_a, 'Periodogram for x_a[n]');
    plot_empiric_values(sxxb, Bb, std_b, RMSE_b, 'Periodogram for x_b[n]');

end

function [sxx_hat, B, std, RMSE] = calc_empiric_values(sxx_array, sxx)
    sxx_hat = empiric_average(sxx_array);
    B = empiric_bias(sxx_hat, sxx);
    V = empiric_variance(sxx_hat, sxx_array);
    std = V.^0.5;
    MSE = empiric_mse(sxx_array, sxx);
    RMSE = MSE.^0.5;
end

function [] = plot_empiric_values(sxx, B, std, RMSE, varname)
    figure;
    title(strcat('Performance:', {' '}, varname, '[n]'))
    hold on
    plot_positive_spectrum_graph_only(sxx);
    plot_positive_spectrum_graph_only(B);
    plot_positive_spectrum_graph_only(std);
    plot_positive_spectrum_graph_only(RMSE);
    hold off
    legend('sxx', 'Bias', 'std', 'RMSE')
end

function [sxx_hat] = empiric_average(sxx_array)
    M = size(sxx_array, 1);
    sxx_hat = sum(sxx_array)';
    sxx_hat = sxx_hat / M;
end

function [B] = empiric_bias(sxx_hat, sxx)
    B = sxx_hat - sxx;
end

function [V] = empiric_variance(sxx_hat, sxx_array)
    M = size(sxx_array, 1);
    v_array = (sxx_array - sxx_hat').^2;
    V = sum(v_array)';
    V = V/ M;
end

function [MSE] = empiric_mse(sxx_array, sxx)
    M = size(sxx_array, 1);
    v_array = (sxx_array - sxx').^2;
    V = sum(v_array)';
    MSE = V/ M;
end

function [] = question3(x_a_mone, x_b_mehane, Sxx_a, Sxx_b)
    % Section a: biased correlogram estimator
    [xa, xb] = generate_xa_xb(x_a_mone, x_b_mehane);
    question3_parts_abc(xa, Sxx_a, 'x_a');
    question3_parts_abc(xb, Sxx_b, 'x_b');
end

function [] = question3_parts_abc(x, Sxx, varname)
    N = 400;
    K = 1024;
    Rxx_biased_correlogram = estimate_biased_correlogram(x, N);
    Sxx_pos_biased_correlogram = get_positive_fft(Rxx_biased_correlogram , K);
    figure;
    plot_positive_spectrum(Sxx_pos_biased_correlogram, 'S_{xx}^B', varname);

    % Section b: periodogram estimator
    X_k = get_positive_fft(x, K);
    Sxx_pos_periodogram = 1 / N * (X_k .* conj(X_k));
    figure;
    plot_positive_spectrum(Sxx_pos_periodogram, 'S_{xx}^P', varname);

    % Section c: a,b with original
    figure;
    % 'S_{xx}^B'
    title(strcat('S_{xx} for', {' '}, varname, '[n]'))
    hold on
    plot_positive_spectrum_graph_only(Sxx_pos_biased_correlogram);
    % 'S_{xx}^P'
    plot_positive_spectrum_graph_only(Sxx_pos_periodogram);
    plot_positive_spectrum_graph_only(Sxx);
    hold off
    legend('biased correlogram', 'periodogram', 'true spectrum')

    % Section d: unbiased correlogram
    Rxx_unbiased_correlogram = estimate_unbiased_correlogram(x, N);
    Sxx_pos_unbiased_correlogram = get_positive_fft(Rxx_unbiased_correlogram , K);
    figure;
    plot_positive_spectrum(Sxx_pos_unbiased_correlogram, 'S_{xx}', varname);

    figure;
    % 'S_{xx}^B'
    title(strcat('S_{xx} for', {' '}, varname, '[n]'));
    hold on;
    plot_positive_spectrum_graph_only(Sxx_pos_biased_correlogram);
    % 'S_{xx}^P'
    plot_positive_spectrum_graph_only(Sxx_pos_unbiased_correlogram);
    plot_positive_spectrum_graph_only(Sxx_pos_periodogram);

    plot_positive_spectrum_graph_only(Sxx);
    hold off;
    legend('biased correlogram', 'unbiased correlogram', 'periodogram', 'true spectrum')


end

function [] = plot_positive_spectrum(Sxx, s_name, varname)
    plot_positive_spectrum_graph_only(Sxx);
    title(strcat('|', s_name, '(e^{j\omega})|(\omega) for', {' '}, varname))
    xlabel('\omega [rad]');
    ylabel(strcat('|', s_name, '(e^{j\omega})| [dB] for', {' '}, varname));
end

function [] = plot_positive_spectrum_graph_only(Sxx)
    w = linspace(0, pi, length(Sxx));
    plot(w, abs(Sxx));
end


function ffted_x_positive = get_positive_fft(x, K)
    x_padded = [x' zeros(1, 2 * K - size(x, 1))];
    ffted_x = fftshift(fft(x_padded));
    ffted_x_positive = ffted_x(K + 1 : 2 * K)';
end

function Rxx = estimate_biased_correlogram(x, N)
    Rxx = zeros((2 * N) - 1, 1);
    for l = (0: N - 1)
        r = estimate_biased_correlogram_l(x, l, N);
        Rxx(N - l) = r;
        Rxx(N + l) = r;
    end
end

function Rxx_l = estimate_biased_correlogram_l(x, l, N)
    abs_l = abs(l);
    Rxx_l = 0;
    for n = (1 : N - abs_l)
        Rxx_l = Rxx_l + x(n) * x(n + abs_l);
    end
    Rxx_l = Rxx_l / N;
end    

function Rxx = estimate_unbiased_correlogram(x, N)
    Rxx = zeros((2 * N) - 1, 1);
    for l = (0: N - 1)
        r = estimate_unbiased_correlogram_l(x, l, N);
        Rxx(N - l) = r;
        Rxx(N + l) = r;
    end
end

function Rxx_l = estimate_unbiased_correlogram_l(x, l, N)
    abs_l = abs(l);
    Rxx_l = 0;
    for n = (1 : N - abs_l)
        Rxx_l = Rxx_l + x(n) * x(n + abs_l);
    end
    Rxx_l = Rxx_l / (N - abs_l);
end

function [sxx, w, sxx_zeros, sxx_poles] = calc_real_spectrum(x_mone, x_mehane, K)
    x_mone_roots = roots(x_mone);
    x_mehane_roots = roots(x_mehane);
    sxx_zeros = add_conj_inverse_roots(x_mone_roots);
    sxx_poles = add_conj_inverse_roots(x_mehane_roots);
    sxx_mone = poly(sxx_zeros);
    sxx_mehane = poly(sxx_poles);
    [sxx, w] = freqz(sxx_mone, sxx_mehane, K);
end

function [sxxa, sxxb] = question1(x_a_mone, x_b_mehane)
    K = 1024;
    % Question 1 section a
    [sxxa, wa, sxxa_zeros, ~] = calc_real_spectrum(x_a_mone, [1], K);
    fprintf('S_xx_a zeros:\n');
    figure;
    disp(sxxa_zeros)
    zplane(sxxa_zeros, []);
    title('Zplane of the MA(4) process, x_{a}[n]');
    
    [sxxb, wb, ~, sxxb_poles] = calc_real_spectrum([1], x_b_mehane, K);

    fprintf('s_xx_b poles:\n');
    figure;
    disp(sxxb_poles)
    zplane([], sxxb_poles);
    title('Zplane of the AR(5) process, x_{b}[n]');


    % Question 1 section b    
    figure;
    plot(wa, abs(sxxa));
    title('Spectrum of the MA(4) process, x_{a}[n]')
    xlabel('\omega [rad]');
    ylabel('|S_{xx_a}(e^{j\omega})|');
    
    figure;
    plot(wb, abs(sxxb));
    title('Spectrum of the AR(5) process, x_{b}[n]')
    xlabel('\omega [rad]');
    ylabel('|S_{xx_b}(e^{j\omega})|');
end

function [] = question2(x_a_mone, x_b_mehane)
    [xa, xb] = generate_xa_xb(x_a_mone, x_b_mehane);
        
    figure;
    plot(xa);
    title('x_a[n]');
    xlabel('n');
    
    figure;
    plot(xb);
    title('x_b[n]');
    xlabel('n');
end

function [xa, xb] = generate_xa_xb(x_a_mone, x_b_mehane)
    wa = randn(404, 1);
    wb = randn(2000, 1);

    xa_padded = filter(x_a_mone, 1, wa);
    xa = xa_padded(5:404);
    xb_padded = filter(1, x_b_mehane, wb);
    xb = xb_padded(1601:2000);
end

function [all_roots] = add_conj_inverse_roots(roots)
    conj_inverse_roots = 1 ./ conj(roots);
    all_roots = transpose([roots' conj_inverse_roots']);
end




