%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSP Project 2021
% MATLAB Exercise
% Part 1
% Dvir Ben Asuli, Assaf Gadish
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
x_a_mone = [1 -0.4 -0.69 -0.964 -0.714];
x_b_mehane = [1 0.6 0.34 -0.034 -0.3807 -0.4212];

%%% Question 1
%[Sxx_a, Sxx_b] = question1(x_a_mone, x_b_mehane);

%%% Question 2
%question2(x_a_mone, x_b_mehane);

%%% Question 3
%question3(x_a_mone, x_b_mehane, Sxx_a, Sxx_b);

%%% Question 4
 question4(x_a_mone, x_b_mehane);

%%% Question 5
%question5(x_a_mone, x_b_mehane);

%%% Question 6
% question6(x_a_mone, x_b_mehane);

function [] = question6_estimate_ma(x, sxx, signal_name)
    rxx = estimate_unbiased_correlogram(x, size(x, 1));
    rxx = zero_pad_center(rxx, 6);
    x_ma_3 = question6_calc_maq(rxx, 3);
    x_ma_4 = question6_calc_maq(rxx, 4);
    x_ma_5 = question6_calc_maq(rxx, 5);

    figure;
    title(strcat('Estimates Sxx using MA(q) models for ', {' '}, signal_name));
    hold on
    plot_positive_spectrum_graph_only(x_ma_3);
    plot_positive_spectrum_graph_only(x_ma_4);
    plot_positive_spectrum_graph_only(x_ma_5);
    plot_positive_spectrum_graph_only(sxx);
    hold off
    legend('q=3', 'q=4', 'q=5', 'Sxx')

end

function [] = question6_estimate_ar(x, sxx, signal_name)
    rxx = estimate_unbiased_correlogram(x, size(x, 1));
    N = size(x, 1);
    K = 1024;

    figure;
    title(strcat('Estimates Sxx using AR(p) models for ', {' '}, signal_name));
    hold on;
    for p = [4 5 6]
        toeplitz_A = toeplitz(rxx(N : N + p - 1));
        results_B = -(rxx(N + 1 : N + p));
        a_estimated = [1; linsolve(toeplitz_A, results_B)];
%% TODO: Figure out why the graph is too scaled up. Maybe wrong factor
        sigma_w_square = a_estimated' * rxx(N : N + p);
        ak_fft = get_positive_fft(a_estimated, K);
        sxx_hat = sigma_w_square ./ (ak_fft .* conj(ak_fft));
        plot_positive_spectrum_graph_only(sxx_hat);
    end

    plot_positive_spectrum_graph_only(sxx);
    hold off;
    legend('p=3', 'p=4', 'p=5', 'Sxx');

end

function [] = question6(x_a_mone, x_b_mehane)
    K = 1024;
    [xa, xb] = generate_xa_xb(x_a_mone, x_b_mehane);
    [sxxa, ~, ~, ~] = calc_real_spectrum(x_a_mone, [1], K);
    [sxxb, ~, ~, ~] = calc_real_spectrum([1], x_b_mehane, K);

    question6_estimate_ma(xa, sxxa, 'x_a[n]');
    question6_estimate_ma(xb, sxxb, 'x_b[n]');

    question6_estimate_ar(xa, sxxa, 'x_a[n]');
    question6_estimate_ar(xb, sxxb, 'x_b[n]');
end

function [cropped_x] = zero_pad_center(x, radius)
    n = size(x, 1);
    nc = (n + 1) / 2;
    cropped_x = [zeros(1, nc - radius) x(nc - radius : nc + radius)' zeros(1, nc - radius)]';
end

function [Sxx] = question6_calc_maq(r_xx, q)
    r_xx = zero_pad_center(r_xx, q);
    Sxx = get_positive_fft(r_xx, 1024);
end

function [b_avg, v_avg, mse_avg] = question5_get_avarages(b, v, mse)
    b_avg = norm(b) / 1024;
    v_avg = mean(v);
    mse_avg = mean(mse);
end

function [] = question5_print_stuff(name, Ba, stda, MSEa, Bb, stdb, MSEb)
    [Ba, stda, MSEa] = question5_get_avarages(Ba, stda, MSEa);
    [Bb, stdb, MSEb] = question5_get_avarages(Bb, stdb, MSEb);
    fprintf("%s: x_a[n] got <B^2>=%d, <V>=%d, <MSE>=%d\n", name, Ba, stda, MSEa);
    fprintf("%s: x_b[n] got <B^2>=%d, <V>=%d, <MSE>=%d\n", name, Bb, stdb, MSEb);
end

function [] = question5(x_a_mone, x_b_mehane)
% Section A
    [~, Baa, stdaa, MSEaa, ~, Bba, stdba, MSEba] = question4_periodogram(x_a_mone, x_b_mehane);
    question5_print_stuff("Periodogram", Baa, stdaa, MSEaa, Bba, stdba, MSEba);
%     [Baa, stdaa, MSEaa] = question5_get_avarages(Baa, stdaa, MSEaa);
%     [Bba, stdba, MSEba] = question5_get_avarages(Bba, stdba, MSEba);
    % Section B
    [~, Bab, stdab, MSEab, ~, Bbb, stdbb, MSEbb] = question4_bartlett(x_a_mone, x_b_mehane, 5, 80);
    question5_print_stuff("Bartlett", Bab, stdab, MSEab, Bbb, stdbb, MSEbb);
%     [Bab, stdab, MSEab] = question5_get_avarages(Bab, stdab, MSEab);
%     [Bbb, stdbb, MSEbb] = question5_get_avarages(Bbb, stdbb, MSEbb);
    
    % Section C
    [~, Bac, stdac, MSEac, ~, Bbc, stdbc, MSEbc] = question4_bartlett(x_a_mone, x_b_mehane, 10, 40);
%     [Bac, stdac, MSEac] = question5_get_avarages(Bac, stdac, MSEac);
%     [Bbc, stdbc, MSEbc] = question5_get_avarages(Bbc, stdbc, MSEbc);
    question5_print_stuff("Bartlett", Bac, stdac, MSEac, Bbc, stdbc, MSEbc);

    % Section D
    [~, Bad, stdad, MSEad, ~, Bbd, stdbd, MSEbd] = question4_welsh(x_a_mone, x_b_mehane, 80, 40);
%     [Bad, stdad, MSEad] = question5_get_avarages(Bad, stdad, MSEad);
%     [Bbd, stdbd, MSEbd] = question5_get_avarages(Bbd, stdbd, MSEbd);
    question5_print_stuff("Welsh", Bad, stdad, MSEad, Bbd, stdbd, MSEbd);

    % Section E
    [~, Bae, stdae, MSEae, ~, Bbe, stdbe, MSEbe] = question4_welsh(x_a_mone, x_b_mehane, 40, 20);
%     [Bae, stdae, MSEae] = question5_get_avarages(Bae, stdae, MSEae);
%     [Bbe, stdbe, MSEbe] = question5_get_avarages(Bbe, stdbe, MSEbe);
    question5_print_stuff("Welsh", Bae, stdae, MSEae, Bbe, stdbe,MSEbe);

    % Section F
    [~, Baf, stdaf, MSEaf, ~, Bbf, stdbf, MSEbf] = question4_blackman_tukey(x_a_mone, x_b_mehane, 40);
%     [Baf, stdaf, MSEaf] = question5_get_avarages(Baf, stdaf, MSEaf);
%     [Bbf, stdbf, MSEbf] = question5_get_avarages(Bbf, stdbf, MSEbf);
    question5_print_stuff("Blackman-Tukey", Baf, stdaf, MSEaf, Bbf, stdbf,MSEbf);

    % Section G
    [~, Bag, stdag, MSEag, ~, Bbg, stdbg, MSEbg] = question4_blackman_tukey(x_a_mone, x_b_mehane, 20);
%     [Bag, stdag, MSEag] = question5_get_avarages(Bag, stdag, MSEag);
%     [Bbg, stdbg, MSEbg] = question5_get_avarages(Bbg, stdbg, MSEbg);
    question5_print_stuff("Blackman-Tukey", Bag, stdag, MSEag, Bbg, stdbg,MSEbg);

end

function [] = question4(x_a_mone, x_b_mehane)
    % Section A
%     question4_periodogram(x_a_mone, x_b_mehane)

    % Section B
%     question4_bartlett(x_a_mone, x_b_mehane, 5, 80);
    
    % Section C
%     question4_bartlett(x_a_mone, x_b_mehane, 10, 40);

    % Section D
%     question4_welsh(x_a_mone, x_b_mehane, 80, 40);

    % Section E
%     question4_welsh(x_a_mone, x_b_mehane, 40, 20);

    % Section F
    question4_blackman_tukey(x_a_mone, x_b_mehane, 40);

    % Section G
    question4_blackman_tukey(x_a_mone, x_b_mehane, 20);
end

function wx = multiply_with_bartlett_window(x)
    N = size(x, 1);
    l = [1 : N];
    window_vector = 1 - (l / N);
    wx = (window_vector') .* x;
end

% function [sxx_bt] = calc_blackman_tukey(x, L)
%     K = 1024;
%     N = size(x, 1);
% 
%     Xk = get_positive_fft(x, K);
%     Sxxa_pos_periodogram = 1 / N * (Xk .* conj(Xk));
%     
%     w = (1 : N);
%     Window = sin(w ./ 2 * (2 * L + 1)) ./ sin(w ./2);
%     sxx_bt_larger = 1 / (2 * pi) * conv(Sxxa_pos_periodogram, Window);
%     sxx_bt = sxx_bt_larger(1:1024);
% %     sxx_bt = abs(get_positive_fft(x, K)) / (2 * pi);
% end

function [sxx_bt] = calc_blackman_tukey(x, L)
    K = 1024;
    N = size(x, 1);
    Rxx_correlogram = estimate_biased_correlogram(x, N);
    % Calculate Rxx barlett
    rxx_bt = multiply_with_bartlett_window(Rxx_correlogram);
%     rxx_bt = Rxx_correlogram;
    % Multiply Rxx with rectangle window
    rxx_bt = [zeros(1, N - L - 1) rxx_bt(N - L : N + L)' zeros(1, N - L - 1)]';
    sxx_bt = abs(get_positive_fft(rxx_bt, K));
end

function [sxxa_hat, Ba, std_a, RMSE_a, sxxb_hat, Bb, std_b, RMSE_b] = question4_blackman_tukey(x_a_mone, x_b_mehane, L)
    K = 1024;
    N = 400;
    M = 1000;

    sxxa_array = zeros(M,K);
    sxxb_array = zeros(M,K);

    [sxxa, wa, ~, ~] = calc_real_spectrum(x_a_mone, [1], K);
    [sxxb, wb, ~, ~] = calc_real_spectrum([1], x_b_mehane, K);

    for m = (1:M)
        [xa, xb] = generate_xa_xb(x_a_mone, x_b_mehane);
        
        sxxa_blackman_tukey = calc_blackman_tukey(xa, L);
        sxxa_array(m,:) = sxxa_blackman_tukey;

        sxxb_blackman_tukey = calc_blackman_tukey(xb, L);
        sxxb_array(m,:) = sxxb_blackman_tukey;
    end

    [sxxa_hat, Ba, std_a, RMSE_a] = calc_empiric_values(sxxa_array, sxxa);
    [sxxb_hat, Bb, std_b, RMSE_b] = calc_empiric_values(sxxb_array, sxxb);

    plot_empiric_values(sxxa, Ba, std_a, RMSE_a, sprintf('Blackman-Tukey for x_a[n] with L=%d', L));
    plot_empiric_values(sxxb, Bb, std_b, RMSE_b, sprintf('Blackman-Tukey for x_b[n] with L=%d', L));
end

function [sxx_welsh] = calc_welsh(x, L, D)
    N_SAMPLES = 1024;
    K = floor((size(x, 1) - L) / D) + 1;
    N = size(x, 1);
    all_sxx_welshes = zeros(K, N_SAMPLES);
    % Go over the k-th block (size L)
    for k = (1 : K)
        x_kth_block = zeros(1024, 1);
        x_kth_block(1:L) = x((k - 1) * D + 1 : (k - 1) * D + L);
        X_k_L = get_positive_fft(x_kth_block, 1024);
        current_sxx_welsh = 1 / L * (X_k_L .* conj(X_k_L));
        all_sxx_welshes(k, :) = current_sxx_welsh;
    end

    sxx_welsh = empiric_average(all_sxx_welshes);
end

function [sxxa_hat, Ba, std_a, RMSE_a, sxxb_hat, Bb, std_b, RMSE_b] = question4_welsh(x_a_mone, x_b_mehane, L, D)
    K = 1024;
    N = 400;
    M = 1000;

    sxxa_array = zeros(M,K);
    sxxb_array = zeros(M,K);

    [sxxa, wa, ~, ~] = calc_real_spectrum(x_a_mone, [1], K);
    [sxxb, wb, ~, ~] = calc_real_spectrum([1], x_b_mehane, K);

    for m = (1:M)
        [xa, xb] = generate_xa_xb(x_a_mone, x_b_mehane);
        
        sxxa_welsh = calc_welsh(xa, L, D);
        sxxa_array(m,:) = sxxa_welsh;

        sxxb_welsh = calc_welsh(xb, L, D);
        sxxb_array(m,:) = sxxb_welsh;
    end

    [sxxa_hat, Ba, std_a, RMSE_a] = calc_empiric_values(sxxa_array, sxxa);
    [sxxb_hat, Bb, std_b, RMSE_b] = calc_empiric_values(sxxb_array, sxxb);

    plot_empiric_values(sxxa, Ba, std_a, RMSE_a, sprintf('Welsh for x_a[n] with L=%d, D=%d', L, D));
    plot_empiric_values(sxxb, Bb, std_b, RMSE_b, sprintf('Welsh for x_b[n] with L=%d, D=%d', L, D));
end

function [sxx_bartlett] = calc_bartlett(x, K, L)
    N_SAMPLES = 1024;
    N = size(x, 1);
    assert(N <= K * L, "Bartlett errors: size of x must be at least K * L");
    all_sxx_bartletts = zeros(K, N_SAMPLES);
    % Go over the k-th block (size L)
    for k = (1 : K)
        x_kth_block = zeros(1024, 1);
        x_kth_block(1:L) = x((k - 1) * L + 1 : k * L);
        X_k_L = get_positive_fft(x_kth_block, 1024);
        current_sxx_bartlett = 1 / L * (X_k_L .* conj(X_k_L));
        all_sxx_bartletts(k, :) = current_sxx_bartlett;
    end

    sxx_bartlett = empiric_average(all_sxx_bartletts);
end


function [sxxa_hat, Ba, std_a, RMSE_a, sxxb_hat, Bb, std_b, RMSE_b] = question4_bartlett(x_a_mone, x_b_mehane, K_bartlett, L)
    K = 1024;
    N = 400;
    M = 1000;

    sxxa_array = zeros(M,K);
    sxxb_array = zeros(M,K);

    [sxxa, wa, ~, ~] = calc_real_spectrum(x_a_mone, [1], K);
    [sxxb, wb, ~, ~] = calc_real_spectrum([1], x_b_mehane, K);

    for m = (1:M)
        [xa, xb] = generate_xa_xb(x_a_mone, x_b_mehane);
        
        sxxa_bartlett = calc_bartlett(xa, K_bartlett, L);
        sxxa_array(m,:) = sxxa_bartlett;

        sxxb_bartlett = calc_bartlett(xb, K_bartlett, L);
        sxxb_array(m,:) = sxxb_bartlett;
    end

    [sxxa_hat, Ba, std_a, RMSE_a] = calc_empiric_values(sxxa_array, sxxa);
    [sxxb_hat, Bb, std_b, RMSE_b] = calc_empiric_values(sxxb_array, sxxb);

    plot_empiric_values(sxxa, Ba, std_a, RMSE_a, sprintf('Bartlett for x_a[n] with K=%d, L=%d', K_bartlett, L));
    plot_empiric_values(sxxb, Bb, std_b, RMSE_b, sprintf('Bartlett for x_b[n] with K=%d, L=%d', K_bartlett, L));

end

function [sxxa_hat, Ba, std_a, RMSE_a, sxxb_hat, Bb, std_b, RMSE_b] = question4_periodogram(x_a_mone, x_b_mehane)
    K = 1024;
    N = 400;
    M = 1000;

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
    plot_spectrum_graph_only(B);
    plot_spectrum_graph_only(std);
    plot_spectrum_graph_only(RMSE);
    hold off
    legend('sxx', 'bias', 'std', 'RMSE')
end

function [sxx_hat] = empiric_average(sxx_array)
    M = size(sxx_array, 1);
    sxx_hat = sum(sxx_array)';
    sxx_hat = sxx_hat / M;
end

function [B] = empiric_bias(sxx_hat, sxx)
    B = sxx_hat - abs(sxx);
end

function [V] = empiric_variance(sxx_hat, sxx_array)
    M = size(sxx_array, 1);
    v_array = (sxx_array - sxx_hat').^2;
    V = sum(v_array)';
    V = V/ M;
end

function [MSE] = empiric_mse(sxx_array, sxx)
    M = size(sxx_array, 1);
    v_array = (sxx_array - abs(sxx)').^2;
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

function [] = plot_spectrum_graph_only(Sxx)
    w = linspace(0, pi, length(Sxx));
    plot(w, Sxx);
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
    factor_mone = prod(abs(x_mone_roots)); 
    factor_mehane = prod(abs(x_mehane_roots)); 
    sxx_zeros = add_conj_inverse_roots(x_mone_roots);
    sxx_poles = add_conj_inverse_roots(x_mehane_roots);
    sxx_mone = poly(sxx_zeros) * factor_mone;
    sxx_mehane = poly(sxx_poles) * factor_mehane;
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




