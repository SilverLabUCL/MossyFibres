function [A_state, Q_state, RL_A_state, RL_Q_state, wsk_A_state, wsk_Q_state] = detect_state(RL_rs, wsk_rs, nnids) 
 
    RL_flat = RL_rs(:);
    wsk_flat = wsk_rs(:);
    N = numel(RL_flat);
    %nnids = nnids(nnids >= 1 & nnids <= N);
    valid_mask = true(N,1);
    valid_mask(nnids) = false;
    
    RL_A_state = zeros(N,1);
    RL_Q_state = zeros(N,1);
    wsk_A_state = zeros(N,1);
    wsk_Q_state = zeros(N,1);
    
    %% === RL===
    RL_valid = RL_flat(valid_mask);
    RL_std = std(RL_valid);
    RL_mean = mean(RL_valid);
    RL_median = median(RL_valid);

    RL_quiet_mask = valid_mask & (...
        (RL_flat < 0.03) | ...
        (RL_flat < RL_mean) | ...
        (RL_flat < RL_median));
    RL_Q_state(RL_quiet_mask) = 1;

    [pk_rl, loc_rl] = findpeaks(RL_flat, ...
        'MinPeakProminence', 0.2 * RL_std, ...
        'MinPeakHeight', RL_mean + 0.3 * RL_std, ...
        'MinPeakDistance', 5);

        for i = 1:length(loc_rl)
            val = pk_rl(i);
            th = 0.15 * val;
            l = loc_rl(i);
            r = loc_rl(i);
            while l > 1 && RL_flat(l) > th && valid_mask(l), l = l - 1; end
            l = l + 1;
            while r < N && RL_flat(r) > th && valid_mask(r), r = r + 1; end
            r = r - 1;
            RL_A_state(l:r) = 1;
        end
    
    %% === Whisker states detect
    wsk_valid = wsk_flat(valid_mask);
    wsk_std = std(wsk_valid);
    wsk_mean = mean(wsk_valid);
    wsk_median = median(wsk_valid);

    wsk_quiet_mask = valid_mask & (...
        (wsk_flat < 0.03) | ...
        (wsk_flat < wsk_mean) | ...
        (wsk_flat < wsk_median));
    wsk_Q_state(wsk_quiet_mask) = 1;

    [pk_wsk, loc_wsk] = findpeaks(wsk_flat, ...
        'MinPeakProminence', 0.1 * wsk_std, ...
        'MinPeakHeight', wsk_mean + 0.5 * wsk_std, ...
        'MinPeakDistance', 3);

    for i = 1:length(loc_wsk)
        val = pk_wsk(i);
        th = max(0.1 * val, wsk_mean + 0.1 * wsk_std);
        l = loc_wsk(i);
        r = loc_wsk(i); 
        while l > 1 && wsk_flat(l) > th && valid_mask(l), l = l - 1; end
        l = l + 1;
        while r < N && wsk_flat(r) > th && valid_mask(r), r = r + 1; end
        r = r - 1;
        wsk_A_state(l:r) = 1; 
    end

    % RL_A_state(RL_Q_state == 1) = 0;
    % wsk_A_state(wsk_Q_state == 1) = 0;

    A_state = zeros(N,1);
    Q_state = zeros(N,1);

    both_active = RL_A_state & wsk_A_state;
    CC_both_a = bwconncomp(both_active);
    for i = 1:CC_both_a.NumObjects
        idx = CC_both_a.PixelIdxList{i};
        if numel(idx) >= 300
            A_state(idx) = 1;
        end
    end

    both_quiet = RL_Q_state & wsk_Q_state & ~A_state;
    CC_both_q = bwconncomp(both_quiet);
    for i = 1:CC_both_q.NumObjects
        idx = CC_both_q.PixelIdxList{i};
        if numel(idx) >= 300
            Q_state(idx) = 1;
        end
    end

    %% === Reshape & Plot ===
    A_state = reshape(A_state, size(RL_rs));
    Q_state = reshape(Q_state, size(RL_rs));

    RL_A_state = reshape(RL_A_state, size(RL_rs));
    RL_Q_state = reshape(RL_Q_state, size(RL_rs));
    wsk_A_state = reshape(wsk_A_state, size(RL_rs));
    wsk_Q_state = reshape(wsk_Q_state, size(RL_rs));

    time = 1:N;
    figure;

    subplot(4,1,1); hold on;
    plot(time, RL_flat, 'k', 'LineWidth', 0.5);
    plot(time(RL_A_state(:)==1), RL_flat(RL_A_state(:)==1), 'r.', 'MarkerSize', 2);
    plot(time(RL_Q_state(:)==1), RL_flat(RL_Q_state(:)==1), 'b.', 'MarkerSize', 2);
    ylabel('RL'); title('RL signal: individual A/Q detection (no length filter)');
    legend('Raw', 'RL-Active', 'RL-Quiet', 'Location', 'best');

    subplot(4,1,2); hold on;
    plot(time, wsk_flat, 'k', 'LineWidth', 0.5);
    plot(time(wsk_A_state(:)==1), wsk_flat(wsk_A_state(:)==1), 'r.', 'MarkerSize', 2);
    plot(time(wsk_Q_state(:)==1), wsk_flat(wsk_Q_state(:)==1), 'b.', 'MarkerSize', 2);
    ylabel('Whisker'); title('Whisker signal: individual A/Q detection (no length filter)');
    legend('Raw', 'Wsk-Active', 'Wsk-Quiet', 'Location', 'best');

    subplot(4,1,3); hold on;
    plot(time, RL_flat, 'k', 'LineWidth', 0.5);
    plot(time, wsk_flat, 'g', 'LineWidth', 0.5);
    plot(time(A_state(:)==1), RL_flat(A_state(:)==1), 'r.', 'MarkerSize', 3);
    plot(time(Q_state(:)==1), RL_flat(Q_state(:)==1), 'b.', 'MarkerSize', 3);
    ylabel('Signals'); title('Final A/Q states (both signals agree & ≥300pts)');
    legend('RL', 'Whisker', 'Final-Active', 'Final-Quiet', 'Location', 'best');

    subplot(4,1,4);
    plot(time, A_state(:), 'r', 'LineWidth', 2); hold on;
    plot(time, Q_state(:), 'b', 'LineWidth', 2);
    xlabel('Time'); ylabel('State'); title('Final Binary A/Q states');
    legend('Active', 'Quiet', 'Location', 'best');
    ylim([-0.1, 1.1]);
end