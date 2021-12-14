%%% Paul Adkisson
%%% 9.6.21
%%% Purpose: Calculate decision time and accuracy from population firing
%%% rates
sim_name = "Test";
sim_path = sprintf("Simulation %s", sim_name);
load(strcat(sim_path, "/bam_constants.mat"))
%t = 0:0.02e-3:t_span;
%stim_amps = [-550, 0, 550]*1e-12;
trials = start_trial:end_trial;
num_trials = length(trials);
for brain = brains
    for k = 1:length(stim_amps)
        stim_amp = stim_amps(k);
        pulse = k <= length(pulse_amps);
        if pulse
            output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
                [sim_name, brain, stim_amp*1e9]);
        else
            output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                [sim_name, brain, stim_amp*1e9]);
        end
        decisions = zeros(num_trials, length(coherences)); %0 for no decision, 1 for C1, 2 for C2
        decision_times = zeros(num_trials, length(coherences));
        avg_dts = zeros(length(coherences), 1);
        std_dts = zeros(length(coherences), 1);
        avg_correct_dts = zeros(length(coherences), 1);
        std_correct_dts = zeros(length(coherences), 1);
        avg_incorrect_dts = zeros(length(coherences), 1);
        std_incorrect_dts = zeros(length(coherences), 1);
        avg_acc = zeros(length(coherences), 1);
        percent_nodec = zeros(length(coherences), 1);
        percent_earlydec = zeros(length(coherences), 1);
        tot_frs = zeros(length(coherences), num_trials, length(t), 4);
        avg_correct_frs = zeros(length(coherences), 1, length(t), 4);
        std_correct_frs = zeros(length(coherences), 1, length(t), 4);
        avg_incorrect_frs = zeros(length(coherences), 1, length(t), 4);
        std_incorrect_frs = zeros(length(coherences), 1, length(t), 4); 
        avg_nodec_frs = zeros(length(coherences), 1, length(t), 4);
        std_nodec_frs = zeros(length(coherences), 1, length(t), 4);
        avg_earlydec_frs = zeros(length(coherences), 1, length(t), 4);
        std_earlydec_frs = zeros(length(coherences), 1, length(t), 4);
        avg_correct_peaks = zeros(length(coherences), 4);
        std_correct_peaks = zeros(length(coherences), 4);
        avg_incorrect_peaks = zeros(length(coherences), 4);
        std_incorrect_peaks = zeros(length(coherences), 4);
        avg_nodec_peaks = zeros(length(coherences), 4);
        std_nodec_peaks = zeros(length(coherences), 4);
        dec_thresholds_lb = zeros(length(coherences), num_trials);
        dec_thresholds_ub = zeros(length(coherences), num_trials);
        totspikes_g1 = zeros(length(coherences), num_group*num_trials);
        totspikes_g2 = zeros(length(coherences), num_group*num_trials);
        totspikes_ns = zeros(length(coherences), (N_E-2*num_group)*num_trials);
        totspikes_int = zeros(length(coherences), N_I*num_trials);
        for i = 1:length(coherences)
            c = coherences(i);
            fprintf("Coherence: %0.1f%% \n", c*100)
            output_coherentpath = strcat(output_stimpath, sprintf("/c=%0.3f", c));
            for trial = start_trial:end_trial
                output_trialpath = strcat(output_coherentpath, sprintf("/trial%0.0f.mat", trial));
                load(output_trialpath, "recspikes")
                relative_trial = trial - start_trial + 1;
                spikes = zeros(length(t), N);
                for nn = 1:N
                    for spike_idx = recspikes(int2str(nn))
                        spikes(spike_idx, nn) = 1;
                    end
                    if nn <= num_group
                        totspikes_g1(i, nn+(num_group)*(relative_trial-1)) = length(recspikes(int2str(nn))) / t_span;
                    elseif nn > num_group && nn <= 2*num_group
                        totspikes_g2(i, (nn - num_group)+(num_group)*(relative_trial-1)) = length(recspikes(int2str(nn))) / t_span;
                    elseif nn > 2*num_group && nn <= N_E
                        totspikes_ns(i, (nn - 2*num_group)+(N_E-2*num_group)*(relative_trial-1)) = length(recspikes(int2str(nn))) / t_span;
                    else
                        totspikes_int(i, (nn - N_E)+N_I*(relative_trial-1)) = length(recspikes(int2str(nn))) / t_span;
                    end
                end
                [pop_frs, fr_vars] = spikes2popfrs(spikes, dt, p, f, N_E);
                save(output_trialpath, "pop_frs", "fr_vars", "recspikes")
                [decision, decision_idx] = get_decision(pop_frs);
                decisions(relative_trial, i) = decision;
                if decision_idx == -1
                    decision_times(relative_trial, i) = Inf;
                else
                    decision_times(relative_trial, i) = t(decision_idx) - t_task;
                end
                tot_frs(i, relative_trial, :, :) = pop_frs;
            end
            %omit trials with no decision & those that decide before onset of task-related input
            coherent_times = decision_times(decisions(:, i)~=0 & decision_times(:, i)>0, i);
            avg_dts(i) = mean(coherent_times);
            std_dts(i) = std(coherent_times);
            %omit trials with no decision & those that decide before onset of task-related input
            coherent_decisions = decisions(decisions(:, i)~=0 & decision_times(:, i)>0, i);
            avg_acc(i) = sum(coherent_decisions == 1) / length(coherent_decisions);
            percent_nodec(i) = sum(decisions(:, i)==0, 'all') / end_trial;
            percent_earlydec(i) = sum(decision_times(:, i)<=0, 'all') / end_trial;

            %Breakdown DTs by outcome
            avg_correct_dts(i) = mean(coherent_times(coherent_decisions==1));
            std_correct_dts(i) = std(coherent_times(coherent_decisions==1));
            avg_incorrect_dts(i) = mean(coherent_times(coherent_decisions==2));
            std_incorrect_dts(i) = std(coherent_times(coherent_decisions==2));

            %Breakdown by outcome (correct, incorrect, no-decision, early-decision)
            correct_frs = tot_frs(i, decisions(:, i)==1, :, :);
            incorrect_frs = tot_frs(i, decisions(:, i)==2, :, :);
            nodec_frs = tot_frs(i, decisions(:, i)==0, :, :);
            earlydec_frs = tot_frs(i, decision_times(:, i)<=0, :, :); 
            avg_correct_frs(i, :, :, :) = mean(correct_frs, 2);
            std_correct_frs(i, :, :, :) = std(correct_frs, 0, 2);
            avg_incorrect_frs(i, :, :, :) = mean(incorrect_frs, 2);
            std_incorrect_frs(i, :, :, :) = std(incorrect_frs, 0, 2);
            avg_incorrect_frs(i, :, :, :) = mean(incorrect_frs, 2);
            std_incorrect_frs(i, :, :, :) = std(incorrect_frs, 0, 2);
            avg_nodec_frs(i, :, :, :) = mean(nodec_frs, 2);
            std_nodec_frs(i, :, :, :) = std(nodec_frs, 0, 2);
            avg_earlydec_frs(i, :, :, :) = mean(earlydec_frs, 2);
            std_nodec_frs(i, :, :, :) = std(earlydec_frs, 0, 2);

            %Record peak frs split by outcome
            peak_frs = max(tot_frs(i, :, :, :), [], 3);
            correct_peaks = peak_frs(:, decisions(:, i)==1, :, :);
            incorrect_peaks = peak_frs(:, decisions(:, i)==2, :, :);
            nodec_peaks = peak_frs(:, decisions(:, i)==0, :, :);
            if ~isempty(correct_peaks)
                avg_correct_peaks(i, :) = mean(correct_peaks, 2);
                std_correct_peaks(i, :) = std(correct_peaks, 0, 2);
            end
            if ~isempty(incorrect_peaks)
                avg_incorrect_peaks(i, :) = mean(incorrect_peaks, 2);
                std_incorrect_peaks(i, :) = std(incorrect_peaks, 0, 2);
            end
            if ~isempty(nodec_peaks)
                avg_nodec_peaks(i, :) = mean(nodec_peaks, 2);
                std_nodec_peaks(i, :) = std(nodec_peaks, 0, 2);
            end

            %Decision Thresholds
            dec_thresholds_lb(i, decisions(:, i)==1) = max(tot_frs(i, decisions(:, i)==1, :, 2), [], 3);
            dec_thresholds_lb(i, decisions(:, i)==2) = max(tot_frs(i, decisions(:, i)==2, :, 1), [], 3);
            dec_thresholds_ub(i, decisions(:, i)==1) = max(tot_frs(i, decisions(:, i)==1, :, 1), [], 3);
            dec_thresholds_ub(i, decisions(:, i)==2) = max(tot_frs(i, decisions(:, i)==2, :, 2), [], 3);
            dec_thresholds_ub(i, decisions(:, i)==0) = max(dec_thresholds_ub(i, :), [], 'all'); %ensure that no-dec trials are not affecting threshold upper bound
        end

        %Pre-stimulus bias
        prestim_bias = tot_frs(:, :, abs(t-t_task)<dt/2, 1) - tot_frs(:, :, abs(t-t_task)<dt/2, 2);

        %Decision Thresholds
        dec_thresh = zeros(2, 1);
        dec_thresh_idx = zeros(2, 2);
        [dec_thresh(1), idx] = max(dec_thresholds_lb, [], 'all', 'linear');
        [c_idx, trial_idx] = ind2sub(size(dec_thresholds_lb), idx);
        dec_thresh_idx(1, :) = [c_idx, trial_idx];
        [dec_thresh(2), idx] = min(dec_thresholds_ub, [], 'all', 'linear');
        [c_idx, trial_idx] = ind2sub(size(dec_thresholds_ub), idx);
        dec_thresh_idx(2, :) = [c_idx, trial_idx];

        %Fit to weibull FN
        %coeffs = lsqcurvefit(@weibull, [1e-9, 0], coherences, avg_acc');
        coeffs = [1e-9, 0];
        avg_frs = mean(tot_frs, 2);
        std_frs = std(tot_frs, 0, 2);
        decisionpath = strcat(output_stimpath, "/decisions.mat");
        save(decisionpath, "decisions", "decision_times", ...
            "avg_dts", "std_dts", "avg_acc", "percent_nodec", "coeffs", ...
            "avg_frs", "std_frs", "avg_correct_frs", "std_correct_frs", ...
            "avg_incorrect_frs", "std_incorrect_frs", ...
            "avg_nodec_frs", "std_nodec_frs", ...
            "percent_earlydec", "avg_earlydec_frs", "std_earlydec_frs", ...
            "avg_correct_dts", "std_correct_dts", "avg_incorrect_dts", "std_incorrect_dts", ...
            "avg_correct_peaks", "std_correct_peaks", ...
            "avg_incorrect_peaks", "std_incorrect_peaks", ...
            "avg_nodec_peaks", "std_nodec_peaks", "prestim_bias", ...
            "dec_thresh", "dec_thresh_idx", "dec_thresholds_lb", "dec_thresholds_ub", ...
            "totspikes_g1", "totspikes_g2", "totspikes_ns", "totspikes_int")
    end
end

function [decision, decision_idx] = get_decision(pop_frs)
    decision_thresh = 15; %Hz
    decision1_idx = find(pop_frs(:, 1)>=decision_thresh, 1);
    decision2_idx = find(pop_frs(:, 2)>=decision_thresh, 1);
    if isempty(decision1_idx) && isempty(decision2_idx) %no decision
        decision = 0;
        decision_idx = -1;
    elseif isempty(decision1_idx)
        decision = 2;
        decision_idx = decision2_idx;
    elseif isempty(decision2_idx) || decision1_idx <= decision2_idx
        decision = 1;
        decision_idx = decision1_idx;
    else
        decision = 2;
        decision_idx = decision2_idx;
    end
end

function [frs, fr_var] = spikes2popfrs(spikes, dt, p, f, N_E)
    win_size = 5e-3;
    avg_win_size = 50e-3;
    num_group = floor(f*N_E);
    w = ones(floor(win_size/dt), 1);
    w = w ./ length(w);
    neuron_frs = filter(w, 1, spikes) ./ dt;
    w = ones(floor(avg_win_size/dt), 1);
    w = w ./ length(w);
    neuron_frs = filter(w, 1, neuron_frs);
    frs = zeros(size(spikes, 1), p+2);
    fr_var = zeros(size(spikes, 1), p+2);
    for i = 1:p
        group_idx = ((i-1)*num_group+1):i*num_group;
        frs(:, i) = mean(neuron_frs(:, group_idx), 2);
        fr_var(:, i) = std(neuron_frs(:, group_idx), [], 2);
    end
    frs(:, end-1) = mean(neuron_frs(:, f*p*N_E+1:N_E), 2);
    fr_var(:, end-1) = std(neuron_frs(:, f*p*N_E+1:N_E), [], 2);
    frs(:, end) = mean(neuron_frs(:, N_E+1:end), 2);
    fr_var(:, end) = std(neuron_frs(:, N_E+1:end), [], 2);
end