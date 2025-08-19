%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sARA-TDC                                                                                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set random seed for reproducibility
rng(19);

% -------------------------------------------------------------------------
% SIMULATION OF RECURRENT EVENTS UNDER sARA1-TDC MODEL (no censoring)
% -------------------------------------------------------------------------

% --- User parameters ---
num_subjects  = 100;                 % Number of subjects
max_events    = 10;                  % Max number of events per subject
eta           = [0, 1, 2, 3, 4, 5];  % Interval boundaries
theta         = [0.5];               % Model parameter (adjust as needed)
rho_set       = [0, 0.7, 1];         % Intervention efficiencies
alpha         = 1.5;                 % Weibull PLP parameter
gamma         = 1.5;                 % Weibull PLP parameter
sigma2        = 0.1;                 % Variance for TDC covariate noise

% --- Baseline cumulative hazard and its inverse ---
Lambda_func     = @(v) alpha * v.^gamma;
Lambda_inv_func = @(x) (x/alpha).^(1/gamma);

% --- Virtual age process and its inverse (sARA1) ---
% V(t) = t - rho * last_event_time
% V_inv(v) = v + rho * last_event_time
V_func     = @(t, last_event, rho) t - rho * last_event;
V_inv_func = @(v, last_event, rho) v + rho * last_event;

% --- Storage ---
all_fail_times = cell(num_subjects, length(rho_set));
all_covariates = cell(num_subjects, length(rho_set));

% --- Simulation loop ---
for rr = 1:length(rho_set)
    rho = rho_set(rr);

    for i = 1:num_subjects
        fail_times = [0];
        covariates = [];
        Z0 = rand * 0.5;

        for j = 1:max_events
            T_j = fail_times(end);
            
            if j == 1
                last_event = 0;
            else
                last_event = fail_times(end-1);
            end

            % Time-dependent covariate
            cov_noise = sqrt(sigma2) * randn;
            Z_t = Z0 + cov_noise;
            covariates = [covariates; Z_t];

            % Step 1: Locate current interval
            l_j = find(T_j >= eta(1:end-1) & T_j < eta(2:end), 1);
            if isempty(l_j)
                break;
            end

            % Step 2: Target quantile
            U = rand;
            Q = -log(1 - U);

            % Step 3: Check if event stays in current interval
            z = theta' * Z_t;
            
            % V_i(t) * exp(theta'*Z_i(t))
            V_T_j = V_func(T_j, last_event, rho);
            V_eta_lj = V_func(eta(l_j + 1), last_event, rho);
            
            % Apply the transformation: V_i(t) * exp(theta'*Z_i(t))
            V_T_j_transformed    = V_T_j    * exp(z);
            V_eta_lj_transformed = V_eta_lj * exp(z);
            
            % Available hazard in current interval
            H_current = Lambda_func(V_eta_lj_transformed) - Lambda_func(V_T_j_transformed);

            if Q <= H_current
                % Event stays in current interval
                Lambda_sum         = Lambda_func(V_T_j_transformed) + Q;
                v_next_transformed = Lambda_inv_func(Lambda_sum);
                v_next             = v_next_transformed / exp(z);   % inverse transform
                T_next             = V_inv_func(v_next, last_event, rho);
            else
                % Step 4: Event crosses intervals
                I = Lambda_func(V_eta_lj_transformed) - Lambda_func(V_T_j_transformed);
                l = l_j + 1;
                L = length(eta) - 1;
                l_prime = []; z_l = [];

                while I < Q && l <= L
                    cov_noise_l = sqrt(sigma2) * randn;
                    Z_l = Z0 + cov_noise_l;
                    z_l = theta' * Z_l;
                    
                    V_eta_l_minus_1 = V_func(eta(l),   last_event, rho) * exp(z_l);
                    V_eta_l         = V_func(eta(l+1), last_event, rho) * exp(z_l);
                    
                    I = I + (Lambda_func(V_eta_l) - Lambda_func(V_eta_l_minus_1));
                    
                    if I >= Q
                        l_prime = l;
                        break;
                    end
                    l = l + 1;
                end
                
                if isempty(l_prime) && l > L
                    T_next = inf;
                else
                    V_eta_l_prime_minus_1 = V_func(eta(l_prime),   last_event, rho) * exp(z_l);
                    V_eta_l_prime         = V_func(eta(l_prime+1), last_event, rho) * exp(z_l);
                    
                    R = Q - (I - (Lambda_func(V_eta_l_prime) - Lambda_func(V_eta_l_prime_minus_1)));
                    
                    Lambda_sum         = Lambda_func(V_eta_l_prime_minus_1) + R;
                    v_next_transformed = Lambda_inv_func(Lambda_sum);
                    v_next             = v_next_transformed / exp(z_l);
                    T_next             = V_inv_func(v_next, last_event, rho);
                end
            end

            if isinf(T_next) || T_next > eta(end) || T_next <= T_j
                break;
            end
            
            fail_times = [fail_times, T_next];
        end

        all_fail_times{i, rr} = fail_times;
        all_covariates{i, rr} = covariates;
    end
end

% --- Diagnostic check ---
% Check monotonicity for first subject
for rr = 1:length(rho_set)
    fail_times = all_fail_times{1, rr};
    if length(fail_times) > 1 && any(diff(fail_times) <= 0)
        fprintf('Warning: Non-monotonic event times detected for rho=%.1f!\n', rho_set(rr));
    else
        fprintf('Event times are strictly increasing for rho=%.1f\n', rho_set(rr));
    end
    
    % Display first subject results
    fprintf('Subject 1, rho=%.1f: %s\n', rho_set(rr), mat2str(fail_times, 4));
end

% --- Summary Statistics ---
fprintf('\n--- Simulation Summary ---\n');
for rr = 1:length(rho_set)
    event_counts = cellfun(@(x) length(x)-1, all_fail_times(:, rr)); % -1 because first element is T_0=0
    fprintf('rho=%.1f: Mean events per subject = %.2f (SD = %.2f)\n', ...
        rho_set(rr), mean(event_counts), std(event_counts));
end

% -------------------------------------------------------------------------
% SIMULATION OF RECURRENT EVENTS UNDER sARA1-TDC MODEL WITH CENSORING
% + counting-process plot
% -------------------------------------------------------------------------

rng(12); % Reproducibility for censored simulation

% Run simulation with censoring
eventsTable = simulate_recurrent_events_with_censoring();

% Inspect first rows
disp('Head of eventsTable:');
disp(eventsTable(1:min(12, height(eventsTable)), :));

% Basic sanity checks
assert(all(ismember(eventsTable.event, [0 1])), 'event must be 0 or 1');
assert(all(eventsTable.time <= eventsTable.tau + 1e-12), 'No times should exceed tau');

% Per-subject checks
uids = unique(eventsTable.id);
nsub = numel(uids);

event_counts = zeros(nsub,1);
max_gap_err  = 0;

for k = 1:nsub
    id = uids(k);
    idx = eventsTable.id == id;
    T   = eventsTable.time(idx);
    G   = eventsTable.gap_time(idx);
    E   = eventsTable.event(idx);
    Tau = eventsTable.tau(idx);

    % times strictly increasing
    assert(all(diff(T) > 0), sprintf('Non-increasing times for id=%d', id));

    % tau consistent within subject
    assert(isscalar(unique(Tau)), sprintf('Multiple tau values for id=%d', id));
    tau_i = Tau(1);

    % last record should be censoring
    assert(E(end) == 0, sprintf('Last row must be censoring for id=%d', id));
    assert(abs(T(end) - tau_i) < 1e-10, sprintf('Last time must equal tau for id=%d', id));

    % gap_time consistency: G(1) == T(1), for m>=2, G(m) == T(m) - T(m-1)
    if ~isempty(T)
        max_gap_err = max(max_gap_err, abs(G(1) - T(1)));
        if numel(T) >= 2
            max_gap_err = max(max_gap_err, max(abs(diff(T) - G(2:end))));
        end
    end

    % count events
    event_counts(k) = sum(E == 1);
end

fprintf('All sanity checks passed.\n');
fprintf('Max absolute gap_time consistency error: %.3e\n', max_gap_err);

% Censoring percentages
total_rows            = height(eventsTable);
censored_rows         = sum(eventsTable.event == 0);
pct_censored_rows     = 100 * censored_rows / total_rows;

subjects_with_no_ev   = sum(event_counts == 0);         % subjects with zero events (only censoring)
pct_subjects_no_ev    = 100 * subjects_with_no_ev / nsub;

subjects_with_events  = sum(event_counts > 0);
pct_subjects_with_ev  = 100 * subjects_with_events / nsub;

% Summary statistics
fprintf('\n--- Summary (with censoring) ---\n');
fprintf('Subjects: %d\n', nsub);
fprintf('Total rows: %d\n', total_rows);
fprintf('Total events: %d\n', sum(eventsTable.event == 1));
fprintf('Mean events per subject: %.2f (SD = %.2f)\n', mean(event_counts), std(event_counts));
fprintf('Censored rows: %d (%.2f%% of rows)\n', censored_rows, pct_censored_rows);
fprintf('Subjects with no events: %d (%.2f%% of subjects)\n', subjects_with_no_ev, pct_subjects_no_ev);
fprintf('Subjects with at least one event: %d (%.2f%% of subjects)\n', subjects_with_events, pct_subjects_with_ev);

% Show a random subject timeline
id_show = uids(randi(nsub));
idx = eventsTable.id == id_show;
T   = eventsTable.time(idx);
E   = eventsTable.event(idx);
tau_i = eventsTable.tau(find(idx,1));

fprintf('\nExample subject id=%d: events=%d, tau=%.3f\n', id_show, sum(E==1), tau_i);
disp(eventsTable(idx, :));

% Plot counting process N_i(t): step function that jumps by 1 at each event
figure('Name', sprintf('Subject %d counting process', id_show), 'Color','w', 'Position',[100 100 1100 650]);
hold on;

T_events = T(E==1);
nEv = numel(T_events);

if nEv > 0
    % Build right-continuous step function from 0 to tau_i
    X = [0; T_events(:); tau_i];
    Y = [0; (1:nEv)'; nEv];
    stairs(X, Y, 'b-', 'LineWidth', 1.8, 'DisplayName', 'N_i(t)');
    % Mark event jump points
    plot(T_events(:), (1:nEv)', 'bo', 'MarkerFaceColor', 'b', 'HandleVisibility','off');
else
    % No events: flat zero line until censoring
    stairs([0; tau_i], [0; 0], 'b-', 'LineWidth', 1.8, 'DisplayName', 'N_i(t)');
end

% Censoring time
xline(tau_i, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Censoring \tau');

ylim([0, max(1, nEv)+0.5]);
xlim([0, tau_i]);
xlabel('Time');
ylabel('Cumulative event count N_i(t)');
title(sprintf('Subject %d: Counting process and censoring', id_show), 'FontSize', 16, 'FontWeight', 'bold');
legend('Location','northwest');
grid on; box on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local function(s) must appear after all script code in MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eventsTable = simulate_recurrent_events_with_censoring()
% SIMULATE_RECURRENT_EVENTS_WITH_CENSORING
% Simulation of recurrent events under sARA1 with the hazard:
%   lambda_i(t;theta) = lambda0(V_i(t)*exp(theta'Z_i(t))) * exp(theta'Z_i(t))
% Piecewise-constant covariates across calendar-time intervals.
% Adds random right-censoring tau ~ Uniform(0, tau_max) per subject and
% returns a table with columns: id, time, gap_time, event, covariate, tau.

% User parameters
n           = 100;                  % Number of subjects
max_events  = 10;                   % Max number of events per subject
eta         = [0, 1, 2, 3, 4, 5];   % Interval boundaries (must cover tau_max)
L           = numel(eta) - 1;       % Number of intervals
theta       = 0.5;                  % Scalar parameter for simplicity
rho         = 0.7;                  % sARA1 efficiency (0 perfect repair, 1 minimal repair)
alpha       = 1.5;                  % Weibull baseline cumulative hazard: Lambda(v)=alpha*v^gamma
gamma       = 1.5;
sigma2      = 0.1;                  % Variance for interval covariates
tau_max     = 3.5;                  % Maximum value for tau

% Censoring times
tau = unifrnd(0, tau_max, [n, 1]);

% Baseline cumulative hazard and its inverse (Weibull/PLP)
Lambda0      = @(v) alpha * v.^gamma;
Lambda0_inv  = @(x) (x ./ alpha).^(1 ./ gamma);

% sARA1 virtual age process and inverse
% After an event at T_prev, V(t) = t - rho*T_prev (unit slope between events)
V      = @(t, T_prev) t - rho * T_prev;
V_inv  = @(v, T_prev) v + rho * T_prev;

% Storage for output table
ids          = [];
times        = [];
gap_times    = [];
events       = [];   % 1 = event, 0 = censored
covariates   = [];
taus_col     = [];

% Simulation loop over subjects
for i = 1:n
    % Draw piecewise-constant covariates per interval for this subject
    % Z_i(l) applies on [eta(l), eta(l+1))
    Zi = sqrt(sigma2) * randn(L, 1);  % scalar covariate per interval
    T_prev = 0;                       % T_0
    last_event_time = 0;              % For sARA1 V(t)
    subject_done = false;

    for j = 1:max_events
        T_j = T_prev;

        % If already beyond censoring, stop (shouldn't happen due to checks)
        if T_j >= tau(i)
            subject_done = true;
            break;
        end

        % Step 1: Locate current interval l_j such that T_j in [eta(l_j), eta(l_j+1))
        l_j = find(T_j >= eta(1:end-1) & T_j < eta(2:end), 1, 'last');
        if isempty(l_j)
            % If T_j is out of bounds relative to eta, censor at tau(i)
            ids        = [ids; i];
            times      = [times; tau(i)];
            gap_times  = [gap_times; tau(i) - T_prev];
            events     = [events; 0];
            % Covariate at censoring interval
            l_tau = find(tau(i) >= eta(1:end-1) & tau(i) < eta(2:end), 1, 'last');
            if isempty(l_tau), l_tau = L; end
            covariates = [covariates; Zi(l_tau)];
            taus_col   = [taus_col; tau(i)];
            subject_done = true;
            break;
        end

        % Step 2: Target quantile
        Q = -log(1 - rand);  % Exponential(1)

        % Step 3: Check if next event stays in current interval
        z_lj = theta * Zi(l_j);             % linear predictor for interval l_j
        V_Tj = V(T_j, last_event_time);
        V_eta_lj = V(eta(l_j+1), last_event_time);

        % Transform virtual ages by exp(z) according to the hazard specification
        V_tilde_Tj    = V_Tj * exp(z_lj);
        V_tilde_etaLj = V_eta_lj * exp(z_lj);

        % Available cumulative hazard within the current interval
        H_current = Lambda0(V_tilde_etaLj) - Lambda0(V_tilde_Tj);

        if Q <= H_current
            % Event occurs within current interval
            v_tilde_next = Lambda0_inv(Lambda0(V_tilde_Tj) + Q);
            v_next       = v_tilde_next / exp(z_lj);
            T_next       = V_inv(v_next, last_event_time);
            z_used       = Zi(l_j);

        else
            % Event crosses intervals: accumulate hazard
            I = H_current;
            l = l_j + 1;
            l_prime = NaN;
            z_l = NaN;

            while I < Q && l <= L
                z_l       = theta * Zi(l);
                V_eta_lm1 = V(eta(l),   last_event_time) * exp(z_l);
                V_eta_l   = V(eta(l+1), last_event_time) * exp(z_l);
                I = I + (Lambda0(V_eta_l) - Lambda0(V_eta_lm1));
                if I >= Q
                    l_prime = l;
                    break;
                end
                l = l + 1;
            end

            if isnan(l_prime)
                % Could not place the event before exceeding the last interval.
                % Censor at tau(i).
                ids        = [ids; i];
                times      = [times; tau(i)];
                gap_times  = [gap_times; tau(i) - T_prev];
                events     = [events; 0];
                l_tau = find(tau(i) >= eta(1:end-1) & tau(i) < eta(2:end), 1, 'last');
                if isempty(l_tau), l_tau = L; end
                covariates = [covariates; Zi(l_tau)];
                taus_col   = [taus_col; tau(i)];
                subject_done = true;
                break;
            else
                % Compute remaining hazard in interval l'
                V_eta_lp1m1 = V(eta(l_prime),   last_event_time) * exp(z_l);
                V_eta_lp1   = V(eta(l_prime+1), last_event_time) * exp(z_l);
                R = Q - (I - (Lambda0(V_eta_lp1) - Lambda0(V_eta_lp1m1)));

                % Final inversion within l'
                v_tilde_next = Lambda0_inv(Lambda0(V_eta_lp1m1) + R);
                v_next       = v_tilde_next / exp(z_l);
                T_next       = V_inv(v_next, last_event_time);
                z_used       = Zi(l_prime);
            end
        end

        % Apply censoring: if the next event exceeds tau(i), record censoring and stop
        if ~(isfinite(T_next)) || T_next <= T_prev || T_next > tau(i)
            % Censor at tau(i)
            ids        = [ids; i];
            times      = [times; tau(i)];
            gap_times  = [gap_times; tau(i) - T_prev];
            events     = [events; 0];
            % Covariate at censoring time interval
            l_tau = find(tau(i) >= eta(1:end-1) & tau(i) < eta(2:end), 1, 'last');
            if isempty(l_tau), l_tau = L; end
            covariates = [covariates; Zi(l_tau)];
            taus_col   = [taus_col; tau(i)];
            subject_done = true;
            break;
        end

        % Record the event row
        ids        = [ids; i];
        times      = [times; T_next];
        gap_times  = [gap_times; T_next - T_prev];
        events     = [events; 1];
        covariates = [covariates; z_used];
        taus_col   = [taus_col; tau(i)];

        % Update for next iteration
        last_event_time = T_next;
        T_prev          = T_next;
    end

    % If the subject had no events and was not censored inside the loop, ensure censoring is recorded
    if ~subject_done
        ids        = [ids; i];
        times      = [times; tau(i)];
        gap_times  = [gap_times; tau(i) - T_prev];
        events     = [events; 0];
        l_tau = find(tau(i) >= eta(1:end-1) & tau(i) < eta(2:end), 1, 'last');
        if isempty(l_tau), l_tau = L; end
        covariates = [covariates; Zi(l_tau)];
        taus_col   = [taus_col; tau(i)];
    end
end

% Create output table
eventsTable = table( ...
    ids, ...
    times, ...
    gap_times, ...
    events, ...
    covariates, ...
    taus_col, ...
    'VariableNames', {'id', 'time', 'gap_time', 'event', 'covariate', 'tau'} ...
);
end
