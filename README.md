# sARA1 Recurrent Events with Time‑Varying Covariates

Simulate and visualize recurrent events under a virtual‑age (sARA1) model with right censoring and time‑varying covariates.

## What this repo provides
- Simulation of recurrent events under sARA1 with a Weibull/PLP baseline.
- Support for time‑varying covariates across calendar‑time intervals.
- Random right‑censoring per subject (uniform on [0, τ_max]).
- Construction of a counting‑process plot N_i(t) that jumps by 1 at each event time.
- Sanity checks and summary statistics, including percent censored.

## How to run
1. Open MATLAB.
2. Run the script:
   ```matlab
   sARA_TDC_Data_Sim
   ```
3. The script will:
   - Run a finite‑sample simulation without censoring across multiple ρ values.
   - Simulate a dataset with right censoring and piecewise time‑varying covariates.
   - Print validation checks and summary stats (including censoring percentages).
   - Display a counting‑process step plot for a randomly chosen subject.

## Save the example figure (optional)
To save the counting‑process plot shown in the script to your repo (e.g., docs/images), add and run:
```matlab
outDir = fullfile('docs','images'); if ~exist(outDir,'dir'), mkdir(outDir); end
exportgraphics(gcf, fullfile(outDir,'recurrent_event_example.png'), 'Resolution', 200);
```

## Notes
- The script defines a local function `simulate_recurrent_events_with_censoring` at the end of the file (MATLAB requires local functions to appear after script code).
- Parameters like `n`, `tau_max`, baseline hazard, and ρ are set near the top of each section and can be adjusted as needed.

## Result

<img width="1096" height="657" alt="Counting Process" src="https://github.com/user-attachments/assets/c57b3e68-c859-42db-abb1-0f76f633d784" />
