/* ================= DATA IMPORT ================= */
proc import
    out = dat
    datafile = "C:\Users\flash\Documents\2026\Masters\MVA\2026 Exam\randmixA.xlsx"
    dbms = xlsx
    replace;
    getnames = yes;
run;

/* ================= PLOT HISTOGRAM FOR X ================= */
title "Distribution of Variable X";
proc sgplot data=dat;
    histogram X / fillattrs=(color=cx4F81BD) transparency=0.3;
    density X / type=normal lineattrs=(color=cxC0504D thickness=2 pattern=dash) legendlabel="Normal Curve";

    xaxis label="X Values";
    yaxis label="Percent" grid;
    keylegend / location=inside position=topright across=1;
run;
title;



/* ================= KDE & CROSS-VALIDATION ================= */
proc iml;
    use dat;
    read all var{X} into x;
    close dat;

    n = nrow(x);
    std_x = std(x);

    print "Sample Size:" n;
    print "Standard Deviation of X:" std_x;

    /* 1. Define the Gaussian KDE CDF function */
    start kde_cdf(eval_points, train_data, h);
        n_train = nrow(train_data);
        n_eval = nrow(eval_points);
        cdf_vals = j(n_eval, 1, 0);
        do i = 1 to n_eval;
            u = (eval_points[i] - train_data) / h;
            cdf_vals[i] = sum(cdf("Normal", u, 0, 1)) / n_train;
        end;
        return(cdf_vals);
    finish;

    /* 2. Define the Empirical CDF (ECDF) function */
    start ecdf(eval_points, val_data);
        n_val = nrow(val_data);
        n_eval = nrow(eval_points);
        ecdf_vals = j(n_eval, 1, 0);
        do i = 1 to n_eval;
            ecdf_vals[i] = sum(val_data <= eval_points[i]) / n_val;
        end;
        return(ecdf_vals);
    finish;

    /* 3. Setup K-Fold Cross Validation */
    k = 5;
    call randseed(12345);

    /* FIX 1: Ensure idx is strictly a column vector */
    idx = sample(1:n, n, "NoRep");
    if nrow(idx) = 1 then idx = T(idx);

    fold_size = floor(n / k);

    /* 4. Define bandwidth grid and include standard deviation */
    h_base = T(do(0.1, 2.0, 0.05));
    h_grid = h_base // std_x;
    call sort(h_grid, 1);

    num_h = nrow(h_grid);
    avg_ks = j(num_h, 1, 0);

    /* 5. Execute the Cross-Validation Loop */
    do i = 1 to num_h;
        h_test = h_grid[i];
        ks_folds = j(k, 1, 0);

        do j = 1 to k;
            /* Isolate validation fold indices */
            val_start = (j-1)*fold_size + 1;
            if j = k then val_end = n;
            else val_end = j*fold_size;

            val_idx = idx[val_start:val_end];

            /* FIX 2: Training fold using vertical concatenation (//) */
            if val_start = 1 then
                train_idx = idx[(val_end+1):n];
            else if val_end = n then
                train_idx = idx[1:(val_start-1)];
            else
                train_idx = idx[1:(val_start-1)] // idx[(val_end+1):n];

            train_data = x[train_idx];
            val_data = x[val_idx];

            /* Evaluate CDFs strictly on the out-of-sample validation data */
            hat_F = kde_cdf(val_data, train_data, h_test);
            F_n   = ecdf(val_data, val_data);

            /* KS Statistic is the maximum absolute difference */
            ks_folds[j] = max(abs(hat_F - F_n));
        end;

        /* Average the KS statistic across all folds for this specific h */
        avg_ks[i] = mean(ks_folds);
    end;

    /* 6. Identify the optimal bandwidth */
    min_ks_idx = loc(avg_ks = min(avg_ks));
    best_h = h_grid[min_ks_idx[1]];

    print "Cross-Validation Complete";
    print "Tested std(x) as bandwidth alongside the grid." ;
    print "Optimal Bandwidth (h) based on KS Statistic:" best_h;

    /* 7. Export the results to a SAS dataset for plotting */
    ks_results = h_grid || avg_ks;
    create cv_results from ks_results[colname={"h" "KS_Statistic"}];
    append from ks_results;
    close cv_results;
quit;


/* ================= PLOT KS STATISTIC VS BANDWIDTH ================= */
title "Out-of-Sample KS Statistic vs. Bandwidth (h)";
proc sgplot data=cv_results;

    /* This tells SGPLOT to ignore anything above 2.0 just for this graph */
    where h <= 2;

    series x=h y=KS_Statistic /
           markers markerattrs=(symbol=circlefilled color=cx4F81BD)
           lineattrs=(color=cx4F81BD thickness=2);

    xaxis label="Bandwidth (h)";
    yaxis label="Average KS Statistic" grid;
run;
title;





proc iml;
    use dat;
    read all var{X} into x;
    close dat;
    n = nrow(x);

    /* FORCE h_list to be a row vector of 5 columns */
    h_list = {0.55 1 2 5 9};

    min_x = min(x) - 3*max(h_list);
    max_x = max(x) + 3*max(h_list);
    grid = T( do(min_x, max_x, (max_x-min_x)/300) );

    /* Pre-allocate total results matrix for speed and reliability */
    /* 1505 rows (301 grid points * 5 bandwidths) by 3 columns */
    all_kde_data = j(1505, 3, 0);

    /* Loop through columns (h_list is 1x5) */
    do k = 1 to ncol(h_list);
        h_val = h_list[k];

        /* Calculate density for this h */
        density = j(nrow(grid), 1, 0);
        do i = 1 to nrow(grid);
            u = (grid[i] - x) / h_val;
            k_u = exp(-0.5 * u##2) / sqrt(2*constant('PI'));
            density[i] = sum(k_u) / (n * h_val);
        end;

        h_col = j(nrow(grid), 1, h_val);

        /* Insert the chunk into the pre-allocated matrix */
        start_row = (k-1)*301 + 1;
        end_row = k*301;
        all_kde_data[start_row:end_row, ] = (grid || density || h_col);
    end;

    create kde_multi from all_kde_data[colname={"X_grid" "Density" "H"}];
    append from all_kde_data;
    close kde_multi;
quit;


proc template;
    define statgraph multi_kde_plot;
        begingraph;
            entrytitle "KDE Comparison across different Bandwidths (h)";

            /* Add a legend for the KDE group */
            legenditem type=line name="kde" /
                label="KDE Bandwidths"
                lineattrs=(thickness=2);

            layout overlay /
                yaxisopts=(label="Percent" griddisplay=on)
                y2axisopts=(label="Density");

                /* 1. Histogram */
                histogram X / scale=percent fillattrs=(color=ligr)
                              datatransparency=0.6 name="hist"
                              legendlabel="Raw Data Histogram";

                /* 2. Multiple series lines grouped by H */
                seriesplot x=X_grid y=Density / group=H
                                                yaxis=y2
                                                lineattrs=(thickness=2)
                                                name="kde_lines";

                /* This legend calls the specific names defined above */
                discretelegend "hist" "kde_lines" /
                    location=inside
                    autoalign=(topright)
                    title="Legend";
            endlayout;
        endgraph;
    end;
run;

/* Re-run the render */
proc sgrender data=combined_multi template=multi_kde_plot;
run;






/* ================= NUMERICAL PROOF OF DENSITY FUNCTION ================= */
proc iml;
    use dat;
    read all var{X} into x;
    close dat;
    n = nrow(x);

    best_h = 0.55;

    /* Define a very fine grid for accurate numerical integration */
    min_x = min(x) - 4*best_h;
    max_x = max(x) + 4*best_h;
    ngrid = 5000;

    /* Calculate the exact width (dx) of each tiny slice of the grid */
    step_size = (max_x - min_x) / ngrid;
    grid = T( do(min_x, max_x, step_size) );

    /* Calculate the density at each grid point */
    density = j(nrow(grid), 1, 0);
    do i = 1 to nrow(grid);
        u = (grid[i] - x) / best_h;
        k_u = exp(-0.5 * u##2) / sqrt(2*constant('PI'));
        density[i] = sum(k_u) / (n * best_h);
    end;

    /* --- THE PROOF --- */

    /* 1. Prove Non-Negativity */
    min_density = min(density);
    is_positive = (min_density >= 0);

    /* 2. Prove it Integrates to 1 (Riemann Sum: Area = sum of (height * width)) */
    total_area = sum(density * step_size);

    print "========== NUMERICAL KDE PROOF (h=0.55) ==========";
    print "Test 1: Is the minimum density value >= 0?";
    print min_density is_positive[format=best1. label="Pass (1=Yes)"];

    print "Test 2: Does the total area under the curve integrate to 1?";
    print total_area[format=8.5 label="Total Area"];
quit;





/* ================= DISCRETE k-NN MODE-SEEKING CLUSTERING ================= */
proc iml;
    use dat;
    read all var{X} into x;
    close dat;

    n = nrow(x);
    h = 0.55; /* The optimal bandwidth from Q2.1 */
    k_nn = 50; /* Set the number of nearest neighbors (adjust if needed) */

    print "Step 1: Calculating KDE density for all data points...";
    f_x = j(n, 1, 0);
    do i = 1 to n;
        u = (x[i] - x) / h;
        k_u = exp(-0.5 * u##2) / sqrt(2*constant('PI'));
        f_x[i] = sum(k_u) / (n * h);
    end;

    print "Step 2: Executing k-NN Hill-Climbing...";
    converged_idx = j(n, 1, 0); /* Store the index of the final mode */

    do i = 1 to n;
        curr_idx = i;
        moved = 1;

        do while(moved = 1);
            /* Find spatial distances to all points */
            dist = abs(x - x[curr_idx]);

            /* Identify the k nearest neighbors using rank */
            r = rank(dist);
            nn_indices = loc(r <= k_nn);

            /* Look at the densities of those specific neighbors */
            nn_densities = f_x[nn_indices];
            max_dens_local = max(nn_densities);

            /* Find which neighbor has that maximum density */
            best_local_idx = nn_indices[loc(nn_densities = max_dens_local)][1];

            /* If the best neighbor is myself, I am at a local peak */
            if best_local_idx = curr_idx then do;
                moved = 0;
                converged_idx[i] = curr_idx;
            end;
            else do;
                curr_idx = best_local_idx; /* Move to the denser neighbor */
            end;
        end;
    end;

    /* Extract the actual X values of the final modes */
    modes = x[converged_idx];
    rounded_modes = round(modes, 0.1);
    unique_modes = unique(rounded_modes);

    /* Assign integer cluster labels based on the converged modes */
    cluster_idx = j(n, 1, 0);
    do i = 1 to n;
        d = abs(unique_modes - modes[i]);
        cluster_idx[i] = loc(d = min(d))[1];
    end;

    print "Clustering Complete.";
    print "Identified Modes:" unique_modes;

    /* Export for Plotting */
    out_data = x || cluster_idx || modes;
    create knn_clusters from out_data[colname={"X" "Cluster" "Converged_Mode"}];
    append from out_data;
    close knn_clusters;
quit;

/* ================= REPRESENTATION OF CLUSTERING RESULTS ================= */
title "Mode-Seeking Clustering Results ( h=0.55)";
proc sgplot data=knn_clusters;
    histogram X / group=Cluster fillattrs=(transparency=0.4) name="hist" binwidth=2;
    xaxis label="X Values";
    yaxis label="Percent" grid;
    discretelegend "hist" / title="Assigned Cluster";
run;
title;
