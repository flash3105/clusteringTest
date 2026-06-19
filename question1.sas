proc import
        out = dat
        datafile =  "C:\Users\flash\Documents\2026\Masters\MVA\2026 Exam\randmixA.xlsx"
        dbms = xlsx
        replace;
        getnames = yes;

run;


/* ================= PLOT HISTOGRAM FOR X ================= */
title "Distribution of Variable X";
proc sgplot data=dat;
    /* Main histogram */
    histogram X / fillattrs=(color=cx4F81BD) transparency=0.3;

    /* Overlay normal density curve */
    density X / type=normal lineattrs=(color=cxC0504D thickness=2 pattern=dash) legendlabel="Normal Curve";

    xaxis label="X Values";
    yaxis label="Percent" grid;
    keylegend / location=inside position=topright across=1;
run;
title; /* clear title */


/* ================= PLOT HISTOGRAM FOR Y ================= */
title "Distribution of Variable Y";
proc sgplot data=dat;
    /* Main histogram */
    histogram Y / fillattrs=(color=cx9BBB59) transparency=0.3;

    /* Overlay normal density curve */
    density Y / type=normal lineattrs=(color=cxC0504D thickness=2 pattern=dash) legendlabel="Normal Curve";

    xaxis label="Y Values";
    yaxis label="Percent" grid;
    keylegend / location=inside position=topright across=1;
run;
title; /* clear title */




proc iml;
    use dat;
    read all var{y} into x; /* Changed from var{y} based on your new code */
    n = nrow(x);
    print "Sample Size:" n;

    /* ================= INITIAL PARAMETERS ================= */
    /* Use spaces to ensure 1x3 row vectors */
    mu = {100 104 111};
    si = {8 1 4.5};
    pi = {0.3 0.3 0.4};

    stop = 0;

    /* Safely initialize an empty matrix for 10 columns */
    free results;

    do i = 1 to 200 while (stop=0);

        /* ===== E-step: Calculate responsibilities (gamma) ===== */
        gamma = pi[1]*PDF("Normal",x,mu[1],si[1]) ||
                pi[2]*PDF("Normal",x,mu[2],si[2]) ||
                pi[3]*PDF("Normal",x,mu[3],si[3]);

        gamma = gamma / gamma[,+];

        /* ===== M-step: Update parameters ===== */
        nk = gamma[+,];
        npi = nk / nk[+];
        nmu = (x#gamma)[+,] / nk;

        /* Split variance calculation to avoid dimensionality crashes */
        nsi2_1 = (gamma[,1]#(x-nmu[1]))`*(x-nmu[1])/nk[1];
        nsi2_2 = (gamma[,2]#(x-nmu[2]))`*(x-nmu[2])/nk[2];
        nsi2_3 = (gamma[,3]#(x-nmu[3]))`*(x-nmu[3])/nk[3];
        nsi = sqrt(nsi2_1 || nsi2_2 || nsi2_3);

        /* ===== Check convergence ===== */
        diff = max(abs(mu - nmu));
        if diff < 0.005 then stop=1;

        /* ===== Update and store ===== */
        mu = nmu;
        si = nsi;
        pi = npi;

        /* Dynamically build the results matrix (1 + 3 + 3 + 3 = 10 cols) */
        if i = 1 then results = (i || npi || nmu || nsi);
        else results = results // (i || npi || nmu || nsi);
    end;

    /* ===== Print Iteration Results ===== */
    nm = {"it" "pi1" "pi2" "pi3" "mu1" "mu2" "mu3" "si1" "si2" "si3"};
    print "EM Convergence Results";
    print results[colname=nm];

    /* ===== Estimated log-likelihood for 3-component GMM ===== */
    ll_gmm = sum(log(
        pi[1]*PDF("Normal",x,mu[1],si[1]) +
        pi[2]*PDF("Normal",x,mu[2],si[2]) +
        pi[3]*PDF("Normal",x,mu[3],si[3])
    ));

    /* ===== Null model: fit single normal ===== */
    mean0 = mean(x);
    std0 = std(x);
    ll_null = sum(log(PDF("Normal",x,mean0,std0)));

    /* ===== Model Metrics (AIC, BIC, LRT) ===== */
    k = 8; /* 3 means + 3 std devs + 2 independent mixing proportions */
    aic = 2*k - 2*ll_gmm;
    bic = log(n)*k - 2*ll_gmm;
    lrt = -2*(ll_null - ll_gmm);

    metrics = ll_gmm || ll_null || aic || bic || lrt;
    names = {"LL_GMM" "LL_Null" "AIC" "BIC" "LRT_Stat"};
    print "Model Fit Statistics";
    print metrics[colname=names];

    /* ===== Create estimated density plot data ===== */
    xplot = T(do(min(x), max(x), (max(x)-min(x))/200));
    gmmd = pi[1]*PDF("Normal",xplot,mu[1],si[1]) +
           pi[2]*PDF("Normal",xplot,mu[2],si[2]) +
           pi[3]*PDF("Normal",xplot,mu[3],si[3]);

    plotdens = xplot || gmmd;
    nm1 = {"y","fy"};
    create densplot from plotdens[colname=nm1];
    append from plotdens;

    create plotd from results[colname=nm];
    append from results;

    /* ===== Create parameter dataset for plotting ===== */
    params = mu || si || pi;
    names = {"mu1" "mu2" "mu3" "si1" "si2" "si3" "pi1" "pi2" "pi3"};
    create gmm_params from params[colname=names];
    append from params;
quit;


/* ================= VISUALIZATIONS ================= */

proc sgplot data=plotd;
    title "Evolution of Mixing Proportions (3 Components)";
    series x=it y=pi1 / legendlabel="pi1" lineattrs=(thickness=2);
    series x=it y=pi2 / legendlabel="pi2" lineattrs=(thickness=2);
    series x=it y=pi3 / legendlabel="pi3" lineattrs=(thickness=2);
    yaxis grid;
run;

proc sgplot data=plotd;
    title "Evolution of Means (3 Components)";
    series x=it y=mu1 / legendlabel="mu1" lineattrs=(thickness=2);
    series x=it y=mu2 / legendlabel="mu2" lineattrs=(thickness=2);
    series x=it y=mu3 / legendlabel="mu3" lineattrs=(thickness=2);
    yaxis grid;
run;

proc sgplot data=plotd;
    title "Evolution of Standard Deviations (3 Components)";
    series x=it y=si1 / legendlabel="si1" lineattrs=(thickness=2);
    series x=it y=si2 / legendlabel="si2" lineattrs=(thickness=2);
    series x=it y=si3 / legendlabel="si3" lineattrs=(thickness=2);
    yaxis grid;
run;

proc sgplot data=densplot;
    title "Estimated Density using 3-Component GMM";
    series x=y y=fy / lineattrs=(color=blue thickness=2);
    xaxis grid;
    yaxis grid;
run;
title; /* clear title */



proc iml;
    use dat;
    /* Read your variable (assuming it is Y based on your snippet) */
    read all var{y} into x;
    n = nrow(x);

    /* ================= CONVERGED PARAMETERS ================= */
    /* Inserted from your final EM iteration */
    pi = {0.2014942  0.3089843  0.4895215};
    mu = {110.19112  100.08072  120.58842};
    si = {2.8090352  2.1794235  7.8177596};

    /* ================= CALCULATE RESPONSIBILITIES ================= */
    /* Calculate the unnormalized weighted density for each component */
    gamma = pi[1]*PDF("Normal", x, mu[1], si[1]) ||
            pi[2]*PDF("Normal", x, mu[2], si[2]) ||
            pi[3]*PDF("Normal", x, mu[3], si[3]);

    /* Normalize to get true probabilities (0 to 1) */
    gamma = gamma / gamma[,+];

    /* ================= HARD CLUSTERING ================= */
    cluster_id = J(n, 1, .);

    do i = 1 to n;
        /* Extract the 3 probabilities for the current observation */
        row_probs = gamma[i,];
        max_prob = max(row_probs);

        /* Assign observation to the cluster with the highest probability */
        if row_probs[1] = max_prob then cluster_id[i] = 1;
        else if row_probs[2] = max_prob then cluster_id[i] = 2;
        else cluster_id[i] = 3;
    end;

    /* ================= EXPORT DATA ================= */
    /* Combine original data, the assigned cluster, and the probabilities */
    out_data = x || cluster_id || gamma;
    col_names = {"Y_Value" "Cluster" "Prob_C1" "Prob_C2" "Prob_C3"};

    create clustered_data from out_data[colname=col_names];
    append from out_data;
    close clustered_data;

    print "Hard Clustering Complete. Data saved to 'clustered_data'.";
quit;


/* ================= VISUALIZE CLUSTERS ================= */

/* 1. Grouped Histogram to see how the distributions overlap */
proc sgplot data=clustered_data;
    title "Hard Clustering Results: 3-Component GMM";
    histogram Y_Value / group=Cluster transparency=0.4 binwidth=2;
    xaxis label="Observed Value";
    yaxis label="Density" grid;
    keylegend / location=inside position=topright;
run;

/* 2. Scatter plot to see exact assignments */
proc sgplot data=clustered_data;
    title "Observation Assignments by Cluster";
    scatter x=Y_Value y=Cluster / group=Cluster markerattrs=(symbol=CircleFilled size=10);
    xaxis label="Observed Value";
    yaxis label="Assigned Cluster (1, 2, or 3)" integer grid;
run;
title; /* clear title */






/* ================= ALL FOR x  ================= */


proc iml;
    use dat;
    read all var{X} into x; /
    n = nrow(x);
    print "Sample Size:" n;

    /* ================= INITIAL PARAMETERS ================= */
    /* Use spaces to ensure 1x3 row vectors */
    mu = {100 104 111};
    si = {8 1 4.5};
    pi = {0.3 0.3 0.4};

    stop = 0;

    /* Safely initialize an empty matrix for 10 columns */
    free results;

    do i = 1 to 200 while (stop=0);

        /* ===== E-step: Calculate responsibilities (gamma) ===== */
        gamma = pi[1]*PDF("Normal",x,mu[1],si[1]) ||
                pi[2]*PDF("Normal",x,mu[2],si[2]) ||
                pi[3]*PDF("Normal",x,mu[3],si[3]);

        gamma = gamma / gamma[,+];

        /* ===== M-step: Update parameters ===== */
        nk = gamma[+,];
        npi = nk / nk[+];
        nmu = (x#gamma)[+,] / nk;

        /* Split variance calculation to avoid dimensionality crashes */
        nsi2_1 = (gamma[,1]#(x-nmu[1]))`*(x-nmu[1])/nk[1];
        nsi2_2 = (gamma[,2]#(x-nmu[2]))`*(x-nmu[2])/nk[2];
        nsi2_3 = (gamma[,3]#(x-nmu[3]))`*(x-nmu[3])/nk[3];
        nsi = sqrt(nsi2_1 || nsi2_2 || nsi2_3);

        /* ===== Check convergence ===== */
        diff = max(abs(mu - nmu));
        if diff < 0.005 then stop=1;

        /* ===== Update and store ===== */
        mu = nmu;
        si = nsi;
        pi = npi;

        /* Dynamically build the results matrix (1 + 3 + 3 + 3 = 10 cols) */
        if i = 1 then results = (i || npi || nmu || nsi);
        else results = results // (i || npi || nmu || nsi);
    end;

    /* ===== Print Iteration Results ===== */
    nm = {"it" "pi1" "pi2" "pi3" "mu1" "mu2" "mu3" "si1" "si2" "si3"};
    print "EM Convergence Results";
    print results[colname=nm];

    /* ===== Estimated log-likelihood for 3-component GMM ===== */
    ll_gmm = sum(log(
        pi[1]*PDF("Normal",x,mu[1],si[1]) +
        pi[2]*PDF("Normal",x,mu[2],si[2]) +
        pi[3]*PDF("Normal",x,mu[3],si[3])
    ));

    /* ===== Null model: fit single normal ===== */
    mean0 = mean(x);
    std0 = std(x);
    ll_null = sum(log(PDF("Normal",x,mean0,std0)));

    /* ===== Model Metrics (AIC, BIC, LRT) ===== */
    k = 8; /* 3 means + 3 std devs + 2 independent mixing proportions */
    aic = 2*k - 2*ll_gmm;
    bic = log(n)*k - 2*ll_gmm;
    lrt = -2*(ll_null - ll_gmm);

    metrics = ll_gmm || ll_null || aic || bic || lrt;
    names = {"LL_GMM" "LL_Null" "AIC" "BIC" "LRT_Stat"};
    print "Model Fit Statistics";
    print metrics[colname=names];

    /* ===== Create estimated density plot data ===== */
    xplot = T(do(min(x), max(x), (max(x)-min(x))/200));
    gmmd = pi[1]*PDF("Normal",xplot,mu[1],si[1]) +
           pi[2]*PDF("Normal",xplot,mu[2],si[2]) +
           pi[3]*PDF("Normal",xplot,mu[3],si[3]);

    plotdens = xplot || gmmd;
    nm1 = {"x","fx"};
    create densplot from plotdens[colname=nm1];
    append from plotdens;

    create plotd from results[colname=nm];
    append from results;

    /* ===== Create parameter dataset for plotting ===== */
    params = mu || si || pi;
    names = {"mu1" "mu2" "mu3" "si1" "si2" "si3" "pi1" "pi2" "pi3"};
    create gmm_params from params[colname=names];
    append from params;
quit;


/* ================= VISUALIZATIONS ================= */

/*proc sgplot data=plotd;
    title "Evolution of Mixing Proportions (3 Components)";
    series x=it y=pi1 / legendlabel="pi1" lineattrs=(thickness=2);
    series x=it y=pi2 / legendlabel="pi2" lineattrs=(thickness=2);
    series x=it y=pi3 / legendlabel="pi3" lineattrs=(thickness=2);
    yaxis grid;
run;

proc sgplot data=plotd;
    title "Evolution of Means (3 Components)";
    series x=it y=mu1 / legendlabel="mu1" lineattrs=(thickness=2);
    series x=it y=mu2 / legendlabel="mu2" lineattrs=(thickness=2);
    series x=it y=mu3 / legendlabel="mu3" lineattrs=(thickness=2);
    yaxis grid;
run;

proc sgplot data=plotd;
    title "Evolution of Standard Deviations (3 Components)";
    series x=it y=si1 / legendlabel="si1" lineattrs=(thickness=2);
    series x=it y=si2 / legendlabel="si2" lineattrs=(thickness=2);
    series x=it y=si3 / legendlabel="si3" lineattrs=(thickness=2);
    yaxis grid;
run; */


proc sgplot data=densplot;
    title "Estimated Density using 3-Component GMM";
    series x=x y=fx / lineattrs=(color=red thickness=2);
    xaxis grid;
    yaxis grid;
run;
title; /* clear title */



proc iml;
    use dat;

    read all var{X} into x;
    n = nrow(x);

    /* ================= CONVERGED PARAMETERS ================= */
    /* Inserted from your final EM iteration */
    pi = {0.1833133 0.2919718 0.5247149};
    mu = {110.09227 99.974039 120.0155};
    si = { 2.7433898 2.0569701 4.0361854};

    /* ================= CALCULATE RESPONSIBILITIES ================= */
    /* Calculate the unnormalized weighted density for each component */
    gamma = pi[1]*PDF("Normal", x, mu[1], si[1]) ||
            pi[2]*PDF("Normal", x, mu[2], si[2]) ||
            pi[3]*PDF("Normal", x, mu[3], si[3]);

    /* Normalize to get true probabilities (0 to 1) */
    gamma = gamma / gamma[,+];

    /* ================= HARD CLUSTERING ================= */
    cluster_id = J(n, 1, .);

    do i = 1 to n;
        /* Extract the 3 probabilities for the current observation */
        row_probs = gamma[i,];
        max_prob = max(row_probs);

        /* Assign observation to the cluster with the highest probability */
        if row_probs[1] = max_prob then cluster_id[i] = 1;
        else if row_probs[2] = max_prob then cluster_id[i] = 2;
        else cluster_id[i] = 3;
    end;

    /* ================= EXPORT DATA ================= */
    /* Combine original data, the assigned cluster, and the probabilities */
    out_data = x || cluster_id || gamma;
    col_names = {"X_Value" "Cluster" "Prob_C1" "Prob_C2" "Prob_C3"};

    create clustered_dataX from out_data[colname=col_names];
    append from out_data;
    close clustered_dataX;

    print "Hard Clustering Complete. Data saved to 'clustered_dataX'.";
quit;


/* ================= VISUALIZE CLUSTERS ================= */

/* 1. Grouped Histogram to see how the distributions overlap */
proc sgplot data=clustered_dataX;
    title "Hard Clustering Results: 3-Component GMM";
    histogram X_Value / group=Cluster transparency=0.4 binwidth=2;
    xaxis label="Observed Value";
    yaxis label="Density" grid;
    keylegend / location=inside position=topright;
run;

/* 2. Scatter plot to see exact assignments */
proc sgplot data=clustered_dataX;
    title "Observation Assignments by Cluster";
    scatter x=X_Value y=Cluster / group=Cluster markerattrs=(symbol=CircleFilled size=10);
    xaxis label="Observed Value";
    yaxis label="Assigned Cluster (1, 2, or 3)" integer grid;
run;
title; /* clear title */





proc iml;
    /* ================= 1. INPUT CONVERGED PARAMETERS ================= */
    pi_y = {0.2014942 0.3089843 0.4895215};
    mu_y = {110.19112 100.08072 120.58842};
    si_y = {2.8090352 2.1794235 7.8177596};

    pi_x = {0.1833133 0.2919718 0.5247149};
    mu_x = {110.09227 99.974039 120.01550};
    si_x = {2.7433898 2.0569701 4.0361854};

    /* ================= 2. CREATE A DENSE EVALUATION GRID ================= */
    min_val = min(mu_x - 4*si_x, mu_y - 4*si_y);
    max_val = max(mu_x + 4*si_x, mu_y + 4*si_y);

    step = (max_val - min_val) / 10000;
    v = do(min_val, max_val, step)`;
    n_points = nrow(v);

    /* ================= 3. CALCULATE THEORETICAL CDFs ================= */
    cdf_x = J(n_points, 1, 0);
    cdf_y = J(n_points, 1, 0);

    do k = 1 to ncol(pi_x);
        cdf_x = cdf_x + pi_x[k] * CDF("Normal", v, mu_x[k], si_x[k]);
    end;

    do k = 1 to ncol(pi_y);
        cdf_y = cdf_y + pi_y[k] * CDF("Normal", v, mu_y[k], si_y[k]);
    end;

    /* ================= 4. COMPUTE KOLMOGOROV-SMIRNOV STATISTIC ================= */
    diff = abs(cdf_x - cdf_y);
    KS_stat = max(diff);

    max_idx = loc(diff = KS_stat);
    v_max = v[max_idx[1]];
    cdf_x_max = cdf_x[max_idx[1]];
    cdf_y_max = cdf_y[max_idx[1]];

    print "========== GMM KOLMOGOROV-SMIRNOV TEST ==========";
    print KS_stat[format=8.5 label="KS Statistic (Max Diff)"]
          v_max[format=8.3 label="Occurs at Value"];

    /* ================= 5. FORMAT EXPORT FOR SGPLOT ================= */
    /* Create columns for the gap line filled with missing values (.) */
    gap_x_col = J(n_points, 1, .);
    gap_y_col = J(n_points, 1, .);

    /* Insert the two coordinates for the maximum gap at the very top of the columns */
    gap_x_col[1] = v_max;
    gap_y_col[1] = cdf_x_max;
    gap_x_col[2] = v_max;
    gap_y_col[2] = cdf_y_max;

    /* Combine everything into ONE dataset */
    plot_data = v || cdf_x || cdf_y || diff || gap_x_col || gap_y_col;
    cnames = {"Grid_Value" "CDF_X" "CDF_Y" "Absolute_Diff" "Gap_X" "Gap_Y"};

    create ks_plot from plot_data[colname=cnames];
    append from plot_data;
    close ks_plot;
quit;


/* ================= VISUALIZE THE KS DISTANCE ================= */
title "Theoretical CDFs for GMM X and GMM Y";
proc sgplot data=ks_plot;
    /* Plot the CDF of X */
    series x=Grid_Value y=CDF_X / lineattrs=(color=cx4F81BD thickness=2) legendlabel="GMM X CDF";

    /* Plot the CDF of Y */
    series x=Grid_Value y=CDF_Y / lineattrs=(color=cx9BBB59 thickness=2) legendlabel="GMM Y CDF";

    /* Overlay the actual maximum difference (KS Statistic) using the merged columns */
    series x=Gap_X y=Gap_Y / lineattrs=(color=red thickness=3 pattern=shortdash) legendlabel="KS Max Distance";

    xaxis label="Variable Value";
    yaxis label="Cumulative Probability";
    keylegend / location=inside position=topleft;
run;
title; /* clear title */


proc iml;
    use dat;
    /* Read raw variables */
    read all var{X} into x;
    read all var{Y} into y;

    nx = nrow(x);
    ny = nrow(y);

    /* ================= 1. FUNCTION: EMPIRICAL KS STATISTIC ================= */
    /* A highly vectorized function to calculate the empirical KS stat quickly */
    start get_ks(x_vec, y_vec);
        n1 = nrow(x_vec);
        n2 = nrow(y_vec);

        /* Combine data and create an ID tag (1 for X, 0 for Y) */
        z = x_vec // y_vec;
        id = J(n1, 1, 1) // J(n2, 1, 0);

        /* Sort the combined array while keeping track of the original IDs */
        call sortndx(idx, z, 1);
        id_sort = id[idx];

        /* Calculate empirical CDFs dynamically using cumulative sums */
        cdf_x = cusum(id_sort) / n1;
        cdf_y = cusum(1 - id_sort) / n2;

        /* Return the maximum absolute distance */
        return (max(abs(cdf_x - cdf_y)));
    finish;

    /* ================= 2. OBSERVED TEST STATISTIC ================= */
    obs_ks = get_ks(x, y);
    print "Observed Empirical KS Statistic:" obs_ks;

    /* Export the observed statistic to a macro variable for plotting later */
    call symputx("OBS_KS_MACRO", obs_ks);

        /* ================= 3. POOLED BOOTSTRAP ================= */
    z_pool = x // y;
    n_pool = nrow(z_pool);

    num_boots = 1000; /* SAFE VARIABLE NAME */
    boot_ks = J(num_boots, 1, .);

    print "Running 1000 Bootstrap Iterations... Please wait.";

    do boot_iter = 1 to num_boots; /* SAFE INDEX */
        /* Draw synthetic samples of size nx and ny from the pooled data */
        idx_x = sample(1:n_pool, nx, "replace");
        idx_y = sample(1:n_pool, ny, "replace");

        x_boot = z_pool[idx_x];
        y_boot = z_pool[idx_y];

        /* Calculate and store the KS statistic for this bootstrap iteration */
        boot_ks[boot_iter] = get_ks(x_boot, y_boot);
    end;

    /* ================= 4. P-VALUE CALCULATION ================= */
    /* The p-value is the proportion of simulated KS stats that are >= our observed stat */
    p_val = sum(boot_ks >= obs_ks) / num_boots;

    print "========== BOOTSTRAP HYPOTHESIS TEST ==========";
    print num_boots[label="Iterations"] obs_ks[label="Observed KS"] p_val[label="P-Value"];

    /* Export p-value to a macro variable for the plot title */
    call symputx("PVAL_MACRO", p_val);

    /* ================= 5. EXPORT FOR VISUALIZATION ================= */
    create boot_res from boot_ks[colname="KS_Star"];
    append from boot_ks;
    close boot_res;
quit;


/* ================= VISUALIZE THE NULL DISTRIBUTION ================= */
title "Bootstrap Null Distribution of the KS Statistic (p-value = &PVAL_MACRO)";
proc sgplot data=boot_res;
    /* Plot the distribution of simulated KS statistics */
    histogram KS_Star / fillattrs=(color=cx4F81BD) transparency=0.3;
    density KS_Star / type=kernel lineattrs=(color=cx385D8A thickness=2) legendlabel="Null Distribution";

    /* Overlay the actual observed statistic */
    refline &OBS_KS_MACRO / axis=x lineattrs=(color=red thickness=3 pattern=shortdash)
                            label="Observed KS" labelloc=inside labelpos=max;

    xaxis label="Simulated KS Statistic (D*)";
    yaxis label="Frequency" grid;
    keylegend / location=inside position=topright;
run;
title; /* clear title */



proc iml;
    /* Set a random seed for reproducibility */
    call randseed(12345);

    /* ================= 1. INPUT CONVERGED PARAMETERS ================= */
    /* Parameters for X */
    pi_x = {0.1833133 0.2919718 0.5247149};
    mu_x = {110.09227 99.974039 120.01550};
    si_x = {2.7433898 2.0569701 4.0361854};

    /* Parameters for Y */
    pi_y = {0.2014942 0.3089843 0.4895215};
    mu_y = {110.19112 100.08072 120.58842};
    si_y = {2.8090352 2.1794235 7.8177596};

    /* ================= 2. MONTE CARLO SETUP ================= */
    N_sim = 1000000; /* 1 Million simulations for extreme precision */

    /* Generate base random numbers */
    u_x = J(N_sim, 1, .); call randgen(u_x, "Uniform");
    z_x = J(N_sim, 1, .); call randgen(z_x, "Normal");

    u_y = J(N_sim, 1, .); call randgen(u_y, "Uniform");
    z_y = J(N_sim, 1, .); call randgen(z_y, "Normal");

    /* ================= 3. SIMULATE GMM FOR X ================= */
    sim_X = J(N_sim, 1, .);

    /* Vectorized assignment based on component probabilities */
    idx1 = loc(u_x <= pi_x[1]);
    if ncol(idx1) > 0 then sim_X[idx1] = mu_x[1] + si_x[1]*z_x[idx1];

    idx2 = loc(u_x > pi_x[1] & u_x <= (pi_x[1]+pi_x[2]));
    if ncol(idx2) > 0 then sim_X[idx2] = mu_x[2] + si_x[2]*z_x[idx2];

    idx3 = loc(u_x > (pi_x[1]+pi_x[2]));
    if ncol(idx3) > 0 then sim_X[idx3] = mu_x[3] + si_x[3]*z_x[idx3];

    /* Calculate Monte Carlo P(X > 100) */
    mc_prob_X = sum(sim_X > 100) / N_sim;

    /* ================= 4. SIMULATE GMM FOR Y ================= */
    sim_Y = J(N_sim, 1, .);

    idx1 = loc(u_y <= pi_y[1]);
    if ncol(idx1) > 0 then sim_Y[idx1] = mu_y[1] + si_y[1]*z_y[idx1];

    idx2 = loc(u_y > pi_y[1] & u_y <= (pi_y[1]+pi_y[2]));
    if ncol(idx2) > 0 then sim_Y[idx2] = mu_y[2] + si_y[2]*z_y[idx2];

    idx3 = loc(u_y > (pi_y[1]+pi_y[2]));
    if ncol(idx3) > 0 then sim_Y[idx3] = mu_y[3] + si_y[3]*z_y[idx3];

    /* Calculate Monte Carlo P(Y < 120) */
    mc_prob_Y = sum(sim_Y < 120) / N_sim;

    /* ================= 5. EXACT THEORETICAL PROBABILITIES ================= */
    /* P(X > 100) = 1 - CDF(100) */
    exact_X = 1 - (pi_x[1]*CDF("Normal", 100, mu_x[1], si_x[1]) +
                   pi_x[2]*CDF("Normal", 100, mu_x[2], si_x[2]) +
                   pi_x[3]*CDF("Normal", 100, mu_x[3], si_x[3]));

    /* P(Y < 120) = CDF(120) */
    exact_Y = (pi_y[1]*CDF("Normal", 120, mu_y[1], si_y[1]) +
               pi_y[2]*CDF("Normal", 120, mu_y[2], si_y[2]) +
               pi_y[3]*CDF("Normal", 120, mu_y[3], si_y[3]));

    /* ================= 6. OUTPUT RESULTS ================= */
    print "========== MONTE CARLO INTEGRATION RESULTS (N = 1,000,000) ==========";

    Result_X = mc_prob_X || exact_X || abs(mc_prob_X - exact_X);
    Result_Y = mc_prob_Y || exact_Y || abs(mc_prob_Y - exact_Y);

    cnames = {"Monte_Carlo_Est" "Exact_Theoretical" "Absolute_Error"};
    rnames_X = {"P(X > 100)"};
    rnames_Y = {"P(Y < 120)"};

    print Result_X[colname=cnames rowname=rnames_X format=8.5];
    print Result_Y[colname=cnames rowname=rnames_Y format=8.5];
quit;







































proc iml;
    use dat;
    read all var{X} into x;
    n = nrow(x);
    print "Sample Size:" n;

    /* ================= INITIAL PARAMETERS ================= */
    /* Use spaces (not commas) to create 1x2 row vectors */
    mu = {100 104};
    si = {8 4.5};
    pi = {0.5 0.5};

    stop = 0;

    /* Safely initialize an empty matrix */
    free results;

    do i = 1 to 200 while (stop=0);

        /* ===== E-step: Calculate responsibilities (gamma) ===== */
        gamma = pi[1]*PDF("Normal",x,mu[1],si[1]) ||
                pi[2]*PDF("Normal",x,mu[2],si[2]);

        gamma = gamma / gamma[,+];

        /* ===== M-step: Update parameters ===== */
        nk = gamma[+,];
        npi = nk / nk[+];
        nmu = (x#gamma)[+,] / nk;

        /* Split variance calculation for exact matrix scaling */
        nsi2_1 = (gamma[,1]#(x-nmu[1]))`*(x-nmu[1])/nk[1];
        nsi2_2 = (gamma[,2]#(x-nmu[2]))`*(x-nmu[2])/nk[2];
        nsi = sqrt(nsi2_1 || nsi2_2);

        /* ===== Check convergence ===== */
        diff = max(abs(mu - nmu));
        if diff < 0.005 then stop=1;

        /* ===== Update and store ===== */
        mu = nmu;
        si = nsi;
        pi = npi;

        /* Dynamically build the results matrix */
        if i = 1 then results = (i || npi || nmu || nsi);
        else results = results // (i || npi || nmu || nsi);
    end;

    /* ===== Print Iteration Results ===== */
    nm = {"it" "pi1" "pi2" "mu1" "mu2" "si1" "si2"};
    print "EM Convergence Results";
    print results[colname=nm];

    /* ===== Estimated log-likelihood for 2-component GMM ===== */
    ll_gmm = sum(log(
        pi[1]*PDF("Normal",x,mu[1],si[1]) +
        pi[2]*PDF("Normal",x,mu[2],si[2])
    ));

    /* ===== Null model: fit single normal ===== */
    mean0 = mean(x);
    std0 = std(x);
    ll_null = sum(log(PDF("Normal",x,mean0,std0)));

    /* ===== Model Metrics (AIC, BIC, LRT) ===== */
    k = 5; /* 2 means + 2 std devs + 1 mixing proportion */
    aic = 2*k - 2*ll_gmm;
    bic = log(n)*k - 2*ll_gmm;
    lrt = -2*(ll_null - ll_gmm);

    metrics = ll_gmm || ll_null || aic || bic || lrt;
    names = {"LL_GMM" "LL_Null" "AIC" "BIC" "LRT_Stat"};
    print "Model Fit Statistics";
    print metrics[colname=names];

    /* ===== Create estimated density plot data ===== */
    xplot = T(do(min(x), max(x), (max(x)-min(x))/200));
    gmmd = pi[1]*PDF("Normal",xplot,mu[1],si[1]) +
           pi[2]*PDF("Normal",xplot,mu[2],si[2]);

    plotdens = xplot || gmmd;
    nm1 = {"x","fx"};
    create densplot from plotdens[colname=nm1];
    append from plotdens;

    create plotd from results[colname=nm];
    append from results;
quit;


/* ================= VISUALIZATIONS ================= */

proc sgplot data=plotd;
    title "Evolution of Mixing Proportions (2 Components)";
    series x=it y=pi1 / legendlabel="pi1" lineattrs=(thickness=2);
    series x=it y=pi2 / legendlabel="pi2" lineattrs=(thickness=2);
    yaxis grid;
run;

proc sgplot data=plotd;
    title "Evolution of Means (2 Components)";
    series x=it y=mu1 / legendlabel="mu1" lineattrs=(thickness=2);
    series x=it y=mu2 / legendlabel="mu2" lineattrs=(thickness=2);
    yaxis grid;
run;

proc sgplot data=plotd;
    title "Evolution of Standard Deviations (2 Components)";
    series x=it y=si1 / legendlabel="si1" lineattrs=(thickness=2);
    series x=it y=si2 / legendlabel="si2" lineattrs=(thickness=2);
    yaxis grid;
run;

proc sgplot data=densplot;
    title "Estimated Density using 2-Component GMM";
    series x=x y=fx / lineattrs=(color=red thickness=2);
    xaxis grid;
    yaxis grid;
run;
title;
