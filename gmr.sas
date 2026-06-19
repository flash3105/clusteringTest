proc import
        out = dat
        datafile =  "C:\Users\flash\Documents\2026\Masters\MVA\2026 Exam\gmrdat2026.xlsx"
        dbms = xlsx
        replace;
        getnames = yes;

run;


proc sgplot data=dat;
    scatter x=x y=y;
    title "Scatter plot of data";
run;

ods graphics off;
proc reg data=dat;
    model y = x;
quit;


proc iml;

use dat;
read all into xy;

X =J(nrow(xy),1,1) || xy[,1];
y = xy[,2];

/*from the plot , it looks like there are two lines , so we will assume two mixtures*/

n=nrow(X);

si=std(y); /*initial std*/
pi ={0.5 , 0.5}; /*naive approach*/
beta={20000 0.5 , 35000 -1.5} ;


mixtures = nrow(pi);


int_pi = pi;
int_si = si;
int_beta = beta;

print "*****************************start em********************";

stop =0;

entropy_store = J(n, 250, 0);


do i = 1 to 250 while (stop = 0);

        gamma = J(n,mixtures,.);   /*responsibilities*/

        do r = 1 to n;
                do j = 1 to mixtures;
                        gamma[r,j] = pi[j,1] * PDF("Normal" , y[r,1],X[r,]*beta[j,]`,si);

                end;
        end;

        gamma = gamma / gamma[,+]; /*normalization*/

        if i =1 then do;
                gamma_print = gamma[1:50,];
                print gamma_print;
        end;
        /* compute entropy each iteration*/

        entropy = J(n,1,0);

        do r=1 to n;
                do c=1 to mixtures;
                        if gamma[r,c] >0 then  do;

                                entropy[r] = entropy[r] -gamma[r,c]*log(gamma[r,c]);
                        end;

                 end;
        end;

        entropy_store[,i] =entropy;

        nk = gamma[+,];
        new_pi = nk/nk[+];
        new_beta = J(1, ncol(beta), 1); /*placeholder*/


        do j = 1 to mixtures;
                W = diag(gamma[, j]);
                beta_k = inv(X`*W*X)*X`*W*y;
                new_beta = new_beta // beta_k`;
        end;

        new_beta = new_beta[2:nrow(new_beta), ];  /*trimout the placeholder*/


        nsi2 = 0;

        do j= 1 to n;
                do c=1 to mixtures;
                v = gamma[j,c]*((y[j,1]-X[j,]*new_beta[c,]`)**2); /*new variance building blocks*/
                nsi2 = nsi2 + v;

                end;
        end;

        nsi2 = nsi2/n;   /*new variance*/

        nsi = sqrt(nsi2);

        diff = max(abs(pi- new_pi`),abs(si-nsi),abs(beta-new_beta)); /*for convergence criterion*/

        if diff < 0.0005 then stop=1;
        si = nsi;
        pi = new_pi`;
        beta= new_beta;
        final_iter =i;
        final_gamma = gamma;

        results =  results // (i || nsi || new_pi || new_beta[1,] || new_beta[2,]);

end;

/* Final E-step using converged parameters */

final_gamma = J(n, mixtures, .);

do r = 1 to n;
    do j = 1 to mixtures;
        final_gamma[r,j] =
            pi[j,1] *
            pdf("Normal",
                y[r,1],
                X[r,]*beta[j,]`,
                si);
    end;
end;

final_gamma = final_gamma / final_gamma[,+];
y_fit1 = X * beta[1,]`;
y_fit2 = X * beta[2,]`;

y_pred =
      final_gamma[,1] # y_fit1
    + final_gamma[,2] # y_fit2;

SST = sum((y - mean(y))##2);
SSE = sum((y - y_pred)##2);

R2 = 1 - SSE/SST;

print  R2 ;

k = 5; /* Number of parameters: 2 intercepts, 2 slopes, 1 variance */

/* Calculate AIC and BIC (for model comparison) */
log_likelihood = sum(log(final_gamma[,+])); /* Approximate log-likelihood */
AIC = -2 * log_likelihood + 2 * k;
BIC = -2 * log_likelihood + k * log(n);

print log_likelihood AIC BIC;




nm = {"it" "si" "pi1" "pi2" "beta1_b0" "beta1_b1" "beta2_b0" "beta2_b1"};
results = (0 || int_si || int_pi` || int_beta[1,] || int_beta[2,]) // results;

final_entropy = entropy_store[, final_iter];

cutoff = mean(final_entropy);
uncertain = (final_entropy > cutoff);

plot_data = xy || final_entropy || uncertain;
nm_plot = {"x" "y" "entropy" "uncertain"};
create plot_data from plot_data[colname=nm_plot];
append from plot_data;
close plot_data;



print results[colname = nm];
create convergence_plot from results[colname = nm];
append from results;



data plot_data;
    set plot_data;
    if entropy < 0.2 then ent_grp = 1;
    else if entropy < 0.5 then ent_grp = 2;
    else if entropy < 0.8 then ent_grp = 3;
    else ent_grp = 4;
run;

proc sgplot data=plot_data;
    scatter x=x y=y / group=ent_grp markerattrs=(symbol=circlefilled);
    title "GMR Model with Shannon Entropy (Uncertainty Levels)";
run;






/* Create combined dataset with both data and fitted lines */
data combined;
    set dat;
    /* Calculate fitted values for both components */
    y_fit1 = 107033.21 - 32.47851 * x;
    y_fit2 = 1923.5811 + 9.3140032 * x;
run;

/* Plot */
proc sgplot data=combined;
    scatter x=x y=y /
        markerattrs=(symbol=circlefilled color=green)
        transparency=0.4;

    series x=x y=y_fit1 /
        lineattrs=(color=red thickness=2)
        legendlabel="Component 1";

    series x=x y=y_fit2 /
        lineattrs=(color=blue thickness=2)
        legendlabel="Component 2";

    xaxis label="X";
    yaxis label="Y";
    title "Final GMR Model: Two Regression Lines";
    keylegend / position=right;
run;


proc sort data=plot_data;
    by x;
run;

proc sort data=combined;
    by x;
run;

/* Now merge them */
data final_plot;
    merge plot_data combined (keep=x y_fit1 y_fit2);
    by x;
run;


/* Merge and create labels */
data final_plot;
    merge plot_data combined (keep=x y_fit1 y_fit2);
    by x;

    /* Create meaningful labels for entropy groups */
    if ent_grp = 1 then entropy_label = "Low (Certain)";
    else if ent_grp = 2 then entropy_label = "Medium";
    else if ent_grp = 3 then entropy_label = "High";
    else if ent_grp = 4 then entropy_label = "Very High";
run;

proc sgplot data=final_plot;
    /* Data points colored by entropy with labels */
    scatter x=x y=y /
        group=entropy_label
        markerattrs=(symbol=circlefilled size=7)
        transparency=0.3
        name="scatter"
        legendlabel="Uncertainty Level";

    /* Fitted lines */
    series x=x y=y_fit1 /
        lineattrs=(color=yellow thickness=3 pattern=solid)
        legendlabel="Component 1"
        name="line1";

    series x=x y=y_fit2 /
        lineattrs=(color=darkblue thickness=3 pattern=solid)
        legendlabel="Component 2"
        name="line2";

    xaxis label="X Variable"
          labelattrs=(weight=bold)
          grid;
    yaxis label="Y Variable"
          labelattrs=(weight=bold)
          grid;

    title "GMR Model with Shannon Entropy (Uncertainty Levels)"
          font='Arial Bold'
          size=14;

    keylegend "scatter" "line1" "line2" /
        position=right
        title="Legend"
        titleattrs=(weight=bold)
        opaque
        border;
run;


    rowSum = final_gamma[,+] + 1e-12;
    final_gamma = final_gamma / rowSum;

    y_fit1 = X_boot * beta[1,]`;
    y_fit2 = X_boot * beta[2,]`;

    y_pred = final_gamma[,1] # y_fit1 + final_gamma[,2] # y_fit2;

    /* ================= SAFE R2 ================= */
    SST = sum((y_boot - mean(y_boot))##2);
    SSE = sum((y_boot - y_pred)##2);

    R2 = 1 - SSE/SST;

    r_boot[b] = R2;

    print "B=" b "Completed, R2=" R2 "Iterations=" i diff stop;

end;  /* end bootstrap loop */

/* ================= OUTPUT ================= */
print r_boot;

quit;




proc iml;
use dat;
read all into xy;

X = J(nrow(xy),1,1) || xy[,1];
y = xy[,2];
n = nrow(X);

/* ================= INITIAL PARAMETERS ================= */
si = std(y);
mix_probs = {0.5, 0.5}; /* renamed from pi */
beta = {20000 0.5, 35000 -1.5};
mixtures = nrow(mix_probs);

/* store initial values */
int_pi = mix_probs;
int_si = si;
int_beta = beta;

/* ================= BOOTSTRAP SETTINGS ================= */
num_boots = 1000; /* safe variable name */
r_boot = J(num_boots, 1, .);

print "************ START BOOTSTRAP EM ************";

/* ================= BOOTSTRAP LOOP ================= */
do boot_iter = 1 to num_boots; /* Safe index and boundary */

   /* print "START BOOTSTRAP ITERATION:" boot_iter; */


    /* reset params */
    mix_probs = int_pi;
    si = int_si;
    beta = int_beta;

    /* bootstrap sample */
    idx = sample(1:n, n, "replace");
    X_boot = X[idx,];
    y_boot = y[idx,];
    n_boot = nrow(X_boot);

    /* ================= EM ================= */
    conv_flag = 0; /* safe flag */
    i = 0;

    do while (conv_flag = 0);
        i = i + 1;

        /* ===== Vectorized E-step ===== */
        gamma = J(n_boot, mixtures, 0);
        do j = 1 to mixtures;
            mu = X_boot * beta[j,]`;
            gamma[,j] = mix_probs[j,1] * pdf("Normal", y_boot, mu, si);
        end;

        /* normalize safely */
        rowSum = gamma[,+] + 1e-12;
        gamma = gamma / rowSum;

        /* ===== Vectorized M-step ===== */
        nk = gamma[+,];
        new_pi = nk / sum(nk);

        /* ================= FIXED BETA UPDATE ================= */
        new_beta = J(mixtures, ncol(X_boot), .);

        do j = 1 to mixtures;
            W_vec = gamma[,j];

            XtWX = X_boot` * (W_vec # X_boot);
            XtWy = X_boot` * (W_vec # y_boot);

            /* ridge regularization */
            beta_j = inv(XtWX + 1e-6*I(ncol(X_boot))) * XtWy;
            new_beta[j, ] = beta_j`;
        end;

        /* ===== Vectorized Variance update ===== */
        nsi2 = 0;
        do j = 1 to mixtures;
            mu = X_boot * new_beta[j,]`;
            nsi2 = nsi2 + sum(gamma[,j] # (y_boot - mu)##2);
        end;
        nsi = sqrt(nsi2 / n_boot);

        /* ===== convergence ===== */
        diff_pi   = max(abs(mix_probs - new_pi`));
        diff_si   = abs(si - nsi);
        diff_beta = max(abs(beta - new_beta));

        diff = max(diff_pi // diff_si // diff_beta);

        /* prevent early stopping and max iterations */
        if (i > 10 & diff < 1e-4) | i >= 250 then do;
            conv_flag = 1;
        end;

        /* update */
        mix_probs = new_pi`;
        beta = new_beta;
        si = max(nsi, 1e-6);

    end; /* end EM loop */

    /* ================= FINAL PREDICTION ================= */
    final_gamma = J(n_boot, mixtures, 0);
    do j = 1 to mixtures;
        mu = X_boot * beta[j,]`;
        final_gamma[,j] = mix_probs[j,1] * pdf("Normal", y_boot, mu, si);
    end;

    rowSum = final_gamma[,+] + 1e-12;
    final_gamma = final_gamma / rowSum;

    y_fit1 = X_boot * beta[1,]`;
    y_fit2 = X_boot * beta[2,]`;
    y_pred = final_gamma[,1] # y_fit1 + final_gamma[,2] # y_fit2;

    /* ================= SAFE R2 ================= */
    SST = sum((y_boot - mean(y_boot))##2);
    SSE = sum((y_boot - y_pred)##2);
    R2 = 1 - SSE/SST;

    r_boot[boot_iter] = R2;

  /* print "Bootstrap" boot_iter "Completed, R2=" R2 "Iterations=" i "Diff=" diff;*/

end; /* end bootstrap loop */


/* ================= OUTPUT ================= */
print "************ FINAL BOOTSTRAP RESULTS ************";

/* 1. Print up to 10 results safely */
n_print = min(nrow(r_boot), 10);
sample_r_boot = r_boot[1:n_print, 1];
print "First" n_print "Bootstrap R2 Values:" sample_r_boot;

/* 2. Calculate summary statistics */
mean_R2 = mean(r_boot);
std_R2 = std(r_boot);

/* 3. Construct 95% Confidence Intervals */
/* Method A: Percentile CI (Recommended for R-squared) */
call qntl(ci_percentile, r_boot, {0.025, 0.975});
ci_perc_lower = ci_percentile[1];
ci_perc_upper = ci_percentile[2];

/* Method B: Normal Approximation CI */
z_val = quantile("Normal", 0.975);
ci_norm_lower = mean_R2 - (z_val * std_R2);
ci_norm_upper = mean_R2 + (z_val * std_R2);

/* 4. Display Final Summary */
print "Bootstrap Mean R2:" mean_R2;
print "Bootstrap Std Error:" std_R2;

print "--- 95% Confidence Intervals ---";
print "Percentile Method:" ci_perc_lower ci_perc_upper;
print "Normal Approx:    " ci_norm_lower ci_norm_upper;

/* ================= SAVE TO DATASET ================= */
/* Create a dataset named 'boot_results' from the r_boot matrix */
create boot_results from r_boot[colname="R2"];
append from r_boot;
close boot_results;

quit;


/* ================= PLOT HISTOGRAM ================= */
title "Distribution of Bootstrap R-Squared Values";
proc sgplot data=boot_results;
    /* The main histogram */
    histogram R2 / fillattrs=(color=cx4F81BD) transparency=0.3;

    /* Overlay a normal distribution curve (Red) */
    density R2 / type=normal lineattrs=(color=cxC0504D thickness=2 pattern=dash) legendlabel="Normal Curve";

    /* Overlay a kernel density estimate curve (Green) */
    density R2 / type=kernel lineattrs=(color=cx9BBB59 thickness=2) legendlabel="Kernel Density";

    /* Axes and formatting */
    xaxis label="Bootstrap R-Squared (R²)";
    yaxis label="Percent" grid;
    keylegend / location=inside position=topright across=1;
run;
title; /* clear title */






proc iml;
use dat;
read all into xy;
close dat; /* Good practice to close the dataset */

/* ================= 1. DEFINE QUADRATIC DESIGN MATRIX ================= */
/* We add the squared term: 1, x, and x^2 */
X = J(nrow(xy),1,1) || xy[,1] || (xy[,1]##2);
y = xy[,2];

n = nrow(X);

/* ================= 2. INITIALIZATION ================= */
si = std(y); /* initial std */
pi = {0.5, 0.5}; /* naive approach */

/* IMPORTANT: We must now initialize three parameters per component: b0, b1, b2 */
/* The x^2 coefficient is initialized close to 0 to prevent immediate divergence */
beta = {20000  0.5  0.0001,
        35000 -1.5 -0.0001};

mixtures = nrow(pi);

int_pi = pi;
int_si = si;
int_beta = beta;

print "***************************** start em ********************";

stop = 0;
entropy_store = J(n, 250, 0);

/* ================= 3. EM ALGORITHM LOOP ================= */
do i = 1 to 250 while (stop = 0);

    /* --- E-Step --- */
    gamma = J(n, mixtures, .); /* responsibilities */

    do r = 1 to n;
        do j = 1 to mixtures;
            gamma[r,j] = pi[j,1] * PDF("Normal" , y[r,1], X[r,]*beta[j,]`, si);
        end;
    end;

    gamma = gamma / gamma[,+]; /* normalization */

    if i = 1 then do;
        gamma_print = gamma[1:50,];
        print gamma_print;
    end;

    /* Compute entropy each iteration */
    entropy = J(n, 1, 0);
    do r = 1 to n;
        do c = 1 to mixtures;
            if gamma[r,c] > 0 then do;
                entropy[r] = entropy[r] - gamma[r,c]*log(gamma[r,c]);
            end;
        end;
    end;
    entropy_store[,i] = entropy;

    /* --- M-Step --- */
    nk = gamma[+,];
    new_pi = nk / nk[+];

    new_beta = J(1, ncol(beta), 1); /* placeholder */

    do j = 1 to mixtures;
        W = diag(gamma[, j]);
        /* The exact same WLS formula automatically handles the 3x3 matrix inversion */
        beta_k = inv(X`*W*X)*X`*W*y;
        new_beta = new_beta // beta_k`;
    end;

    new_beta = new_beta[2:nrow(new_beta), ]; /* trim out the placeholder */

    /* Variance Update */
    nsi2 = 0;
    do j = 1 to n;
        do c = 1 to mixtures;
            v = gamma[j,c]*((y[j,1] - X[j,]*new_beta[c,]`)**2);
            nsi2 = nsi2 + v;
        end;
    end;
    nsi2 = nsi2 / n;
    nsi = sqrt(nsi2);

    /* --- Convergence Check --- */
    diff = max(abs(pi - new_pi`), abs(si - nsi), abs(beta - new_beta));

    if diff < 0.0005 then stop = 1;
    si = nsi;
    pi = new_pi`;
    beta = new_beta;
    final_iter = i;
    final_gamma = gamma;

    /* Notice: Added new_beta[1,3] and new_beta[2,3] for the output log */
    results = results // (i || nsi || new_pi || new_beta[1,] || new_beta[2,]);

end;

/* ================= 4. FINAL EVALUATION ================= */
/* Final E-step using converged parameters */
final_gamma = J(n, mixtures, .);
do r = 1 to n;
    do j = 1 to mixtures;
        final_gamma[r,j] = pi[j,1] * pdf("Normal", y[r,1], X[r,]*beta[j,]`, si);
    end;
end;
final_gamma = final_gamma / final_gamma[,+];

y_fit1 = X * beta[1,]`;
y_fit2 = X * beta[2,]`;

y_pred = final_gamma[,1] # y_fit1 + final_gamma[,2] # y_fit2;

SST = sum((y - mean(y))##2);
SSE = sum((y - y_pred)##2);
R2 = 1 - SSE/SST;

print R2;

/* Number of parameters: 2 intercepts, 2 linear slopes, 2 quadratic slopes, 1 variance */
k = 7;

/* Calculate AIC and BIC */
log_likelihood = sum(log(final_gamma[,+]));
AIC = -2 * log_likelihood + 2 * k;
BIC = -2 * log_likelihood + k * log(n);

print log_likelihood AIC BIC;

/* ================= 5. PLOTTING DATA ================= */
nm = {"it" "si" "pi1" "pi2" "beta1_b0" "beta1_b1" "beta1_b2" "beta2_b0" "beta2_b1" "beta2_b2"};
results = (0 || int_si || int_pi` || int_beta[1,] || int_beta[2,]) // results;

final_entropy = entropy_store[, final_iter];

cutoff = mean(final_entropy);
uncertain = (final_entropy > cutoff);

plot_data = xy || final_entropy || uncertain || y_fit1 || y_fit2 || y_pred;
nm_plot = {"x" "y" "entropy" "uncertain" "y_fit1" "y_fit2" "y_pred"};

create plot_data from plot_data[colname=nm_plot];
append from plot_data;
close plot_data;

print results[colname = nm];
create convergence_plot from results[colname = nm];
append from results;
close convergence_plot;
quit;

/* ================= SGPLOT FOR QUADRATIC GMR ================= */
title "Quadratic Gaussian Mixture Regression (K=2)";
proc sgplot data=plot_data;
    scatter x=x y=y / colorresponse=entropy colormodel=(red green) name="data";
    series x=x y=y_fit1 / lineattrs=(color=red thickness=2) name="comp1" legendlabel="Component 1 (Quadratic)";
    series x=x y=y_fit2 / lineattrs=(color=blue thickness=2) name="comp2" legendlabel="Component 2 (Quadratic)";
    xaxis label="X Variable";
    yaxis label="Y Variable";
run;
title;




/* 1. Sort the data by X first to ensure the series plot draws smooth curves */
proc sort data=dat out=dat_sorted;
    by x;
run;

/* 2. Create combined dataset with both data and fitted quadratic curves */
data combined;
    set dat_sorted;

    /* Calculate fitted values for both components using the quadratic parameters */
    y_fit1 = -343099.43 + 303.09743 * x - 0.062401 * (x**2);
    y_fit2 = 236115.03 - 156.4234 * x + 0.029287 * (x**2);
run;

/* 3. Plot */
proc sgplot data=combined;
    scatter x=x y=y /
        markerattrs=(symbol=circlefilled color=green)
        transparency=0.4
        name="scatter";

    series x=x y=y_fit1 /
        lineattrs=(color=red thickness=2)
        legendlabel="Component 1 (Quadratic)"
        name="comp1";

    series x=x y=y_fit2 /
        lineattrs=(color=blue thickness=2)
        legendlabel="Component 2 (Quadratic)"
        name="comp2";

    xaxis label="X";
    yaxis label="Y";
    title "Final Quadratic GMR Model: Two Parabolic Curves";
    keylegend "comp1" "comp2" / position=right title="Components";
run;
title;






proc iml;
use dat;
read all into xy;

/* 1. ADD SQUARED TERM TO DESIGN MATRIX */
X = J(nrow(xy),1,1) || xy[,1] || (xy[,1]##2);
y = xy[,2];
n = nrow(X);

/* ================= INITIAL PARAMETERS ================= */
si = std(y);
mix_probs = {0.5, 0.5};

/* 2. EXPAND INITIAL BETA MATRIX (Adding small quadratic starting weights) */
beta = {20000  0.5  0.0001,
        35000 -1.5 -0.0001};
mixtures = nrow(mix_probs);

/* store initial values */
int_pi = mix_probs;
int_si = si;
int_beta = beta;

/* ================= BOOTSTRAP SETTINGS ================= */
num_boots = 1000;
r_boot = J(num_boots, 1, .);

print "************ START QUADRATIC BOOTSTRAP EM ************";

/* ================= BOOTSTRAP LOOP ================= */
do boot_iter = 1 to num_boots;

    /* reset params */
    mix_probs = int_pi;
    si = int_si;
    beta = int_beta;

    /* bootstrap sample */
    idx = sample(1:n, n, "replace");
    X_boot = X[idx,];
    y_boot = y[idx,];
    n_boot = nrow(X_boot);

    /* ================= EM ================= */
    conv_flag = 0;
    i = 0;

    do while (conv_flag = 0);
        i = i + 1;

        /* ===== Vectorized E-step ===== */
        gamma = J(n_boot, mixtures, 0);
        do j = 1 to mixtures;
            mu = X_boot * beta[j,]`;
            gamma[,j] = mix_probs[j,1] * pdf("Normal", y_boot, mu, si);
        end;

        /* normalize safely */
        rowSum = gamma[,+] + 1e-12;
        gamma = gamma / rowSum;

        /* ===== Vectorized M-step ===== */
        nk = gamma[+,];
        new_pi = nk / sum(nk);

        /* ================= FIXED BETA UPDATE ================= */
        new_beta = J(mixtures, ncol(X_boot), .);

        do j = 1 to mixtures;
            W_vec = gamma[,j];
            XtWX = X_boot` * (W_vec # X_boot);
            XtWy = X_boot` * (W_vec # y_boot);

            /* ridge regularization to prevent singular matrices */
            beta_j = inv(XtWX + 1e-6*I(ncol(X_boot))) * XtWy;
            new_beta[j, ] = beta_j`;
        end;

        /* ===== Vectorized Variance update ===== */
        nsi2 = 0;
        do j = 1 to mixtures;
            mu = X_boot * new_beta[j,]`;
            nsi2 = nsi2 + sum(gamma[,j] # (y_boot - mu)##2);
        end;
        nsi = sqrt(nsi2 / n_boot);

        /* ===== convergence ===== */
        diff_pi   = max(abs(mix_probs - new_pi`));
        diff_si   = abs(si - nsi);
        diff_beta = max(abs(beta - new_beta));

        diff = max(diff_pi // diff_si // diff_beta);

        /* prevent early stopping and max iterations */
        if (i > 10 & diff < 1e-4) | i >= 250 then do;
            conv_flag = 1;
        end;

        /* update */
        mix_probs = new_pi`;
        beta = new_beta;
        si = max(nsi, 1e-6);

    end; /* end EM loop */

    /* ================= FINAL PREDICTION ================= */
    final_gamma = J(n_boot, mixtures, 0);
    do j = 1 to mixtures;
        mu = X_boot * beta[j,]`;
        final_gamma[,j] = mix_probs[j,1] * pdf("Normal", y_boot, mu, si);
    end;

    rowSum = final_gamma[,+] + 1e-12;
    final_gamma = final_gamma / rowSum;

    y_fit1 = X_boot * beta[1,]`;
    y_fit2 = X_boot * beta[2,]`;
    y_pred = final_gamma[,1] # y_fit1 + final_gamma[,2] # y_fit2;

    /* ================= SAFE R2 ================= */
    SST = sum((y_boot - mean(y_boot))##2);
    SSE = sum((y_boot - y_pred)##2);
    R2 = 1 - SSE/SST;

    r_boot[boot_iter] = R2;

end; /* end bootstrap loop */

/* ================= OUTPUT ================= */
print "************ FINAL BOOTSTRAP RESULTS ************";

/* 1. Print up to 10 results safely */
n_print = min(nrow(r_boot), 10);
sample_r_boot = r_boot[1:n_print, 1];
print "First" n_print "Bootstrap R2 Values:" sample_r_boot;

/* 2. Calculate summary statistics */
mean_R2 = mean(r_boot);
std_R2 = std(r_boot);

/* 3. Construct 95% Confidence Intervals */
/* Method A: Percentile CI */
call qntl(ci_percentile, r_boot, {0.025, 0.975});
ci_perc_lower = ci_percentile[1];
ci_perc_upper = ci_percentile[2];

/* Method B: Normal Approximation CI */
z_val = quantile("Normal", 0.975);
ci_norm_lower = mean_R2 - (z_val * std_R2);
ci_norm_upper = mean_R2 + (z_val * std_R2);

/* 4. Display Final Summary */
print "Bootstrap Mean R2:" mean_R2;
print "Bootstrap Std Error:" std_R2;
print "--- 95% Confidence Intervals ---";
print "Percentile Method:" ci_perc_lower ci_perc_upper;
print "Normal Approx:    " ci_norm_lower ci_norm_upper;

/* ================= SAVE TO DATASET ================= */
create boot_results from r_boot[colname="R2"];
append from r_boot;
close boot_results;

quit;

/* ================= PLOT HISTOGRAM ================= */
title "Distribution of Quadratic Bootstrap R-Squared Values";
proc sgplot data=boot_results;
    histogram R2 / fillattrs=(color=cx4F81BD) transparency=0.3;
    density R2 / type=normal lineattrs=(color=cxC0504D thickness=2 pattern=dash) legendlabel="Normal Curve";
    density R2 / type=kernel lineattrs=(color=cx9BBB59 thickness=2) legendlabel="Kernel Density";
    xaxis label="Bootstrap R-Squared (R²)";
    yaxis label="Percent" grid;
    keylegend / location=inside position=topright across=1;
run;
title;
