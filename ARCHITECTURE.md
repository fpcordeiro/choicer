# choicer — System Architecture

**Run:** O3 MXL Hessian BLAS-3 Restructure **Date:** 2026-06-20
**Commit:** 79a4310 (worktree: agent-a19d36074ea8df700)

------------------------------------------------------------------------

## System Architecture

### Module Structure

``` mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    subgraph API["API Layer — R"]
        MNL[mnlogit_utils.R<br/>run_mnlogit / prepare_mnl_data]
        MXL[mxlogit_utils.R<br/>run_mxlogit / prepare_mxl_data]
        NL[nestlogit_utils.R<br/>run_nestlogit / prepare_nl_data]
        MNP[mnprobit_utils.R<br/>run_mnprobit / prepare_mnp_data]
    end

    subgraph Core["Core Layer — R"]
        CLS[classes.R<br/>constructors / run_optimizer<br/>compute_hessian / ensure_vcov]
        MTH[methods.R<br/>S3 generics + methods<br/>coef / vcov / predict<br/>elasticities / blp / diversion]
        WTP[wtp.R<br/>wtp — delta-method SEs]
        SRP[surplus.R<br/>logsum / consumer_surplus]
        GOF[gof.R<br/>gof — McFadden R2 / hit rate]
        PND[predict_newdata.R<br/>resolve_predict_newdata]
        SMP[sampling.R<br/>wesml_weights<br/>sample_choice_based]
        SIM[simulation.R<br/>simulate_*_data — DGPs]
        REC[recovery.R<br/>recovery_table / monte_carlo<br/>mc_asymptotics]
        UTL[utils.R<br/>shared helpers]
    end

    subgraph CPP["C++ Kernel Layer — src/"]
        MNLC[mnlogit.cpp<br/>loglik / gradient / Hessian<br/>BLP / elasticities]
        MXLC[mxlogit.cpp<br/>loglik / gradient<br/>Hessian BLAS-3<br/>BHHH / Halton]
        NLC[nestlogit.cpp<br/>loglik / gradient<br/>analytical Hessian]
        MNPC[mnprobit.cpp<br/>Gibbs sampler]
        UTLC[utils.cpp<br/>softmax / shared kernels]
    end

    subgraph HDR["C++ Headers — src/"]
        CH[choicer.h<br/>public header]
        CIH[choicer_internal.h<br/>internal structs]
        BSH[bayes_samplers.h<br/>MNP samplers]
        RNG[rng.h<br/>thread-safe RNG]
        HAL[halton.h<br/>HaltonGen — on-the-fly draws]
    end

    MNL --> CLS
    MXL --> CLS
    NL  --> CLS
    MNP --> CLS

    CLS --> MTH
    CLS --> MNLC
    CLS --> MXLC
    CLS --> NLC
    CLS --> MNPC

    MTH --> WTP
    MTH --> SRP
    MTH --> GOF
    MTH --> PND

    MNLC --> UTLC
    MXLC --> UTLC
    NLC  --> UTLC
    MNPC --> BSH
    MNPC --> RNG
    MXLC --> HAL

    MNLC --> CIH
    MXLC --> CIH
    NLC  --> CIH
    MNPC --> CIH
    CIH  --> CH

    SIM --> REC
    SMP --> MNL

    style MXLC fill:#1e90ff,stroke:#1565c0,color:#fff
    style MXL fill:#1e90ff,stroke:#1565c0,color:#fff
```

| Module | Purpose | Key Dependencies | Changed in This Run |
|----|----|----|----|
| `R/mnlogit_utils.R` | Data prep and estimation wrapper for MNL | classes.R, mnlogit.cpp | No |
| `R/mxlogit_utils.R` | Data prep and estimation wrapper for MXL; BHHH SE | classes.R, mxlogit.cpp | No |
| `R/nestlogit_utils.R` | Data prep and estimation wrapper for NL | classes.R, nestlogit.cpp | No |
| `R/mnprobit_utils.R` | Data prep and estimation wrapper for Bayesian MNP | classes.R, mnprobit.cpp | No |
| `R/classes.R` | S3 constructors, optimizer dispatch, Hessian inversion | nloptr, all C++ kernels | No |
| `R/methods.R` | S3 generics: coef, vcov, predict, elasticities, blp, diversion | classes.R | No |
| `R/wtp.R` | WTP with delta-method SEs | methods.R | No |
| `R/surplus.R` | Logsum and consumer surplus | methods.R | No |
| `R/gof.R` | Goodness-of-fit statistics | methods.R | No |
| `R/predict_newdata.R` | Counterfactual prediction helper | methods.R | No |
| `R/sampling.R` | WESML weights and choice-based sampling | utils.R | No |
| `R/simulation.R` | DGPs for MNL/MXL/NL/MNP simulation | recovery.R | No |
| `R/recovery.R` | Parameter-recovery diagnostics, Monte Carlo | simulation.R | No |
| `R/utils.R` | Shared R helpers | — | No |
| `src/mnlogit.cpp` | MNL likelihood, gradient, Hessian (block-decomposition), post-estimation | choicer_internal.h, utils.cpp | No |
| `src/mxlogit.cpp` | MXL likelihood, gradient, **Hessian (BLAS-3 restructured)**, BHHH, Halton | choicer_internal.h, halton.h, utils.cpp | **YES — `mxl_hessian_parallel` O3** |
| `src/nestlogit.cpp` | NL likelihood, gradient, analytical Hessian | choicer_internal.h, utils.cpp | No |
| `src/mnprobit.cpp` | Bayesian MNP Gibbs sampler | bayes_samplers.h, rng.h | No |
| `src/utils.cpp` | Shared C++ kernels: softmax, log-sum-exp | choicer_internal.h | No |
| `src/choicer.h` | Public C++ API header | — | No |
| `src/choicer_internal.h` | Internal structs and validation helpers | choicer.h | No |
| `src/bayes_samplers.h` | MNP-specific sampling routines | rng.h | No |
| `src/rng.h` | Thread-safe per-draw RNG | — | No |
| `src/halton.h` | On-the-fly Halton draw generator (HaltonGen) | — | No |
| `docs/mixed_logit_math.md` | MXL mathematical reference | — | **YES — §4.4 added** |

------------------------------------------------------------------------

### Function Call Graph

#### Main Pipeline: MXL Hessian (this run’s change)

``` mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    A[run_mxlogit] --> B[prepare_mxl_data]
    A --> C[run_optimizer]
    C --> D[mxl_loglik_parallel<br/>returns -ell]
    D --> C
    C --> CONV{Converged?}
    CONV -- No --> C
    CONV -- Yes --> E[compute_hessian<br/>se_method dispatch]
    E --> F[mxl_hessian_parallel]
    style F fill:#1e90ff,stroke:#1565c0,color:#fff
    F --> F0[Parse theta<br/>validate inputs<br/>WGamma batch Cholesky]
    F0 --> F1[N-loop per individual]
    F1 --> F2[S-loop: stable_softmax<br/>P_choice_s per draw]
    F2 --> F3[Alt-loop: build zc_a<br/>sum_Pzz_cc / sum_Pz_c / g_c]
    F3 --> F4[End-of-draw Identity C<br/>buf_Pzz_cc += P_s x sum_Pzz_cc<br/>buf_diff_HV_cc += P_s x sum_diff_HV]
    F3 --> F5[End-of-draw Identities A+B<br/>G_stash col s = sqrt_Ps x g_is<br/>F_stash col s = sqrt_Ps x sum_Pz_s]
    F4 --> F6[Post-S BLAS-3<br/>GF = join_rows G_stash F_stash<br/>opg_pz = GF x GF_t]
    F5 --> F6
    F6 --> F7[Block assembly hess_t1<br/>cc / cd / dc / dd blocks]
    F7 --> F8[Finalize per individual<br/>H_i = hess_t1/P_hat - g_i g_i_t<br/>local_hess += w_i H_i]
    F8 --> F9[OpenMP critical merge<br/>global_hess += local_hess]
    F9 --> G[invert_hessian -> vcov / se]
    G --> H[choicer_mxl S3 object]
    H --> POST[Post-estimation<br/>elasticities / WTP / BLP / CS / GoF]
```

| Function | Purpose | Key Dependencies | Changed in This Run |
|----|----|----|----|
| `run_mxlogit` | Public MXL estimation entry point | prepare_mxl_data, run_optimizer | No |
| `prepare_mxl_data` | Validates inputs, builds X/W matrices, Halton draws | data.table, randtoolbox | No |
| `run_optimizer` | Unified optimizer dispatch (nloptr / optim / custom) | nloptr | No |
| `mxl_loglik_parallel` | MXL log-likelihood (OpenMP over individuals) | stable_softmax | No |
| `mxl_gradient_parallel` | MXL gradient (OpenMP; per-draw score accumulation) | stable_softmax | No |
| `mxl_hessian_parallel` | **O3: BLAS-3 restructured Hessian** — block buffers + G/F stash + fused syrk | stable_softmax, arma BLAS | **YES** |
| `mxl_bhhh_parallel` | OPG/BHHH variance estimator (outer product of scores) | stable_softmax | No |
| `invert_hessian` | Invert Hessian to get vcov; fallback to pseudo-inverse | arma::inv_sympd | No |
| `new_choicer_mxl` | S3 constructor for choicer_mxl objects | — | No |

------------------------------------------------------------------------

### Data Flow

``` mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    IN[Input: data.frame<br/>id/alt/choice/covariates] --> PREP[prepare_mxl_data<br/>validate + build X/W matrices<br/>Halton draws or HaltonGen seed]
    PREP --> XMX[X / W stacked matrices<br/>eta_draws cube or gen_seed<br/>alt_idx / choice_idx / M]
    XMX --> OPT{Optimizer loop<br/>nloptr / optim}
    OPT --> GRAD[mxl_loglik_parallel<br/>returns -ell]
    GRAD --> OPT
    OPT --> CONV{Converged?}
    CONV -- No --> OPT
    CONV -- Yes --> HESS[mxl_hessian_parallel<br/>O3 BLAS-3 restructure]
    style HESS fill:#1e90ff,stroke:#1565c0,color:#fff
    HESS --> SLOOP[S-loop per individual<br/>stable_softmax / P_choice_s]
    SLOOP --> BUFACC[Identity C: buf_Pzz_cc/cd/dd<br/>buf_diff_HV_cc<br/>scalar-times-matrix accumulation]
    SLOOP --> STASH[Identities A+B: G_stash F_stash<br/>col s = sqrt_Ps x g_is or sum_Pz_s]
    BUFACC --> BLAS3[BLAS-3: GF x GF_t<br/>= G x G_t + F x F_t]
    STASH --> BLAS3
    BLAS3 --> ASSEM[Block assembly hess_t1<br/>cc/cd/dc/dd finalize H_i]
    ASSEM --> VCOV[invert_hessian -> vcov / se]
    VCOV --> OBJ[choicer_mxl S3 object]
    OBJ --> POST[Post-estimation<br/>elasticities / WTP / BLP / CS / GoF]
```

------------------------------------------------------------------------

## Key Architectural Change — O3 MXL Hessian BLAS-3 Restructure

### Before

`mxl_hessian_parallel` maintained two
$`n\_\text{params} \times n\_\text{params}`$ allocations per individual:
a per-draw dense matrix `H_is` and a running accumulator
`hess_term1_numerator`. At each draw $`s`$ it assembled the full `H_is`
from `sum_Pzz`, `sum_Pz`, and `sum_diff_H_V`, then computed
`hess_term1_numerator += P_s * (g_is * g_isᵀ + H_is)` — $`S`$ separate
rank-1 outer-product additions inside the draw loop.

### After

The restructure eliminates both per-draw allocations by splitting the
numerator into four independently handled pieces:

**Identity C (linear pieces):** `buf_Pzz_cc/cd/dd` and `buf_diff_HV_cc`
accumulate the two terms of $`H_{is}`$ that are linear in $`P_{is}`$ via
direct scalar-times-matrix additions, avoiding the dense `H_is`
assembly.

**Identities A + B (outer-product pieces):** `G_stash` and `F_stash`
(both $`n\_\text{params} \times S`$) are filled column-by-column during
the draw loop. After the loop a single BLAS level-3 multiply
`[G|F][G|F]ᵀ` delivers
$`G G^\top + F F^\top = \sum_s P_{is}(g_{is} g_{is}^\top + \mathrm{sum\_Pz}_{is}\,\mathrm{sum\_Pz}_{is}^\top)`$,
replacing $`S`$ separate rank-1 outer products with one cache-efficient
kernel.

**Identity D (optional):** Per-draw `dsyrk` over $`\tilde{Z}_c^s`$ to
batch `sum_Pzz_cc_s` — available but deferred; per-alt rank-1 loop
retained.

**No API change.** Signature, return type, sign convention, and WESML
weight handling are all unchanged. Worst-case numerical deviation
vs. oracle: $`4.588 \times 10^{-13}`$ (floating-point summation order
only), 2,000 times inside the $`10^{-6}`$ gate. All 1,371 existing tests
pass; CRAN: 0 errors / 0 warnings.

**Six buffer types replacing two allocations:**

| Buffer / Stash | Shape | Role | Identity |
|----|----|----|----|
| `buf_Pzz_cc` | $`K_c \times K_c`$ | $`\sum_s P_s\,\mathrm{sum\_Pzz}_{cc}`$ | C |
| `buf_Pzz_cd` | $`K_c \times J_d`$ | $`\sum_s P_s\,\mathrm{sum\_Pzz}_{cd}`$ | C |
| `buf_Pzz_dd` | $`J_d`$-vector | $`\sum_s P_s\,\mathrm{sum\_Pzz}_{dd}`$ diagonal | C |
| `buf_diff_HV_cc` | $`K_c \times K_c`$ | $`\sum_s P_s\,\mathrm{sum\_diff\_HV}`$ | C |
| `G_stash` | $`n\_\text{params} \times S`$ | $`G_{[:,s]} = \sqrt{P_s}\,g_{is}`$ | A |
| `F_stash` | $`n\_\text{params} \times S`$ | $`F_{[:,s]} = \sqrt{P_s}\,\mathrm{sum\_Pz}_s`$ | B |
