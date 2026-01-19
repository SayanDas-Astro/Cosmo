
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import traceback

print("DEBUG: Imports done.")

try:
    # 1. DATA ENTRY
    print("DEBUG: Defining data...")
    bcg_data = [
        ("Phoenix A",      1.00e11, 2.50e12, 2.40e15),
        ("NGC 4889",       2.10e10, 1.00e12, 1.20e15),
        ("NGC 3842",       9.70e9,  3.50e11, 1.80e15),
        ("M87",            6.50e9,  6.00e11, 6.40e14),
        ("Cygnus A",       2.50e9,  4.00e11, 5.00e14),
        ("NGC 1399",       8.80e8,  3.00e11, 3.00e14),
        ("Abell 1835",     3.00e10, 1.20e12, 1.10e15),
        ("Hydra A",        1.00e9,  1.00e12, 5.50e14),
        ("MS0735",         1.00e10, 1.10e12, 9.00e14),
        ("Abell 2029",     1.00e10, 1.00e12, 8.00e14),
        ("Perseus",        3.40e8,  2.50e11, 6.00e14),
        ("Abell 478",      8.00e9,  1.00e12, 7.00e14),
        ("Abell 2199",     1.50e9,  8.00e11, 4.00e14),
        ("PKS 0745",       5.00e9,  9.00e11, 8.50e14)
    ]

    field_data = [
        ("Sombrero",   1.00e9, 1.40e11, 1.00e13),
        ("M60",        4.50e9, 5.50e11, 8.00e13),
        ("M49",        2.40e9, 6.00e11, 1.00e14),
        ("M31",        1.40e8, 1.00e11, 1.50e12),
        ("Milky Way",  4.10e6, 5.00e10, 1.20e12),
        ("NGC 3377",   1.80e8, 3.00e10, 5.00e11),
        ("NGC 3115",   2.00e9, 2.00e11, 5.00e12),
        ("Cen A",      5.50e7, 1.00e11, 2.00e12)
    ]

    print("DEBUG: Arrays...")
    bcg_names = np.array([x[0] for x in bcg_data])
    bcg_mbh   = np.array([x[1] for x in bcg_data])
    bcg_mstar = np.array([x[2] for x in bcg_data])
    bcg_mcl   = np.array([x[3] for x in bcg_data])
    bcg_ratio = bcg_mbh / bcg_mstar

    field_names = np.array([x[0] for x in field_data])
    field_mbh   = np.array([x[1] for x in field_data])
    field_mstar = np.array([x[2] for x in field_data])
    field_mhalo = np.array([x[3] for x in field_data])
    field_ratio = field_mbh / field_mstar

    # 2. VALIDATION
    print("DEBUG: Running t-test...")
    mean_bcg = np.mean(bcg_ratio)
    mean_field = np.mean(field_ratio)
    t_stat, p_ttest = stats.ttest_ind(bcg_ratio, field_ratio, equal_var=False)
    print(f"t-test p: {p_ttest}")

    print("DEBUG: Running Spearman...")
    rho, p_spearman = stats.spearmanr(bcg_mcl, bcg_ratio)
    print(f"Spearman p: {p_spearman}")

    # 3. FIGURES
    print("DEBUG: Plotting Fig 1...")
    plt.figure(figsize=(10, 7))
    sc = plt.scatter(bcg_mcl, bcg_ratio, c=np.log10(bcg_mbh), cmap='plasma', s=300, edgecolor='k', zorder=5)
    plt.xscale('log')
    plt.yscale('log')
    plt.title(f'Correlation (rho={rho:.2f})')
    plt.savefig('fig1_correlation.png')
    print("Saved FIG 1")

    print("DEBUG: Plotting Fig 2...")
    plt.figure(figsize=(8, 6))
    plt.boxplot([field_ratio, bcg_ratio], labels=['Field', 'BCG'])
    plt.title(f'T-test p={p_ttest:.2f}')
    plt.savefig('fig2_comparison.png')
    print("Saved FIG 2")

    print("DEBUG: Plotting Fig 3...")
    all_mbh = np.concatenate([bcg_mbh, field_mbh])
    all_mdm = np.concatenate([bcg_mcl, field_mhalo])
    slope, intercept, r_val, p_val, std_err = stats.linregress(np.log10(all_mdm), np.log10(all_mbh))
    
    print(f"Slope: {slope:.2f} +/- {std_err:.2f}")

    plt.figure(figsize=(10, 7))
    plt.scatter(all_mdm, all_mbh, label='Data Points')
    
    # Plot Fit Line
    x_range = np.linspace(min(all_mdm), max(all_mdm), 100)
    y_fit = (10**intercept) * (x_range**slope)
    
    plt.plot(x_range, y_fit, 'k--', label=f'Fit: $\\alpha={slope:.2f} \pm {std_err:.2f}$ ($R^2={r_val**2:.2f}$)')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Total Gravitational Mass $M_{\rm DM}$ [$M_\odot$]')
    plt.ylabel(r'Black Hole Mass $M_{\rm BH}$ [$M_\odot$]')
    plt.title(f'Universal Relation: $M_{{BH}} \propto M_{{DM}}^{{{slope:.2f} \pm {std_err:.2f}}}$')
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.savefig('fig3_universal.png')
    print("Saved FIG 3")

    print("DONE ALL")

except Exception as e:
    print("CRASHED!")
    traceback.print_exc()

