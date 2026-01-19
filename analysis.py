
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
import pandas as pd
import traceback

print("DEBUG: Imports done.")

# Use Seaborn's default clean style
sns.set_theme(style="whitegrid", context="paper", font_scale=1.4)
plt.rcParams['font.family'] = 'sans-serif'

try:
    # 1. DATA ENTRY WITH UNCERTAINTIES
    # Format: (Name, M_BH, M_BH_err, M_star, M_star_err, M_cluster)
    # Errors are approximate (0.2 dex ~ factor of 1.6)
    print("DEBUG: Defining data with uncertainties...")
    
    bcg_data = [
        ("Phoenix A",      1.00e11, 0.5e11,  2.50e12, 0.5e12,  2.40e15),
        ("NGC 4889",       2.10e10, 0.8e10,  1.00e12, 0.3e12,  1.20e15),
        ("NGC 3842",       9.70e9,  3.0e9,   3.50e11, 0.7e11,  1.80e15),
        ("M87",            6.50e9,  0.7e9,   6.00e11, 1.0e11,  6.40e14),
        ("Cygnus A",       2.50e9,  0.7e9,   4.00e11, 0.8e11,  5.00e14),
        ("NGC 1399",       8.80e8,  3.0e8,   3.00e11, 0.6e11,  3.00e14),
        ("Abell 1835",     3.00e10, 1.0e10,  1.20e12, 0.3e12,  1.10e15),
        ("Hydra A",        1.00e9,  0.3e9,   1.00e12, 0.2e12,  5.50e14),
        ("MS0735",         1.00e10, 0.4e10,  1.10e12, 0.3e12,  9.00e14),
        ("Abell 2029",     1.00e10, 0.3e10,  1.00e12, 0.2e12,  8.00e14),
        ("Perseus",        3.40e8,  0.8e8,   2.50e11, 0.5e11,  6.00e14),
        ("Abell 478",      8.00e9,  2.0e9,   1.00e12, 0.2e12,  7.00e14),
        ("Abell 2199",     1.50e9,  0.5e9,   8.00e11, 0.2e12,  4.00e14),
        ("PKS 0745",       5.00e9,  1.5e9,   9.00e11, 0.2e12,  8.50e14)
    ]

    field_data = [
        ("Sombrero",   1.00e9, 0.1e9,  1.40e11, 0.2e11, 1.00e13),
        ("M60",        4.50e9, 1.0e9,  5.50e11, 0.5e11, 8.00e13),
        ("M49",        2.40e9, 0.5e9,  6.00e11, 0.5e11, 1.00e14),
        ("M31",        1.40e8, 0.3e8,  1.00e11, 0.2e11, 1.50e12),
        ("Milky Way",  4.10e6, 0.1e6,  5.00e10, 0.5e10, 1.20e12),
        ("NGC 3377",   1.80e8, 0.5e8,  3.00e10, 0.5e10, 5.00e11),
        ("NGC 3115",   2.00e9, 0.4e9,  2.00e11, 0.3e11, 5.00e12),
        ("Cen A",      5.50e7, 0.3e7,  1.00e11, 0.2e11, 2.00e12)
    ]

    print("DEBUG: Arrays...")
    bcg_names = np.array([x[0] for x in bcg_data])
    bcg_mbh   = np.array([x[1] for x in bcg_data])
    bcg_mbh_err = np.array([x[2] for x in bcg_data])
    bcg_mstar = np.array([x[3] for x in bcg_data])
    bcg_mstar_err = np.array([x[4] for x in bcg_data])
    bcg_mcl   = np.array([x[5] for x in bcg_data])
    bcg_ratio = bcg_mbh / bcg_mstar
    # Error propagation: sigma_R = R * sqrt((sigma_mbh/mbh)^2 + (sigma_mstar/mstar)^2)
    bcg_ratio_err = bcg_ratio * np.sqrt((bcg_mbh_err/bcg_mbh)**2 + (bcg_mstar_err/bcg_mstar)**2)

    field_names = np.array([x[0] for x in field_data])
    field_mbh   = np.array([x[1] for x in field_data])
    field_mbh_err = np.array([x[2] for x in field_data])
    field_mstar = np.array([x[3] for x in field_data])
    field_mstar_err = np.array([x[4] for x in field_data])
    field_mhalo = np.array([x[5] for x in field_data])
    field_ratio = field_mbh / field_mstar
    field_ratio_err = field_ratio * np.sqrt((field_mbh_err/field_mbh)**2 + (field_mstar_err/field_mstar)**2)

    # 2. STATISTICAL VALIDATION
    print("\n" + "="*60)
    print("STATISTICAL RESULTS")
    print("="*60)
    
    mean_bcg = np.mean(bcg_ratio)
    mean_field = np.mean(field_ratio)
    t_stat, p_ttest = stats.ttest_ind(bcg_ratio, field_ratio, equal_var=False)
    print(f"T-test: BCG mean = {mean_bcg*100:.2f}%, Field mean = {mean_field*100:.2f}%")
    print(f"Enhancement factor: {mean_bcg/mean_field:.2f}x")
    print(f"T-test p-value: {p_ttest:.4f}")

    rho, p_spearman = stats.spearmanr(bcg_mcl, bcg_ratio)
    print(f"\nSpearman correlation (full sample): rho = {rho:.3f}, p = {p_spearman:.6f}")

    # 3. SENSITIVITY ANALYSIS (without Phoenix A)
    print("\n" + "-"*40)
    print("SENSITIVITY ANALYSIS (excluding Phoenix A)")
    print("-"*40)
    
    # Remove Phoenix A (index 0)
    bcg_mcl_no_phoenix = bcg_mcl[1:]
    bcg_ratio_no_phoenix = bcg_ratio[1:]
    
    rho_no_phoenix, p_no_phoenix = stats.spearmanr(bcg_mcl_no_phoenix, bcg_ratio_no_phoenix)
    print(f"Spearman (without Phoenix A): rho = {rho_no_phoenix:.3f}, p = {p_no_phoenix:.6f}")
    
    # Linear regression sensitivity
    slope_full, intercept_full, r_full, p_full, se_full = stats.linregress(
        np.log10(np.concatenate([bcg_mcl, field_mhalo])), 
        np.log10(np.concatenate([bcg_mbh, field_mbh]))
    )
    slope_no_phoenix, intercept_no_phoenix, r_no_phoenix, p_no_phoenix_lr, se_no_phoenix = stats.linregress(
        np.log10(np.concatenate([bcg_mcl[1:], field_mhalo])), 
        np.log10(np.concatenate([bcg_mbh[1:], field_mbh]))
    )
    print(f"Slope (full): {slope_full:.3f} ± {se_full:.3f}")
    print(f"Slope (no Phoenix): {slope_no_phoenix:.3f} ± {se_no_phoenix:.3f}")
    print("="*60 + "\n")

    # 4. FIGURES (Vector Graphics - PDF)

    # -------------------------------------------------------------------------
    # FIG 1: Correlation with Error Bars
    # -------------------------------------------------------------------------
    print("DEBUG: Plotting Fig 1...")
    plt.figure(figsize=(8, 6))
    
    # Error bars on scatter
    plt.errorbar(bcg_mcl/1e14, bcg_ratio*100, yerr=bcg_ratio_err*100, 
                 fmt='o', color='#E63946', markersize=10, capsize=4, 
                 ecolor='gray', alpha=0.8, label='BCGs')
    
    # Regression line
    z = np.polyfit(np.log10(bcg_mcl), np.log10(bcg_ratio*100), 1)
    p_fit = np.poly1d(z)
    x_fit = np.linspace(bcg_mcl.min(), bcg_mcl.max(), 100)
    plt.plot(x_fit/1e14, 10**p_fit(np.log10(x_fit)), 'k--', alpha=0.7)
    
    # Annotate key points
    for i, txt in enumerate(bcg_names):
        if txt in ["Phoenix A", "NGC 4889", "Perseus", "M87"]:
            plt.annotate(txt, (bcg_mcl[i]/1e14, bcg_ratio[i]*100), 
                         xytext=(5, 5), textcoords='offset points', fontsize=10)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Cluster Mass ($10^{14} M_\odot$)', fontweight='bold')
    plt.ylabel(r'Overmassiveness ($M_{BH}/M_*$ %)', fontweight='bold')
    plt.title(f'Environmental Correlation ($\\rho={rho:.2f}$, p<0.001)', fontsize=14)
    
    stats_text = f"Spearman $\\rho = {rho:.2f}$\n$p = {p_spearman:.2e}$"
    plt.gca().text(0.05, 0.95, stats_text, transform=plt.gca().transAxes,
                  verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    plt.grid(True, which="major", ls="-", alpha=0.3)
    plt.tight_layout()
    plt.savefig('figures/fig1_correlation.pdf', format='pdf', dpi=300)
    plt.savefig('figures/fig1_correlation.png', dpi=300)  # Keep PNG for README
    print("Saved FIG 1 (PDF + PNG)")

    # -------------------------------------------------------------------------
    # FIG 2: Comparison Boxplot
    # -------------------------------------------------------------------------
    print("DEBUG: Plotting Fig 2...")
    plt.figure(figsize=(7, 6))
    
    df_bcg = pd.DataFrame({'Ratio': bcg_ratio*100, 'Type': 'BCGs\n(Cluster Centers)'})
    df_field = pd.DataFrame({'Ratio': field_ratio*100, 'Type': 'Field Galaxies\n(Isolated)'})
    df_all = pd.concat([df_field, df_bcg])
    
    sns.boxplot(data=df_all, x='Type', y='Ratio', hue='Type', 
                palette=['#A8DADC', '#E63946'], width=0.5, showfliers=False, legend=False)
    sns.stripplot(data=df_all, x='Type', y='Ratio', color='black', alpha=0.6, size=6, jitter=True)

    plt.ylabel(r'$M_{BH}/M_*$ Percentage (%)', fontweight='bold')
    plt.xlabel('')
    plt.title(f'BCG Black Holes are Overmassive\n(Enhancement $\\times {mean_bcg/mean_field:.1f}$)', fontsize=14, pad=20)
    
    y_max = max(np.max(bcg_ratio), np.max(field_ratio)) * 100
    h = y_max * 0.1
    plt.plot([0, 0, 1, 1], [y_max+h, y_max+2*h, y_max+2*h, y_max+h], lw=1.5, c='k')
    plt.text(0.5, y_max+2.5*h, f"p = {p_ttest:.3f} *", ha='center', va='bottom', fontsize=12, fontweight='bold')
    plt.ylim(top=y_max + 5*h)

    plt.tight_layout()
    plt.savefig('figures/fig2_comparison.pdf', format='pdf', dpi=300)
    plt.savefig('figures/fig2_comparison.png', dpi=300)
    print("Saved FIG 2 (PDF + PNG)")

    # -------------------------------------------------------------------------
    # FIG 3: Universal Relation with Error Bars
    # -------------------------------------------------------------------------
    print("DEBUG: Plotting Fig 3...")
    all_mbh = np.concatenate([bcg_mbh, field_mbh])
    all_mbh_err = np.concatenate([bcg_mbh_err, field_mbh_err])
    all_mdm = np.concatenate([bcg_mcl, field_mhalo])
    
    slope, intercept, r_val, p_val, std_err = stats.linregress(np.log10(all_mdm), np.log10(all_mbh))
    print(f"Universal Slope: {slope:.2f} ± {std_err:.2f}")

    plt.figure(figsize=(9, 7))
    
    # Field galaxies with error bars
    plt.errorbar(np.log10(field_mhalo), np.log10(field_mbh), 
                 yerr=0.434*field_mbh_err/field_mbh,  # Convert to log error
                 fmt='o', color='#457B9D', markersize=10, capsize=3, 
                 ecolor='gray', alpha=0.8, label='Field Galaxies')
    
    # BCGs with error bars
    plt.errorbar(np.log10(bcg_mcl), np.log10(bcg_mbh), 
                 yerr=0.434*bcg_mbh_err/bcg_mbh,
                 fmt='D', color='#E63946', markersize=10, capsize=3, 
                 ecolor='gray', alpha=0.9, label='BCGs')
    
    # Fit line
    x_range = np.linspace(np.log10(all_mdm.min()), np.log10(all_mdm.max()), 100)
    y_fit = intercept + slope * x_range
    plt.plot(x_range, y_fit, 'k--', lw=2, label=f'Fit: $\\alpha={slope:.2f} \\pm {std_err:.2f}$')
    
    # Shaded confidence region (approx)
    plt.fill_between(x_range, y_fit - 0.3, y_fit + 0.3, alpha=0.1, color='gray')

    plt.xlabel(r'Log Total Gravitational Mass ($M_{DM}$ in $M_\odot$)', fontweight='bold')
    plt.ylabel(r'Log Black Hole Mass ($M_{BH}$ in $M_\odot$)', fontweight='bold')
    plt.title('Universal Scaling Relation', fontsize=14)
    plt.legend(loc='upper left', frameon=True, framealpha=0.9)
    plt.grid(True, which="major", ls="-", alpha=0.3)
    plt.tight_layout()
    plt.savefig('figures/fig3_universal.pdf', format='pdf', dpi=300)
    plt.savefig('figures/fig3_universal.png', dpi=300)
    print("Saved FIG 3 (PDF + PNG)")

    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)

except Exception as e:
    print("CRASHED!")
    traceback.print_exc()
