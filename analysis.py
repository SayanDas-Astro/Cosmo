
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
import traceback

print("DEBUG: Imports done.")

# Use Seaborn's default clean style
sns.set_theme(style="whitegrid", context="paper", font_scale=1.4)
plt.rcParams['font.family'] = 'sans-serif' # clean sans-serif like Arial/Helvetica

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

    # -------------------------------------------------------------------------
    # FIG 1: Correlation (Simplified)
    # -------------------------------------------------------------------------
    print("DEBUG: Plotting Fig 1...")
    plt.figure(figsize=(8, 6))
    
    # Regression plot with confidence interval helps visualize the trend immediately
    # We need a dataframe for seaborn regplot usually, but arrays work too
    sns.regplot(x=bcg_mcl/1e14, y=bcg_ratio*100, color="#E63946", scatter_kws={'s': 100, 'edgecolor':'k'}, line_kws={'color': 'black', 'ls':'--'})
    
    # Annotate limited set to avoid clutter
    for i, txt in enumerate(bcg_names):
        if txt in ["Phoenix A", "NGC 4889", "Perseus", "M87"]:
            plt.annotate(txt, (bcg_mcl[i]/1e14, bcg_ratio[i]*100), 
                         xytext=(5, 5), textcoords='offset points', fontsize=11, fontweight='bold')

    plt.xscale('log')
    plt.yscale('log')
    
    # Clear Axis Labels
    plt.xlabel(r'Cluster Mass ($10^{14} M_\odot$)', fontweight='bold')
    plt.ylabel(r'Overmassiveness ($M_{BH}/M_*$ %)', fontweight='bold')
    plt.title(f'Cluster Environment drives Black Hole Growth\n(Correlation $\\rho={rho:.2f}$)', fontsize=14)
    
    # Add stats box
    stats_text = f"Spearman $\\rho = {rho:.2f}$\n$p = {p_spearman:.4f}$"
    plt.gca().text(0.05, 0.95, stats_text, transform=plt.gca().transAxes,
                  verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    plt.grid(True, which="major", ls="-", alpha=0.5)
    plt.tight_layout()
    plt.savefig('figures/fig1_correlation.png', dpi=300)
    print("Saved FIG 1")

    # -------------------------------------------------------------------------
    # FIG 2: Comparison (Boxplot + Stripplot)
    # -------------------------------------------------------------------------
    print("DEBUG: Plotting Fig 2...")
    plt.figure(figsize=(7, 6))
    
    # Create dataset for seaborn
    import pandas as pd
    df_bcg = pd.DataFrame({'Ratio': bcg_ratio*100, 'Type': 'BCGs\n(Cluster Centers)'})
    df_field = pd.DataFrame({'Ratio': field_ratio*100, 'Type': 'Field Galaxies\n(Isolated)'})
    df_all = pd.concat([df_field, df_bcg])
    
    # Violin plot nicely shows distribution shape + data points
    sns.boxplot(data=df_all, x='Type', y='Ratio', hue='Type', palette=['#A8DADC', '#E63946'], width=0.5, showfliers=False, legend=False)
    sns.stripplot(data=df_all, x='Type', y='Ratio', color='black', alpha=0.6, size=6, jitter=True)

    plt.ylabel(r'$M_{BH}/M_*$ Percentage (%)', fontweight='bold')
    plt.xlabel('')
    plt.title(f'BCG Black Holes are Overmassive\n(Enhancement $\\times {mean_bcg/mean_field:.1f}$)', fontsize=14, pad=20)
    
    # Add p-value bracket
    y_max = max(np.max(bcg_ratio), np.max(field_ratio)) * 100
    h = y_max * 0.1
    plt.plot([0, 0, 1, 1], [y_max+h, y_max+2*h, y_max+2*h, y_max+h], lw=1.5, c='k')
    plt.text(0.5, y_max+2.5*h, f"p = {p_ttest:.3f} *", ha='center', va='bottom', color='k', fontsize=12, fontweight='bold')

    # Extend y-axis to make room for the bracket
    plt.ylim(top=y_max + 5*h)

    plt.tight_layout()
    plt.savefig('figures/fig2_comparison.png', dpi=300)
    print("Saved FIG 2")

    # -------------------------------------------------------------------------
    # FIG 3: Universal Relation (Log-Log Regression)
    # -------------------------------------------------------------------------
    print("DEBUG: Plotting Fig 3...")
    all_mbh = np.concatenate([bcg_mbh, field_mbh])
    all_mdm = np.concatenate([bcg_mcl, field_mhalo])
    
    slope, intercept, r_val, p_val, std_err = stats.linregress(np.log10(all_mdm), np.log10(all_mbh))
    print(f"Slope: {slope:.2f} +/- {std_err:.2f}")

    plt.figure(figsize=(9, 7))
    
    # Use regplot for the fit + confidence interval (95%)
    # Combine data for regression
    sns.regplot(x=np.log10(all_mdm), y=np.log10(all_mbh), scatter=False, 
                color="black", line_kws={'linestyle':'--'}, label=f'Universal Fit (Slope = {slope:.2f} $\\pm$ {std_err:.2f})')
    
    # Plot points separately to distinguish types
    plt.scatter(np.log10(field_mhalo), np.log10(field_mbh), 
                c='#457B9D', s=120, alpha=0.8, edgecolors='w', label='Field Galaxies (Halo Mass)')
    plt.scatter(np.log10(bcg_mcl), np.log10(bcg_mbh), 
                c='#E63946', marker='D', s=120, alpha=0.9, edgecolors='k', label='BCGs (Cluster Mass)')

    plt.xlabel(r'Log Total Gravitational Mass ($M_{DM}/M_{500}$ in $M_\odot$)', fontweight='bold')
    plt.ylabel(r'Log Black Hole Mass ($M_{BH}$ in $M_\odot$)', fontweight='bold')
    plt.title('Universal Scaling across 4 orders of magnitude', fontsize=14)
    plt.legend(loc='upper left', frameon=True, framealpha=0.9)
    
    plt.grid(True, which="major", ls="-", alpha=0.5)
    plt.tight_layout()
    plt.savefig('figures/fig3_universal.png', dpi=300)
    print("Saved FIG 3")

    print("DONE ALL")

except Exception as e:
    print("CRASHED!")
    traceback.print_exc()
