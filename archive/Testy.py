import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# Same dataset as before
data = [
    ["Phoenix A*", "BCG", 1.0e11, 2.5e15, 'C'],
    ["Holmberg 15A","BCG", 4.0e10, 1.0e15, 'B'],
    ["IC 1101",    "BCG", 4.0e10, 6.0e14, 'C'],
    ["NGC 4889",   "BCG", 2.1e10, 9.0e14, 'A'],
    ["NGC 1600",   "BCG", 1.7e10, 3.0e14, 'B'],
    ["NGC 3842",   "BCG", 9.7e9,  4.0e14, 'A'],
    ["M87",        "BCG", 6.5e9,  1.0e14, 'A'],
    ["Cygnus A",   "BCG", 2.5e9,  6.0e14, 'B'],
    ["NGC 1399",   "BCG", 8.8e8,  8.0e13, 'B'],
    ["NGC 7768",   "BCG", 1.3e9,  2.0e14, 'B'],
    ["M60",        "E",   4.5e9,  5.0e13, 'B'],
    ["M49",        "E",   2.4e9,  4.5e13, 'B'],
    ["M84",        "E",   1.5e9,  4.0e13, 'B'],
    ["NGC 4649",   "E",   4.7e9,  6.0e13, 'B'],
    ["NGC 5846",   "E",   1.1e9,  3.0e13, 'B'],
    ["NGC 3115",   "E",   2.0e9,  2.5e12, 'B'],
    ["Centaurus A","E",   5.5e7,  2.0e13, 'B'],
    ["Sombrero",   "E",   1.0e9,  1.0e13, 'B'],
    ["Andromeda",  "S",   1.4e8,  2.0e12, 'A'],
    ["Milky Way",  "S",   4.1e6,  1.5e12, 'A'],
    ["NGC 4594",   "S",   1.0e9,  2.0e12, 'B'],
    ["M81",        "S",   7.0e7,  1.0e12, 'B'],
    ["NGC 4258",   "S",   4.0e7,  8.0e11, 'B'],
    ["NGC 1023",   "S",   4.4e7,  6.0e11, 'B'],
    ["Circinus",   "S",   1.7e6,  5.0e11, 'B'],
    ["NGC 7457",   "S",   9.0e6,  4.0e11, 'B'],
    ["M32",        "E",   2.5e6,  4.0e11, 'B'],
    ["NGC 3377",   "E",   1.8e8,  6.0e11, 'B'],
    ["NGC 3379",   "E",   4.0e8,  9.0e11, 'B'],
    ["NGC 821",    "E",   4.0e7,  7.0e11, 'B'],
    ["NGC 2778",   "E",   1.5e7,  3.0e11, 'B'],
]

# Split by type
bcg_data = [d for d in data if d[1] == 'BCG']
field_data = [d for d in data if d[1] in ['E', 'S']]

# Extract BCG data
bcg_mbh = np.array([d[2] for d in bcg_data])
bcg_mhalo = np.array([d[3] for d in bcg_data])
bcg_log_x = np.log10(bcg_mhalo)
bcg_log_y = np.log10(bcg_mbh)

# Extract field data
field_mbh = np.array([d[2] for d in field_data])
field_mhalo = np.array([d[3] for d in field_data])
field_log_x = np.log10(field_mhalo)
field_log_y = np.log10(field_mbh)

# Fit each separately
bcg_slope, bcg_int, bcg_r, bcg_p, bcg_err = stats.linregress(bcg_log_x, bcg_log_y)
field_slope, field_int, field_r, field_p, field_err = stats.linregress(field_log_x, field_log_y)

# Combined fit (for comparison)
all_log_x = np.log10(np.array([d[3] for d in data]))
all_log_y = np.log10(np.array([d[2] for d in data]))
all_slope, all_int, all_r, all_p, all_err = stats.linregress(all_log_x, all_log_y)

print("=" * 70)
print("CRITICAL TEST: DO BCGs AND FIELD GALAXIES HAVE DIFFERENT RELATIONS?")
print("=" * 70)
print("\n1. BCGs ONLY (using CLUSTER masses):")
print(f"   Slope:     {bcg_slope:.3f} ± {bcg_err:.3f}")
print(f"   R²:        {bcg_r**2:.3f}")
print(f"   P-value:   {bcg_p:.2e}")
print(f"   N:         {len(bcg_data)} galaxies")

print("\n2. FIELD + ELLIPTICALS ONLY (using GALAXY halos):")
print(f"   Slope:     {field_slope:.3f} ± {field_err:.3f}")
print(f"   R²:        {field_r**2:.3f}")
print(f"   P-value:   {field_p:.2e}")
print(f"   N:         {len(field_data)} galaxies")

print("\n3. ALL GALAXIES COMBINED:")
print(f"   Slope:     {all_slope:.3f} ± {all_err:.3f}")
print(f"   R²:        {all_r**2:.3f}")
print(f"   P-value:   {all_p:.2e}")

print("\n" + "=" * 70)
print("INTERPRETATION:")
print("=" * 70)
slope_diff = abs(bcg_slope - field_slope)
combined_err = np.sqrt(bcg_err**2 + field_err**2)
sigma_diff = slope_diff / combined_err

if sigma_diff > 2:
    print(f"⚠️  DIFFERENT RELATIONS! ({sigma_diff:.1f}σ difference)")
    print("   → BCGs follow M_BH ∝ M_cluster (cluster mass)")
    print("   → Field galaxies follow M_BH ∝ M_galaxy (galaxy halo)")
    print("   → These are DIFFERENT physical systems!")
else:
    print(f"✓  SAME RELATION ({sigma_diff:.1f}σ difference)")
    print("   → Universal M_BH-M_halo relation across all environments")
    print("   → But be careful: BCG 'halos' are actually clusters!")

print("\n" + "=" * 70)
print("WHAT DOES THIS MEAN FOR YOUR HYPOTHESIS?")
print("=" * 70)
print("If slopes are SIMILAR:")
print("  → M_BH scales with total dark matter (galaxy OR cluster)")
print("  → Suggests fundamental M_BH-halo connection")
print("  → But physical mechanism unclear (why same for different systems?)")
print("\nIf slopes are DIFFERENT:")
print("  → BCGs grow through different physics (mergers + cooling flows)")
print("  → Field galaxies grow through gas accretion")
print("  → Environment matters!")

# --- PLOT ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# Panel 1: Split view
ax1.scatter(bcg_mhalo, bcg_mbh, c='red', marker='o', s=200, 
            edgecolors='k', alpha=0.8, label='BCGs (Cluster Mass)', linewidths=2)
ax1.scatter(field_mhalo, field_mbh, c='blue', marker='s', s=100, 
            edgecolors='k', alpha=0.7, label='Field+Ellipticals (Galaxy Halo)', linewidths=1.5)

# Fit lines
bcg_x = np.linspace(min(bcg_mhalo), max(bcg_mhalo), 100)
bcg_y = 10**(bcg_slope * np.log10(bcg_x) + bcg_int)
ax1.plot(bcg_x, bcg_y, 'r--', linewidth=2.5, 
         label=f'BCG fit: slope={bcg_slope:.2f}±{bcg_err:.2f}')

field_x = np.linspace(min(field_mhalo), max(field_mhalo), 100)
field_y = 10**(field_slope * np.log10(field_x) + field_int)
ax1.plot(field_x, field_y, 'b--', linewidth=2.5, 
         label=f'Field fit: slope={field_slope:.2f}±{field_err:.2f}')

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel(r'Halo/Cluster Mass ($M_{\odot}$)', fontsize=13, fontweight='bold')
ax1.set_ylabel(r'Black Hole Mass ($M_{\odot}$)', fontsize=13, fontweight='bold')
ax1.set_title('SPLIT ANALYSIS: Are These Different Relations?', fontsize=14, fontweight='bold')
ax1.legend(loc='lower right', fontsize=10)
ax1.grid(True, which="both", ls="-", alpha=0.2)

# Panel 2: Residuals
# Calculate residuals from combined fit
bcg_predicted = 10**(all_slope * bcg_log_x + all_int)
bcg_residuals = bcg_mbh / bcg_predicted

field_predicted = 10**(all_slope * field_log_x + all_int)
field_residuals = field_mbh / field_predicted

ax2.scatter(bcg_mhalo, bcg_residuals, c='red', marker='o', s=200, 
            edgecolors='k', alpha=0.8, label='BCGs', linewidths=2)
ax2.scatter(field_mhalo, field_residuals, c='blue', marker='s', s=100, 
            edgecolors='k', alpha=0.7, label='Field+Ellipticals', linewidths=1.5)
ax2.axhline(1.0, color='k', linestyle='--', linewidth=2, alpha=0.5)
ax2.axhline(2.0, color='gray', linestyle=':', linewidth=1, alpha=0.5)
ax2.axhline(0.5, color='gray', linestyle=':', linewidth=1, alpha=0.5)

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel(r'Halo/Cluster Mass ($M_{\odot}$)', fontsize=13, fontweight='bold')
ax2.set_ylabel(r'Residual: $M_{BH,obs} / M_{BH,predicted}$', fontsize=13, fontweight='bold')
ax2.set_title('RESIDUAL TEST: Do BCGs Deviate Systematically?', fontsize=14, fontweight='bold')
ax2.legend(loc='best', fontsize=11)
ax2.grid(True, which="both", ls="-", alpha=0.2)
ax2.set_ylim(0.1, 10)

plt.tight_layout()
plt.show()

print("\n📊 Look at the residual plot (right panel):")
print("   - If BCGs cluster ABOVE 1.0 → They're overluminous (more mass than expected)")
print("   - If BCGs cluster BELOW 1.0 → They're underluminous")
print("   - If scattered around 1.0 → They follow the same relation")