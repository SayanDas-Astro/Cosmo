import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# DATA: [Name, Type, M_BH (Solar), M_Stellar (Solar)]
data_tiebreaker = [
    # --- BCGs (Cluster Centers) ---
    ["Phoenix A*", "BCG", 1.0e11, 2.5e12],
    ["NGC 4889",   "BCG", 2.1e10, 1.0e12],
    ["NGC 3842",   "BCG", 9.7e9,  3.5e11],
    ["M87",        "BCG", 6.5e9,  6.0e11],
    ["NGC 1399",   "BCG", 8.8e8,  3.0e11],
    ["Cygnus A",   "BCG", 2.5e9,  4.0e11],
    
    # --- Field / Ellipticals ---
    ["Sombrero",   "E",   1.0e9,  1.4e11],
    ["M60",        "E",   4.5e9,  5.5e11],
    ["M49",        "E",   2.4e9,  6.0e11],
    ["Andromeda",  "S",   1.4e8,  1.0e11],
    ["Milky Way",  "S",   4.1e6,  5.0e10],
    ["NGC 3377",   "E",   1.8e8,  3.0e10],
    ["NGC 3115",   "E",   2.0e9,  2.0e11],
    ["Centaurus A","E",   5.5e7,  1.0e11]
]

# Extract data
names = [row[0] for row in data_tiebreaker]
types = [row[1] for row in data_tiebreaker]
m_bh = np.array([row[2] for row in data_tiebreaker])
m_star = np.array([row[3] for row in data_tiebreaker])

# Create masks
mask_field = np.array([t != "BCG" for t in types])
mask_bcg = np.array([t == "BCG" for t in types])

# Calculate M_BH/M_* ratios
ratios = m_bh / m_star

# Split by type
ratio_field = ratios[mask_field]
ratio_bcg = ratios[mask_bcg]

print("=" * 70)
print("DISCOVERY: BCG BLACK HOLES ARE OVERMASSIVE")
print("=" * 70)

print(f"\n{'Field/Elliptical Galaxies:'}")
print(f"  Sample size   = {np.sum(mask_field)}")
print(f"  Mean M_BH/M_* = {np.mean(ratio_field):.6f} ({np.mean(ratio_field)*100:.3f}%)")
print(f"  Median        = {np.median(ratio_field):.6f} ({np.median(ratio_field)*100:.3f}%)")
print(f"  Std Dev       = {np.std(ratio_field):.6f}")
print(f"  Range         = {np.min(ratio_field):.6f} to {np.max(ratio_field):.6f}")

print(f"\n{'Brightest Cluster Galaxies (BCGs):'}")
print(f"  Sample size   = {np.sum(mask_bcg)}")
print(f"  Mean M_BH/M_* = {np.mean(ratio_bcg):.6f} ({np.mean(ratio_bcg)*100:.3f}%)")
print(f"  Median        = {np.median(ratio_bcg):.6f} ({np.median(ratio_bcg)*100:.3f}%)")
print(f"  Std Dev       = {np.std(ratio_bcg):.6f}")
print(f"  Range         = {np.min(ratio_bcg):.6f} to {np.max(ratio_bcg):.6f}")

# Statistical test: Are BCGs different from field?
t_stat, p_value = stats.ttest_ind(ratio_bcg, ratio_field)

print(f"\n{'Statistical Test (Student\'s t-test):'}")
print(f"  Null hypothesis: BCGs and Field have same M_BH/M_* ratio")
print(f"  t-statistic = {t_stat:.3f}")
print(f"  p-value     = {p_value:.2e}")
print(f"  Significance: ", end="")
if p_value < 0.001:
    print(f"✓✓✓ EXTREMELY SIGNIFICANT (p < 0.001)")
elif p_value < 0.01:
    print(f"✓✓ HIGHLY SIGNIFICANT (p < 0.01)")
elif p_value < 0.05:
    print(f"✓ SIGNIFICANT (p < 0.05)")
else:
    print(f"✗ NOT SIGNIFICANT (p > 0.05)")

# Calculate overmassive factor
overmass_factor = np.mean(ratio_bcg) / np.mean(ratio_field)
median_overmass = np.median(ratio_bcg) / np.median(ratio_field)

print(f"\n{'Overmassive Factor:'}")
print(f"  Mean-based:   BCGs are {overmass_factor:.1f}× more massive than field galaxies")
print(f"  Median-based: BCGs are {median_overmass:.1f}× more massive than field galaxies")

print("\n" + "=" * 70)
print("INDIVIDUAL BCG OVERMASSIVENESS:")
print("=" * 70)
field_mean = np.mean(ratio_field)
for i, name in enumerate(names):
    if types[i] == "BCG":
        factor = ratios[i] / field_mean
        print(f"  {name:15s}: {ratios[i]*100:.2f}% (expected: {field_mean*100:.2f}%) → {factor:.1f}× overmassive")

print("=" * 70)

# Additional non-parametric test (doesn't assume normal distribution)
from scipy.stats import mannwhitneyu
u_stat, p_value_mw = mannwhitneyu(ratio_bcg, ratio_field, alternative='greater')
print(f"\nMann-Whitney U Test (non-parametric):")
print(f"  Tests if BCGs have HIGHER ratios than field")
print(f"  U-statistic = {u_stat:.1f}")
print(f"  p-value     = {p_value_mw:.2e}")
print("=" * 70)