import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

# For your BCGs, check if overmassiveness increases with cluster mass
bcg_data = [
    # [Name, M_BH/M_*, M_cluster (M_sun)]
    ["Phoenix A*", 0.0400, 2.4e15],  # From your earlier data
    ["NGC 4889",   0.0210, 1.2e15],
    ["NGC 3842",   0.0277, 1.8e15],
    ["M87",        0.0108, 6.4e14],
    ["NGC 1399",   0.0029, 3.0e14],
    ["Cygnus A",   0.0062, 5.0e14],
]

names_bcg = [row[0] for row in bcg_data]
overmass_ratio = np.array([row[1] for row in bcg_data])
m_cluster = np.array([row[2] for row in bcg_data])

# Correlation test
corr, p_val = spearmanr(m_cluster, overmass_ratio)

print("=" * 70)
print("CORRELATION TEST: Does Overmassiveness Scale with Cluster Mass?")
print("=" * 70)
print(f"\nSpearman Rank Correlation:")
print(f"  ρ (rho)  = {corr:.3f}")
print(f"  p-value  = {p_val:.3f}")
print(f"\nInterpretation:")
if p_val < 0.05:
    print(f"  ✓ SIGNIFICANT correlation (p < 0.05)!")
    print(f"  → More massive clusters have more overmassive black holes")
elif p_val < 0.10:
    print(f"  ~ MARGINAL correlation (p < 0.10)")
    print(f"  → Suggestive but not conclusive")
else:
    print(f"  ✗ No significant correlation (p > 0.10)")
    print(f"  → Overmassiveness doesn't strongly correlate with cluster mass")

# Show individual data points
print("\n" + "=" * 70)
print("INDIVIDUAL BCG DATA:")
print("=" * 70)
print(f"{'Galaxy':<15} {'M_BH/M_*':<12} {'M_cluster':<15} {'Overmass Factor'}")
print("-" * 70)
field_mean = 0.00467  # From earlier analysis
for i, name in enumerate(names_bcg):
    overmass_factor = overmass_ratio[i] / field_mean
    print(f"{name:<15} {overmass_ratio[i]*100:>6.2f}%      {m_cluster[i]:.2e}     {overmass_factor:.1f}×")

print("=" * 70)

# Visual plot
plt.figure(figsize=(10, 6))
plt.scatter(m_cluster, overmass_ratio * 100, s=200, c='red', edgecolors='black', linewidths=2, zorder=10)

# Add galaxy labels
for i, name in enumerate(names_bcg):
    plt.annotate(name, (m_cluster[i], overmass_ratio[i] * 100), 
                xytext=(10, 5), textcoords='offset points', fontsize=10)

# Add reference line for field galaxies
plt.axhline(y=0.467, color='blue', linestyle='--', linewidth=2, label='Field Average (0.47%)', alpha=0.7)

# Formatting
plt.xscale('log')
plt.xlabel(r'Cluster Mass $M_{\rm cluster}$ ($M_{\odot}$)', fontsize=13)
plt.ylabel(r'Black Hole Mass Fraction $M_{\rm BH}/M_*$ (%)', fontsize=13)
plt.title(f'Does Overmassiveness Correlate with Cluster Mass?\n(Spearman ρ = {corr:.3f}, p = {p_val:.3f})', fontsize=14)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print(f"\n📊 See the plot above to visualize the trend!")