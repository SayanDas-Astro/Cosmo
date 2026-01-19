import numpy as np
from scipy.stats import linregress, spearmanr
import matplotlib.pyplot as plt

# ============================================================================
# DATA
# ============================================================================

bcg_data = [
    # [Name, M_BH/M_*, M_cluster (M_sun), M_stellar (M_sun)]
    ["Phoenix A*", 0.0400, 2.4e15, 2.5e12],
    ["NGC 4889",   0.0210, 1.2e15, 1.0e12],
    ["NGC 3842",   0.0277, 1.8e15, 3.5e11],
    ["M87",        0.0108, 6.4e14, 6.0e11],
    ["NGC 1399",   0.0029, 3.0e14, 3.0e11],
    ["Cygnus A",   0.0062, 5.0e14, 4.0e11],
]

names_bcg = np.array([row[0] for row in bcg_data])
overmass_ratio = np.array([row[1] for row in bcg_data])
m_cluster = np.array([row[2] for row in bcg_data])
m_stellar = np.array([row[3] for row in bcg_data])

log_m_cluster = np.log10(m_cluster)
log_ratio = np.log10(overmass_ratio)
log_m_stellar = np.log10(m_stellar)

field_mean = 0.00467

# ============================================================================
# TEST 1: LEVERAGE TEST - IS PHOENIX A* DRIVING THE STEEP SLOPE?
# ============================================================================

print("=" * 80)
print("TEST 1: LEVERAGE TEST - Is Phoenix A* Driving the Steep Slope?")
print("=" * 80)

# Fit WITH Phoenix A*
slope_full, intercept_full, r_full, p_full, stderr_full = linregress(log_m_cluster, log_ratio)

# Fit WITHOUT Phoenix A*
mask_no_phoenix = names_bcg != "Phoenix A*"
slope_reduced, intercept_reduced, r_reduced, p_reduced, stderr_reduced = linregress(
    log_m_cluster[mask_no_phoenix], log_ratio[mask_no_phoenix]
)

print(f"\n{'WITH Phoenix A* (N=6):'}")
print(f"  α = {slope_full:.3f} ± {stderr_full:.3f}")
print(f"  R² = {r_full**2:.3f}")
print(f"  p = {p_full:.2e}")

print(f"\n{'WITHOUT Phoenix A* (N=5):'}")
print(f"  α = {slope_reduced:.3f} ± {stderr_reduced:.3f}")
print(f"  R² = {r_reduced**2:.3f}")
print(f"  p = {p_reduced:.2e}")

delta_alpha = slope_full - slope_reduced
print(f"\n{'Change in slope:'} Δα = {delta_alpha:.3f}")

print("\n" + "-" * 80)
print("INTERPRETATION:")
print("-" * 80)
if abs(delta_alpha) > 0.3:
    print("  ⚠️  WARNING: Phoenix A* has STRONG LEVERAGE!")
    print("  → Slope changes significantly when removed")
    print("  → Need more massive BCGs to confirm steep slope")
    print("  → Be cautious about extrapolating to larger clusters")
elif 0.5 < slope_reduced < 0.9:
    print("  ✓ WITHOUT Phoenix A*, slope drops to theory range!")
    print("  → Phoenix A* was artificially pulling slope up")
    print("  → True relation likely α ≈ 0.6-0.8 (matches cooling flows)")
else:
    print("  → Slope remains steep even without Phoenix A*")
    print("  → Suggests genuine super-linear scaling")

# ============================================================================
# TEST 2: STELLAR MASS SCALING - WHY IS α SO STEEP?
# ============================================================================

print("\n" + "=" * 80)
print("TEST 2: STELLAR MASS SCALING - Why Is α So Steep?")
print("=" * 80)

slope_ms, intercept_ms, r_ms, p_ms, stderr_ms = linregress(log_m_cluster, log_m_stellar)

print(f"\nM_* vs M_cluster scaling in YOUR data:")
print(f"  M_* ∝ M_cluster^({slope_ms:.3f} ± {stderr_ms:.3f})")
print(f"  R² = {r_ms**2:.3f}")
print(f"  p = {p_ms:.2e}")

print(f"\nTheory prediction:")
print(f"  M_* ∝ M_cluster^(0.5-0.8)  [Standard galaxy formation]")

expected_alpha_min = 1.3 - 0.8  # L_cool ∝ M^1.3, M_* ∝ M^0.8
expected_alpha_max = 1.3 - 0.5  # L_cool ∝ M^1.3, M_* ∝ M^0.5
implied_alpha = 1.3 - slope_ms

print(f"\nExpected α range: {expected_alpha_min:.2f} to {expected_alpha_max:.2f}")
print(f"Implied α from YOUR M_* scaling: {implied_alpha:.2f}")
print(f"Observed α: {slope_full:.2f}")

print("\n" + "-" * 80)
print("INTERPRETATION:")
print("-" * 80)
if slope_ms < 0.5:
    print(f"  ✓ Your M_* scaling is SHALLOWER than expected!")
    print(f"  → This EXPLAINS why M_BH/M_* has steep slope")
    print(f"  → Your BCGs have less stellar mass than typical for their cluster")
    print(f"  → Black holes grew normally, but stars didn't form efficiently")
elif 0.5 < slope_ms < 0.8:
    print(f"  ✓ Your M_* scaling matches theory")
    print(f"  → But α = {slope_full:.2f} is STILL too steep!")
    print(f"  → Suggests: (1) New physics, OR (2) Sample bias/errors")
else:
    print(f"  → Your M_* scaling is STEEPER than expected")
    print(f"  → This makes α = {slope_full:.2f} even MORE surprising")

# ============================================================================
# TEST 3: WHAT'S THE M_BH vs M_CLUSTER SLOPE? (NOT M_BH/M_*)
# ============================================================================

print("\n" + "=" * 80)
print("TEST 3: Direct M_BH vs M_cluster Scaling")
print("=" * 80)

m_bh = overmass_ratio * m_stellar
log_m_bh = np.log10(m_bh)

slope_bh, intercept_bh, r_bh, p_bh, stderr_bh = linregress(log_m_cluster, log_m_bh)

print(f"\nDirect scaling:")
print(f"  M_BH ∝ M_cluster^({slope_bh:.3f} ± {stderr_bh:.3f})")
print(f"  R² = {r_bh**2:.3f}")
print(f"  p = {p_bh:.2e}")

print(f"\nCheck: α_BH/M* = α_BH - α_M*")
print(f"  {slope_full:.3f} = {slope_bh:.3f} - {slope_ms:.3f}")
print(f"  Difference: {abs(slope_full - (slope_bh - slope_ms)):.3f}")
if abs(slope_full - (slope_bh - slope_ms)) < 0.05:
    print("  ✓ Math checks out!")

print("\n" + "-" * 80)
print("INTERPRETATION:")
print("-" * 80)
print(f"  Theory: M_BH ∝ L_cool ∝ M_cluster^1.3")
print(f"  Your data: M_BH ∝ M_cluster^{slope_bh:.2f}")
if 1.2 < slope_bh < 1.4:
    print(f"  ✓ MATCHES cooling flow prediction!")
elif slope_bh < 1.2:
    print(f"  → Shallower than cooling flows predict")
else:
    print(f"  → Steeper than cooling flows predict")

# ============================================================================
# VISUALIZATION
# ============================================================================

fig = plt.figure(figsize=(18, 6))

# Panel 1: Leverage test
ax1 = plt.subplot(131)

# Fits
m_fit = np.logspace(14.3, 15.5, 100)
ratio_fit_full = 10**(intercept_full + slope_full * np.log10(m_fit))
ratio_fit_reduced = 10**(intercept_reduced + slope_reduced * np.log10(m_fit))

ax1.plot(m_fit, ratio_fit_full * 100, 'r--', linewidth=3, alpha=0.8, 
        label=f'With Phoenix: α={slope_full:.2f}±{stderr_full:.2f}')
ax1.plot(m_fit, ratio_fit_reduced * 100, 'b--', linewidth=3, alpha=0.8,
        label=f'Without Phoenix: α={slope_reduced:.2f}±{stderr_reduced:.2f}')

# Data
ax1.scatter(m_cluster[~mask_no_phoenix], overmass_ratio[~mask_no_phoenix] * 100, 
          s=500, c='red', edgecolors='black', linewidths=3, zorder=10, marker='*',
          label='Phoenix A*')
ax1.scatter(m_cluster[mask_no_phoenix], overmass_ratio[mask_no_phoenix] * 100,
          s=250, c='orange', edgecolors='black', linewidths=2, zorder=10,
          label='Other BCGs')

for i, name in enumerate(names_bcg):
    ax1.annotate(name, (m_cluster[i], overmass_ratio[i] * 100),
               xytext=(8, 5), textcoords='offset points', fontsize=9)

ax1.axhline(y=0.467, color='blue', linestyle=':', linewidth=2, alpha=0.5)
ax1.set_xscale('log')
ax1.set_xlabel(r'$M_{\rm cluster}$ ($M_{\odot}$)', fontsize=12)
ax1.set_ylabel(r'$M_{\rm BH}/M_*$ (%)', fontsize=12)
ax1.set_title('Leverage Test: Phoenix A* Effect', fontsize=13, fontweight='bold')
ax1.legend(fontsize=9, loc='upper left')
ax1.grid(True, alpha=0.3)

# Panel 2: M_* vs M_cluster
ax2 = plt.subplot(132)

m_stellar_fit = 10**(intercept_ms + slope_ms * np.log10(m_fit))
ax2.plot(m_fit, m_stellar_fit, 'g--', linewidth=3, alpha=0.8,
        label=f'Fit: M* ∝ M^{slope_ms:.2f}')

# Theory lines
theory_low = 10**(10) * (m_fit / 1e15)**0.5
theory_high = 10**(10) * (m_fit / 1e15)**0.8
ax2.fill_between(m_fit, theory_low, theory_high, alpha=0.2, color='gray', label='Theory: α=0.5-0.8')

ax2.scatter(m_cluster, m_stellar, s=250, c='green', edgecolors='black', linewidths=2, zorder=10)

for i, name in enumerate(names_bcg):
    ax2.annotate(name, (m_cluster[i], m_stellar[i]),
               xytext=(8, 5), textcoords='offset points', fontsize=9)

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel(r'$M_{\rm cluster}$ ($M_{\odot}$)', fontsize=12)
ax2.set_ylabel(r'$M_*$ ($M_{\odot}$)', fontsize=12)
ax2.set_title(f'Stellar Mass Scaling: α={slope_ms:.2f}', fontsize=13, fontweight='bold')
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

# Panel 3: M_BH vs M_cluster
ax3 = plt.subplot(133)

m_bh_fit = 10**(intercept_bh + slope_bh * np.log10(m_fit))
ax3.plot(m_fit, m_bh_fit, 'purple', linestyle='--', linewidth=3, alpha=0.8,
        label=f'Fit: M_BH ∝ M^{slope_bh:.2f}')

# Theory
theory_bh = 10**(7.5) * (m_fit / 1e15)**1.3
ax3.plot(m_fit, theory_bh, 'k:', linewidth=2, alpha=0.6, label='Theory: α=1.3 (cooling)')

ax3.scatter(m_cluster, m_bh, s=250, c='purple', edgecolors='black', linewidths=2, zorder=10)

for i, name in enumerate(names_bcg):
    ax3.annotate(name, (m_cluster[i], m_bh[i]),
               xytext=(8, 5), textcoords='offset points', fontsize=9)

ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlabel(r'$M_{\rm cluster}$ ($M_{\odot}$)', fontsize=12)
ax3.set_ylabel(r'$M_{\rm BH}$ ($M_{\odot}$)', fontsize=12)
ax3.set_title(f'Black Hole Scaling: α={slope_bh:.2f}', fontsize=13, fontweight='bold')
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# ============================================================================
# FINAL VERDICT
# ============================================================================

print("\n" + "=" * 80)
print("🎯 FINAL VERDICT: WHAT SHOULD YOU PUBLISH?")
print("=" * 80)

print(f"\n{'KEY RESULTS:'}")
print(f"  1. M_BH/M_* ∝ M_cluster^{slope_full:.2f} (with Phoenix)")
print(f"  2. M_BH/M_* ∝ M_cluster^{slope_reduced:.2f} (without Phoenix)")
print(f"  3. M_BH ∝ M_cluster^{slope_bh:.2f}")
print(f"  4. M_* ∝ M_cluster^{slope_ms:.2f}")

print(f"\n{'PUBLISHABILITY ASSESSMENT:'}")

if 1.2 < slope_bh < 1.4 and abs(delta_alpha) > 0.3:
    print("  ✅ PUBLISHABLE with caveat")
    print("  → M_BH scales correctly with M_cluster (matches cooling flows)")
    print("  → But M_BH/M_* slope is driven by shallow M_* scaling")
    print("  → Phoenix A* has leverage - need larger sample")
    print("\n  PAPER MESSAGE:")
    print("  'BCG black holes follow M_BH ∝ M_cluster^1.3 as predicted")
    print("   by cooling flow models. The steep M_BH/M_* ∝ M_cluster^1.2")
    print("   arises because stellar mass grows slowly with cluster mass.'")
    
elif 0.5 < slope_reduced < 0.9 and abs(delta_alpha) > 0.3:
    print("  ✅ PUBLISHABLE - Phoenix A* is special")
    print("  → Without Phoenix, α drops to theory-predicted range")
    print("  → Phoenix A* is genuinely extreme outlier")
    print("\n  PAPER MESSAGE:")
    print("  'BCG black holes scale as M_BH/M_* ∝ M_cluster^0.6-0.8,")
    print("   consistent with cooling flow regulation. Phoenix A*")
    print("   shows enhanced growth, possibly indicating a transition")
    print("   to a different accretion regime at M_cluster > 2×10^15 M_☉.'")
    
else:
    print("  ⚠️  PROCEED WITH CAUTION")
    print("  → Slope interpretation unclear")
    print("  → Need more BCGs to confirm")
    print("\n  RECOMMENDATION:")
    print("  Add 3-5 more BCGs before publishing")

print("\n" + "=" * 80)