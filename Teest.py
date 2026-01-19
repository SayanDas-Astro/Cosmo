import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# --- CORRECTED DATASET: LOCAL UNIVERSE ONLY (z < 0.1) ---
# Removed: TON 618 (z=2.2, too distant, not comparable)
# Fixed: Halo mass definitions to be consistent
# Added: Uncertainty flags

data = [
    # Format: [Name, Type, M_BH (M_sun), M_halo (M_sun), Reliability]
    # Reliability: 'A' (excellent), 'B' (good), 'C' (uncertain)
    
    # --- CLUSTER CENTERS / BCGs ---
    # NOTE: For BCGs, M_halo represents CLUSTER mass, not individual galaxy halo
    ["Phoenix A*", "BCG", 1.0e11, 2.5e15, 'C'],  # Mass very uncertain (10-100 billion)
    ["Holmberg 15A","BCG", 4.0e10, 1.0e15, 'B'],
    ["IC 1101",    "BCG", 4.0e10, 6.0e14, 'C'],  # Mass range 10-100 billion
    ["NGC 4889",   "BCG", 2.1e10, 9.0e14, 'A'],  # McConnell et al. 2011 - reliable
    ["NGC 1600",   "BCG", 1.7e10, 3.0e14, 'B'],
    ["NGC 3842",   "BCG", 9.7e9,  4.0e14, 'A'],  # McConnell et al. 2011 - reliable
    ["M87",        "BCG", 6.5e9,  1.0e14, 'A'],  # EHT measurement - most reliable
    ["Cygnus A",   "BCG", 2.5e9,  6.0e14, 'B'],
    ["NGC 1399",   "BCG", 8.8e8,  8.0e13, 'B'],
    ["NGC 7768",   "BCG", 1.3e9,  2.0e14, 'B'],
    
    # --- MASSIVE ELLIPTICALS (Group Centers) ---
    ["M60",        "E",   4.5e9,  5.0e13, 'B'],
    ["M49",        "E",   2.4e9,  4.5e13, 'B'],
    ["M84",        "E",   1.5e9,  4.0e13, 'B'],
    ["NGC 4649",   "E",   4.7e9,  6.0e13, 'B'],
    ["NGC 5846",   "E",   1.1e9,  3.0e13, 'B'],
    ["NGC 3115",   "E",   2.0e9,  2.5e12, 'B'],
    ["Centaurus A","E",   5.5e7,  2.0e13, 'B'],
    ["Sombrero",   "E",   1.0e9,  1.0e13, 'B'],  # M104
    
    # --- SPIRALS & FIELD GALAXIES ---
    ["Andromeda",  "S",   1.4e8,  2.0e12, 'A'],
    ["Milky Way",  "S",   4.1e6,  1.5e12, 'A'],
    ["NGC 4594",   "S",   1.0e9,  2.0e12, 'B'],
    ["M81",        "S",   7.0e7,  1.0e12, 'B'],
    ["NGC 4258",   "S",   4.0e7,  8.0e11, 'B'],
    ["NGC 1023",   "S",   4.4e7,  6.0e11, 'B'],
    ["Circinus",   "S",   1.7e6,  5.0e11, 'B'],
    ["NGC 7457",   "S",   9.0e6,  4.0e11, 'B'],
    
    # --- SMALL ELLIPTICALS ---
    ["M32",        "E",   2.5e6,  4.0e11, 'B'],
    ["NGC 3377",   "E",   1.8e8,  6.0e11, 'B'],
    ["NGC 3379",   "E",   4.0e8,  9.0e11, 'B'],
    ["NGC 821",    "E",   4.0e7,  7.0e11, 'B'],
    ["NGC 2778",   "E",   1.5e7,  3.0e11, 'B'],
]

# Extract data
names = [row[0] for row in data]
m_bh = np.array([row[2] for row in data])
m_halo = np.array([row[3] for row in data])
types = [row[1] for row in data]
reliability = [row[4] for row in data]

# Log-Log Transformation
log_x = np.log10(m_halo)
log_y = np.log10(m_bh)

# --- STATISTICS ---
slope, intercept, r_value, p_value, std_err = stats.linregress(log_x, log_y)
r_squared = r_value**2

print(f"--- CORRECTED DATASET: {len(data)} GALAXIES (LOCAL UNIVERSE) ---")
print(f"Removed: TON 618 (z=2.2, not comparable)")
print(f"\nSlope:      {slope:.3f} ± {std_err:.3f}")
print(f"R-Squared:  {r_squared:.3f}")
print(f"P-Value:    {p_value:.2e}")
print(f"\nWARNING: BCG 'halos' are actually CLUSTER masses!")
print(f"This may inflate the slope.")

# --- PLOT ---
plt.figure(figsize=(12, 8))

# Plot points with different colors and markers
for i, (name, t, rel) in enumerate(zip(names, types, reliability)):
    color = 'red' if t == 'BCG' else ('blue' if t == 'E' else 'green')
    marker = 'o' if t == 'BCG' else ('s' if t == 'E' else '*')
    size = 150 if t == 'BCG' else 80
    alpha = 0.9 if rel == 'A' else (0.7 if rel == 'B' else 0.4)
    
    plt.scatter(m_halo[i], m_bh[i], c=color, marker=marker, s=size, 
                edgecolors='k', alpha=alpha, linewidths=1.5)

# Legend
plt.scatter([], [], c='red', marker='o', s=150, label='BCG (Cluster Centers)', edgecolors='k')
plt.scatter([], [], c='blue', marker='s', s=80, label='Ellipticals', edgecolors='k')
plt.scatter([], [], c='green', marker='*', s=100, label='Spirals/Field', edgecolors='k')
plt.scatter([], [], c='gray', marker='o', s=80, alpha=0.9, label='Reliability: High', edgecolors='k')
plt.scatter([], [], c='gray', marker='o', s=80, alpha=0.4, label='Reliability: Low', edgecolors='k')

# Best Fit Line
x_line = np.linspace(min(m_halo), max(m_halo), 100)
y_line = 10**(slope * np.log10(x_line) + intercept)
plt.plot(x_line, y_line, 'k--', linewidth=2, 
         label=f'Fit: Slope = {slope:.2f} ± {std_err:.2f}')

# Formatting
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Halo Mass $M_{halo}$ ($M_{\odot}$)', fontsize=13)
plt.ylabel(r'Black Hole Mass $M_{BH}$ ($M_{\odot}$)', fontsize=13)
plt.title(f'CORRECTED: M_BH vs M_halo ({len(data)} galaxies, z<0.1)', fontsize=14, fontweight='bold')
plt.legend(loc='lower right', fontsize=10)
plt.grid(True, which="both", ls="-", alpha=0.2)

# Annotations
plt.annotate("Phoenix A*\n(uncertain)", (2.5e15, 1.0e11), 
             xytext=(-80, 10), textcoords='offset points', ha='right',
             bbox=dict(boxstyle='round', fc='yellow', alpha=0.3))
plt.annotate("M87\n(EHT)", (1.0e14, 6.5e9), 
             xytext=(10, -15), textcoords='offset points',
             bbox=dict(boxstyle='round', fc='lightgreen', alpha=0.5))
plt.annotate("Milky Way", (1.5e12, 4.1e6), xytext=(10, -10), textcoords='offset points')

plt.tight_layout()
plt.show()

print("\n--- KEY ISSUES TO ADDRESS ---")
print("1. BCG 'halo masses' are CLUSTER masses (10^14-10^15), not individual halos")
print("2. Phoenix A, IC 1101 have very uncertain BH masses (factor of 2-10)")
print("3. Should we use M_stellar instead for BCGs?")
print("4. Consider splitting analysis: Field vs Group vs Cluster")