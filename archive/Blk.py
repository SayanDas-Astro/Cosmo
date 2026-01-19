import matplotlib.pyplot as plt
import numpy as np

# DATA: [Name, Type, M_BH (Solar), M_Stellar (Solar)]
# Note: M_Stellar is the mass of the STARS only (Bulge/Galaxy), not the Halo.
# Sources: Kormendy & Ho (2013), McConnell & Ma (2013), and recent BCG papers.

data_tiebreaker = [
    # --- BCGs (Cluster Centers) ---
    ["Phoenix A*", "BCG", 1.0e11, 2.5e12],  # Massive stellar core
    ["NGC 4889",   "BCG", 2.1e10, 1.0e12],
    ["NGC 3842",   "BCG", 9.7e9,  3.5e11],
    ["M87",        "BCG", 6.5e9,  6.0e11],
    ["NGC 1399",   "BCG", 8.8e8,  3.0e11],
    ["Cygnus A",   "BCG", 2.5e9,  4.0e11],
    
    # --- Field / Ellipticals ---
    ["Sombrero",   "E",   1.0e9,  1.4e11],
    ["M60",        "E",   4.5e9,  5.5e11],
    ["M49",        "E",   2.4e9,  6.0e11],
    ["Andromeda",  "S",   1.4e8,  1.0e11], # Bulge mass approx
    ["Milky Way",  "S",   4.1e6,  5.0e10], # Bulge mass approx
    ["NGC 3377",   "E",   1.8e8,  3.0e10],
    ["NGC 3115",   "E",   2.0e9,  2.0e11],
    ["Centaurus A","E",   5.5e7,  1.0e11]
]

# Extract
names = [row[0] for row in data_tiebreaker]
types = [row[1] for row in data_tiebreaker]
m_bh = np.array([row[2] for row in data_tiebreaker])
m_star = np.array([row[3] for row in data_tiebreaker])

# Plot
plt.figure(figsize=(10, 7))

# Plot Field/Ellipticals (Blue)
mask_field = [t != "BCG" for t in types]
plt.scatter(m_star[mask_field], m_bh[mask_field], c='blue', s=80, label='Normal Galaxies (Field/E)', zorder=5)

# Plot BCGs (Red)
mask_bcg = [t == "BCG" for t in types]
plt.scatter(m_star[mask_bcg], m_bh[mask_bcg], c='red', s=150, edgecolors='k', label='Cluster Centers (BCGs)', zorder=10)

# Add "Standard" Kormendy & Ho Relation Line (approx slope 1.0 for view)
# Standard view: BH is ~0.2% to 0.5% of Bulge Mass
x_line = np.linspace(min(m_star), max(m_star), 100)
plt.plot(x_line, x_line * 0.002, 'k--', alpha=0.5, label='Standard Ratio (0.2%)')
plt.plot(x_line, x_line * 0.05, 'r--', alpha=0.5, label='Extreme Ratio (5.0%)')

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Stellar Mass $M_{*}$ ($M_{\odot}$)', fontsize=12)
plt.ylabel(r'Black Hole Mass $M_{BH}$ ($M_{\odot}$)', fontsize=12)
plt.title('The Tie-Breaker: Are BCG Black Holes "Normal" for their Stars?', fontsize=14)
plt.legend()
plt.grid(True, which="both", alpha=0.2)

for i, txt in enumerate(names):
    plt.annotate(txt, (m_star[i], m_bh[i]), xytext=(0, 10), textcoords='offset points', ha='center', fontsize=9)

plt.show()