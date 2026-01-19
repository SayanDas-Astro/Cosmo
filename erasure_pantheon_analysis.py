import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp, quad
from scipy.optimize import minimize, differential_evolution
import matplotlib.pyplot as plt
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PART 1: COSMOLOGY MODELS
# ============================================================================

class CosmologyModel:
    """Base class for cosmology models"""
    def __init__(self, H0=70.0):
        self.H0 = H0
        self.c = 299792.458  # km/s
        
    def E(self, z):
        """Dimensionless Hubble parameter E(z) = H(z)/H0"""
        raise NotImplementedError
        
    def luminosity_distance(self, z):
        """Calculate luminosity distance in Mpc"""
        if z < 1e-5:
            return self.c / self.H0 * z  # Linear approximation for small z
        
        try:
            # Comoving distance: D_C = (c/H0) * integral(dz'/E(z'))
            integral, _ = quad(lambda zp: 1.0/self.E(zp), 0, z, 
                              limit=100, epsabs=1e-8, epsrel=1e-8)
            D_C = (self.c / self.H0) * integral
            
            # Luminosity distance
            D_L = D_C * (1 + z)
            return D_L
        except:
            return np.nan
    
    def distance_modulus(self, z):
        """Calculate distance modulus μ = 5*log10(D_L/Mpc) + 25"""
        D_L = self.luminosity_distance(z)
        if D_L <= 0 or np.isnan(D_L):
            return np.nan
        return 5 * np.log10(D_L) + 25


class LambdaCDM(CosmologyModel):
    """Standard ΛCDM cosmology"""
    def __init__(self, Om0=0.3, H0=70.0):
        super().__init__(H0)
        self.Om0 = Om0
        self.Ode0 = 1.0 - Om0
        self.name = "ΛCDM"
        
    def E(self, z):
        """E(z) = sqrt(Ω_m(1+z)^3 + Ω_Λ)"""
        return np.sqrt(self.Om0 * (1+z)**3 + self.Ode0)


class PowerLawErasure(CosmologyModel):
    """
    Power-Law Interacting Dark Energy Model
    Q = β H ρ_m (1+z)^α = β H ρ_m a^(-α)
    """
    def __init__(self, Om0=0.3, beta=0.05, alpha=-2.0, H0=70.0):
        super().__init__(H0)
        self.Om0 = Om0
        self.Ode0 = 1.0 - Om0
        self.beta = beta
        self.alpha = alpha
        self.name = f"Erasure(β={beta:.3f},α={alpha:.2f})"
        self.failed = False  # Track if model failed
        
        # Pre-compute evolution
        self._compute_evolution()
        
    def _compute_evolution(self):
        """Solve coupled ODEs to get Ω_m(z) and Ω_DE(z)"""
        def derivatives(a, y):
            Om, Ode = y
            
            # Prevent negative densities and enforce physical bounds
            Om = np.clip(Om, 1e-10, 10.0)
            Ode = np.clip(Ode, 1e-10, 10.0)
            
            # E(a) = H(a)/H0
            E_squared = Om / a**3 + Ode
            if E_squared < 1e-10:
                return [0, 0]
            E = np.sqrt(E_squared)
            
            # Power-Law Interaction: Q = β H ρ_m a^(-α)
            # The actual physical density is ρ_m = ρ_crit * Ω_m / a^3
            # So Q/ρ_crit = β H Ω_m a^(-3-α)
            
            # Matter continuity: d(ρ_m a^3)/da = -Q a^2 / H
            # This becomes: dΩ_m/da = -3Ω_m/a - (Q/ρ_crit) a^3 / (H a)
            #                        = -3Ω_m/a - β Ω_m a^(-α) a^(-1) / E
            
            # Protect against extreme values of a^(-alpha)
            if self.alpha < 0 and a < 0.01:
                # At very small a (high z), a^(-alpha) can explode
                # Cap the interaction strength
                a_alpha_term = min((a ** (-self.alpha)), 1e6)
            else:
                a_alpha_term = a ** (-self.alpha)
            
            Q_term = (self.beta * Om * a_alpha_term) / (a**4 * E)
            
            # Limit Q_term to prevent numerical instabilities
            Q_term = np.clip(Q_term, -100 * Om / a, 100 * Om / a)
            
            dOm_da = -(3 * Om / a) - Q_term
            dOde_da = Q_term  # Energy conservation
            
            return [dOm_da, dOde_da]
        
        # Solve from high redshift to today
        a_span = [0.01, 1.0]  # z=99 to z=0 (avoid extreme high-z instabilities)
        a_eval = np.logspace(-2, 0, 500)  # Less dense but more stable
        
        # Initial conditions at z=99 (not z=999 to avoid instabilities)
        z_init = 99
        a_init = 1.0 / (1 + z_init)
        
        # At high z, universe should be matter-dominated
        Om_init = 0.99
        Ode_init = 0.01
        
        y0 = [Om_init, Ode_init]
        
        try:
            sol = solve_ivp(derivatives, a_span, y0, 
                           t_eval=a_eval, method='LSODA',  # More robust for stiff equations
                           rtol=1e-6, atol=1e-8,  # Relaxed tolerances for stability
                           max_step=0.01)  # Limit step size
            
            if not sol.success:
                raise ValueError(f"ODE solver failed: {sol.message}")
            
            # Store as functions of redshift
            self.z_grid = 1.0/sol.t - 1
            self.Om_grid = sol.y[0]
            self.Ode_grid = sol.y[1]
            
            # Check for negative densities
            if np.any(self.Om_grid < 0) or np.any(self.Ode_grid < 0):
                raise ValueError("Negative densities encountered")
            
            # Check for energy conservation (in physical density)
            # Total density should remain close to ρ_crit
            total_omega = self.Om_grid / (1 + self.z_grid)**3 + self.Ode_grid
            max_deviation = np.max(np.abs(total_omega - 1.0))
            
            # Only warn if deviation is severe
            if max_deviation > 0.5:
                raise ValueError(f"Energy not conserved: max deviation = {max_deviation:.2f}")
            
        except Exception as e:
            # If model fails, it's unphysical - return huge chi-squared
            # by setting to ΛCDM but with wrong parameters
            self.z_grid = 1.0/a_eval - 1
            self.Om_grid = np.full_like(a_eval, 0.01)  # Wrong parameters
            self.Ode_grid = np.full_like(a_eval, 0.01)  # Will give bad chi2
            self.failed = True
            return
        
    def E(self, z):
        """Interpolate E(z) from pre-computed evolution"""
        try:
            # Find Ω_m(z) and Ω_DE(z) by interpolation
            Om_z = np.interp(z, self.z_grid, self.Om_grid, left=self.Om_grid[0], right=self.Om_grid[-1])
            Ode_z = np.interp(z, self.z_grid, self.Ode_grid, left=self.Ode_grid[0], right=self.Ode_grid[-1])
            
            return np.sqrt(Om_z * (1+z)**3 + Ode_z)
        except:
            # Fallback to ΛCDM if interpolation fails
            return np.sqrt(self.Om0 * (1+z)**3 + self.Ode0)


# ============================================================================
# PART 2: DATA LOADING
# ============================================================================

def load_pantheon_data(filepath):
    """Load Pantheon+ SH0ES data"""
    print("Loading Pantheon+ data...")
    
    try:
        data = pd.read_csv(filepath, sep=r'\s+', comment='#')
        
        # Apply quality cuts
        # 1. Exclude Cepheid calibrators (used to set distance scale)
        # 2. Require z > 0.01 (avoid peculiar velocity contamination)
        # 3. Require valid errors
        
        mask = (
            (data['IS_CALIBRATOR'] == 0) &
            (data['zHD'] > 0.01) &
            (data['MU_SH0ES_ERR_DIAG'] > 0) &
            (data['MU_SH0ES_ERR_DIAG'] < 1.0)  # Exclude outliers
        )
        
        df = data[mask].copy()
        
        print(f"  Total SNe in file: {len(data)}")
        print(f"  After quality cuts: {len(df)}")
        print(f"  Redshift range: {df['zHD'].min():.3f} - {df['zHD'].max():.3f}")
        
        return df['zHD'].values, df['MU_SH0ES'].values, df['MU_SH0ES_ERR_DIAG'].values
        
    except Exception as e:
        print(f"ERROR loading data: {e}")
        raise


# ============================================================================
# PART 3: CHI-SQUARED ANALYSIS
# ============================================================================

def chi_squared(model, z_data, mu_data, mu_err):
    """Calculate chi-squared for a model"""
    chi2 = 0.0
    n_valid = 0
    
    for z, mu_obs, err in zip(z_data, mu_data, mu_err):
        mu_model = model.distance_modulus(z)
        
        if np.isnan(mu_model):
            continue
            
        chi2 += ((mu_obs - mu_model) / err)**2
        n_valid += 1
    
    return chi2, n_valid


def fit_lcdm(z_data, mu_data, mu_err):
    """Fit ΛCDM model"""
    print("\nFitting ΛCDM...")
    
    def objective(params):
        Om0, H0 = params
        model = LambdaCDM(Om0=Om0, H0=H0)
        chi2, _ = chi_squared(model, z_data, mu_data, mu_err)
        return chi2
    
    # Use differential_evolution for global optimization
    result = differential_evolution(
        objective,
        bounds=[(0.1, 0.5), (60.0, 80.0)],
        seed=42,
        maxiter=100,
        workers=1,
        disp=True,
        atol=0.01,
        tol=0.01
    )
    
    best_Om0, best_H0 = result.x
    best_model = LambdaCDM(Om0=best_Om0, H0=best_H0)
    chi2, n_valid = chi_squared(best_model, z_data, mu_data, mu_err)
    
    dof = n_valid - 2  # 2 free parameters
    
    print(f"  Best fit: Ωm = {best_Om0:.4f}, H0 = {best_H0:.2f} km/s/Mpc")
    print(f"  χ² = {chi2:.2f} ({n_valid} SNe)")
    print(f"  DoF = {dof}")
    print(f"  χ²/DoF = {chi2/dof:.4f}")
    
    return best_model, chi2, dof


def fit_erasure(z_data, mu_data, mu_err):
    """Fit Power-Law Erasure model"""
    print("\nFitting Power-Law Erasure Model...")
    print("  (This may take several minutes...)")
    
    def objective(params):
        Om0, beta, alpha, H0 = params
        
        try:
            model = PowerLawErasure(Om0=Om0, beta=beta, alpha=alpha, H0=H0)
            
            # If model construction failed, return huge chi-squared
            if model.failed:
                return 1e10
            
            chi2, _ = chi_squared(model, z_data, mu_data, mu_err)
            
            if np.isnan(chi2) or np.isinf(chi2):
                return 1e10
            return chi2
        except:
            return 1e10
    
    # Use differential_evolution for global search
    result = differential_evolution(
        objective,
        bounds=[
            (0.1, 0.5),      # Om0
            (-0.5, 0.5),     # beta
            (-6.0, 0.0),     # alpha (must be negative for late-time dominance)
            (60.0, 80.0)     # H0
        ],
        seed=42,
        maxiter=200,
        workers=1,
        disp=True,
        atol=0.01,
        tol=0.01
    )
    
    best_Om0, best_beta, best_alpha, best_H0 = result.x
    best_model = PowerLawErasure(Om0=best_Om0, beta=best_beta, alpha=best_alpha, H0=best_H0)
    chi2, n_valid = chi_squared(best_model, z_data, mu_data, mu_err)
    
    dof = n_valid - 4  # 4 free parameters
    
    print(f"  Best fit: Ωm = {best_Om0:.4f}, β = {best_beta:.4f}, α = {best_alpha:.2f}, H0 = {best_H0:.2f}")
    print(f"  χ² = {chi2:.2f} ({n_valid} SNe)")
    print(f"  DoF = {dof}")
    print(f"  χ²/DoF = {chi2/dof:.4f}")
    
    return best_model, chi2, dof


# ============================================================================
# PART 4: STATISTICAL COMPARISON
# ============================================================================

def model_comparison(lcdm_model, lcdm_chi2, lcdm_dof,
                    erasure_model, erasure_chi2, erasure_dof):
    """Compare models using multiple criteria"""
    
    print("\n" + "="*70)
    print("STATISTICAL COMPARISON")
    print("="*70)
    
    # Chi-squared comparison
    delta_chi2 = lcdm_chi2 - erasure_chi2
    delta_params = 2  # Erasure has 2 extra parameters (beta, alpha)
    
    print(f"\n1. Chi-Squared Test:")
    print(f"   ΛCDM:    χ² = {lcdm_chi2:.2f}, χ²/DoF = {lcdm_chi2/lcdm_dof:.4f}")
    print(f"   Erasure: χ² = {erasure_chi2:.2f}, χ²/DoF = {erasure_chi2/erasure_dof:.4f}")
    print(f"   Δχ² = {delta_chi2:.2f} (positive favors Erasure)")
    
    # Rule of thumb: Δχ² > 2*Δparams is "significant"
    # Δχ² > 10*Δparams is "strong evidence"
    threshold_weak = 2 * delta_params
    threshold_strong = 10 * delta_params
    
    if delta_chi2 > threshold_strong:
        verdict = f"✓ STRONG evidence for Erasure (Δχ² > {threshold_strong})"
    elif delta_chi2 > threshold_weak:
        verdict = f"→ WEAK evidence for Erasure (Δχ² > {threshold_weak})"
    elif delta_chi2 > 0:
        verdict = "→ Models comparable (extra parameters not justified)"
    else:
        verdict = "✗ ΛCDM preferred (simpler is better)"
    
    print(f"   Interpretation: {verdict}")
    
    # AIC (Akaike Information Criterion): AIC = χ² + 2k
    # Lower AIC is better
    k_lcdm = 2
    k_erasure = 4
    aic_lcdm = lcdm_chi2 + 2 * k_lcdm
    aic_erasure = erasure_chi2 + 2 * k_erasure
    delta_aic = aic_lcdm - aic_erasure
    
    print(f"\n2. Akaike Information Criterion (AIC):")
    print(f"   ΛCDM:    AIC = {aic_lcdm:.2f}")
    print(f"   Erasure: AIC = {aic_erasure:.2f}")
    print(f"   ΔAIC = {delta_aic:.2f} (positive favors Erasure)")
    
    if delta_aic > 10:
        aic_verdict = "✓ Erasure strongly favored"
    elif delta_aic > 4:
        aic_verdict = "→ Erasure moderately favored"
    elif delta_aic > 0:
        aic_verdict = "→ Models comparable"
    else:
        aic_verdict = "✗ ΛCDM favored"
    
    print(f"   Interpretation: {aic_verdict}")
    
    # BIC (Bayesian Information Criterion): BIC = χ² + k*ln(n)
    # Lower BIC is better (penalizes extra parameters more than AIC)
    n = lcdm_dof + k_lcdm  # Total number of data points
    bic_lcdm = lcdm_chi2 + k_lcdm * np.log(n)
    bic_erasure = erasure_chi2 + k_erasure * np.log(n)
    delta_bic = bic_lcdm - bic_erasure
    
    print(f"\n3. Bayesian Information Criterion (BIC):")
    print(f"   ΛCDM:    BIC = {bic_lcdm:.2f}")
    print(f"   Erasure: BIC = {bic_erasure:.2f}")
    print(f"   ΔBIC = {delta_bic:.2f} (positive favors Erasure)")
    
    if delta_bic > 10:
        bic_verdict = "✓ Erasure strongly favored"
    elif delta_bic > 6:
        bic_verdict = "→ Erasure moderately favored"
    elif delta_bic > 0:
        bic_verdict = "→ Models comparable"
    else:
        bic_verdict = "✗ ΛCDM favored"
    
    print(f"   Interpretation: {bic_verdict}")
    
    return delta_chi2, delta_aic, delta_bic


# ============================================================================
# PART 5: VISUALIZATION
# ============================================================================

def plot_results(z_data, mu_data, mu_err, lcdm_model, erasure_model, 
                lcdm_chi2, lcdm_dof, erasure_chi2, erasure_dof):
    """Create comprehensive visualization"""
    
    # Compute model predictions
    z_plot = np.linspace(0.01, max(z_data), 300)
    mu_lcdm = [lcdm_model.distance_modulus(z) for z in z_plot]
    mu_erasure = [erasure_model.distance_modulus(z) for z in z_plot]
    
    # Compute residuals
    residuals_lcdm = []
    residuals_erasure = []
    z_residuals = []
    
    for z, mu_obs in zip(z_data, mu_data):
        mu_l = lcdm_model.distance_modulus(z)
        mu_e = erasure_model.distance_modulus(z)
        
        if not np.isnan(mu_l) and not np.isnan(mu_e):
            residuals_lcdm.append(mu_obs - mu_l)
            residuals_erasure.append(mu_obs - mu_e)
            z_residuals.append(z)
    
    # Create figure
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(3, 2, height_ratios=[2, 1, 1], hspace=0.3, wspace=0.3)
    
    # Main distance modulus plot
    ax1 = fig.add_subplot(gs[0, :])
    ax1.errorbar(z_data, mu_data, yerr=mu_err, fmt='o', 
                 color='gray', alpha=0.3, markersize=2, linewidth=0.5,
                 label=f'Pantheon+ (N={len(z_data)})')
    ax1.plot(z_plot, mu_lcdm, 'k--', linewidth=2.5,
             label=f'ΛCDM (χ²/DoF={lcdm_chi2/lcdm_dof:.3f})')
    ax1.plot(z_plot, mu_erasure, 'm-', linewidth=2.5,
             label=f'Erasure (χ²/DoF={erasure_chi2/erasure_dof:.3f})')
    ax1.set_xlabel('Redshift z', fontsize=12)
    ax1.set_ylabel('Distance Modulus μ', fontsize=12)
    ax1.set_title('Power-Law Erasure Theory vs ΛCDM: Full Pantheon+ Dataset', 
                  fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11, loc='upper left')
    ax1.grid(True, alpha=0.3)
    
    # Residual plots
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.scatter(z_residuals, residuals_lcdm, c='black', alpha=0.4, s=10)
    ax2.axhline(0, color='red', linestyle='--', linewidth=1.5)
    ax2.set_xlabel('Redshift z', fontsize=11)
    ax2.set_ylabel('Residual (obs - ΛCDM)', fontsize=11)
    ax2.set_title('ΛCDM Residuals', fontsize=12)
    ax2.grid(True, alpha=0.3)
    
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.scatter(z_residuals, residuals_erasure, c='magenta', alpha=0.4, s=10)
    ax3.axhline(0, color='red', linestyle='--', linewidth=1.5)
    ax3.set_xlabel('Redshift z', fontsize=11)
    ax3.set_ylabel('Residual (obs - Erasure)', fontsize=11)
    ax3.set_title('Erasure Model Residuals', fontsize=12)
    ax3.grid(True, alpha=0.3)
    
    # Residual histograms
    ax4 = fig.add_subplot(gs[2, 0])
    ax4.hist(residuals_lcdm, bins=50, color='black', alpha=0.6, edgecolor='black')
    ax4.axvline(0, color='red', linestyle='--', linewidth=1.5)
    ax4.set_xlabel('Residual', fontsize=11)
    ax4.set_ylabel('Count', fontsize=11)
    ax4.set_title(f'ΛCDM: σ = {np.std(residuals_lcdm):.3f}', fontsize=12)
    ax4.grid(True, alpha=0.3)
    
    ax5 = fig.add_subplot(gs[2, 1])
    ax5.hist(residuals_erasure, bins=50, color='magenta', alpha=0.6, edgecolor='black')
    ax5.axvline(0, color='red', linestyle='--', linewidth=1.5)
    ax5.set_xlabel('Residual', fontsize=11)
    ax5.set_ylabel('Count', fontsize=11)
    ax5.set_title(f'Erasure: σ = {np.std(residuals_erasure):.3f}', fontsize=12)
    ax5.grid(True, alpha=0.3)
    
    plt.savefig('erasure_v2_analysis.png', dpi=150, bbox_inches='tight')
    plt.show()


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    # Path to data file
    filepath = r'C:\Users\sayan\Downloads\Pantheon\Cosmos\Pantheon+SH0ES.dat'
    
    # Load data
    z_data, mu_data, mu_err = load_pantheon_data(filepath)
    
    # Fit models
    lcdm_model, lcdm_chi2, lcdm_dof = fit_lcdm(z_data, mu_data, mu_err)
    erasure_model, erasure_chi2, erasure_dof = fit_erasure(z_data, mu_data, mu_err)
    
    # Statistical comparison
    delta_chi2, delta_aic, delta_bic = model_comparison(
        lcdm_model, lcdm_chi2, lcdm_dof,
        erasure_model, erasure_chi2, erasure_dof
    )
    
    # Visualization
    plot_results(z_data, mu_data, mu_err, lcdm_model, erasure_model,
                lcdm_chi2, lcdm_dof, erasure_chi2, erasure_dof)
    
    # Final verdict
    print("\n" + "="*70)
    print("FINAL VERDICT")
    print("="*70)
    
    if delta_chi2 > 20 and delta_aic > 4:
        print("✓ The Power-Law Erasure Model is FAVORED by the data.")
        print(f"  α = {erasure_model.alpha:.2f} suggests late-time interaction.")
        print("  This warrants further investigation with CMB and BAO data.")
    elif delta_chi2 > 4:
        print("→ The models are COMPARABLE.")
        print("  Extra parameters are not strongly justified.")
        print("  Would need independent confirmation (CMB, BAO, etc.)")
    else:
        print("✗ ΛCDM is PREFERRED.")
        print("  The Power-Law Erasure Model does not improve the fit.")
        print("  The theory is falsified by current supernova data.")


if __name__ == "__main__":
    main()