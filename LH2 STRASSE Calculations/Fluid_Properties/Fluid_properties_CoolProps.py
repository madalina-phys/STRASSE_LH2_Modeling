import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from scipy.optimize import curve_fit
import warnings

# Suppress potential warnings from CoolProp for out-of-range values at the edges
warnings.filterwarnings("ignore", category=UserWarning)

# Define the temperature and pressure ranges
T_range = np.linspace(13.9, 30, 100)
P_range = np.linspace(7500, 140000, 100)
fluid = 'ParaHydrogen'

def plot_and_fit(x_data, y_data, title, xlabel, ylabel, fit_type='poly', poly_deg=3):
    """
    Generates a plot, fits the data, and displays the fit parameters.
    Saves the plot to a file.
    """
    plt.figure(figsize=(12, 7))
    plt.plot(x_data, y_data, 'bo', label='CoolProp Data', markersize=4)

    params = None
    equation = "Fit failed"

    try:
        # Fit the data
        if fit_type == 'poly':
            # Use numpy's polynomial fitting
            params = np.polyfit(x_data, y_data, poly_deg)
            fit_func = np.poly1d(params)
            y_fit = fit_func(x_data)
            
            # Create equation string for the plot
            equation = 'y = '
            for i, p in enumerate(params):
                power = poly_deg - i
                if abs(p) > 1e-25: # Only include significant terms
                    if i > 0 and p > 0:
                        equation += ' + '
                    elif p < 0:
                        equation += ' - ' if i > 0 else '-'
                    
                    equation += f'{abs(p):.4g}'
                    if power > 0:
                        equation += f'$\\cdot x^{power}$' if power > 1 else '$\\cdot x$'
            
        elif fit_type == 'log':
            # Use scipy's curve_fit for non-polynomial functions
            def log_func(x, a, b):
                return a * np.log(x) + b
            
            # Provide an initial guess to help the solver
            initial_guess = [1.0, np.mean(y_data)]
            params, _ = curve_fit(log_func, x_data, y_data, p0=initial_guess)
            a, b = params
            y_fit = log_poly_func(x_data, a, b)
            sign = '+' if b > 0 else '-'
            equation = f'y = {a:.4g} $\\cdot$ ln(x) {sign} {abs(b):.4g}'
        
        elif fit_type == 'log_poly':
            def log_poly_func(x, a, b, c):
                return a * np.log(x)*np.log(x) + b * np.log(x) + c
                
            # Provide an initial guess to help the solver
            initial_guess = [0.01,1.0, np.mean(y_data)]
            params, _ = curve_fit(log_poly_func, x_data, y_data, p0=initial_guess)
            a, b, c = params
            y_fit = log_poly_func(x_data, a, b, c)
            sign_b = '+' if b > 0 else '-'
            sign_c = '+' if c > 0 else '-'
            equation = f'y = {a:.4g} $\\cdot$ ln(x)^2 {sign_b} {abs(b):.4g} $\\cdot$ ln(x) {sign_c} {abs(c):.4g}'
                
        plt.plot(x_data, y_fit, 'r-', label=f'Fitted Curve')

    except Exception as e:
        print(f"Could not fit data for '{title}': {e}")
        y_fit = y_data # Plot data without fit line if fitting fails

    # Print parameters to terminal
    print(f"--- {title} ---")
    if params is not None:
        print(f"Fit Parameters for {fit_type} (from highest order term):", params)
    
    # Add details to plot
    plt.title(title, fontsize=16)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    # Place equation text box on the plot
    plt.text(0.05, 0.95, equation, transform=plt.gca().transAxes, fontsize=11,
             verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.7))
    
    # Save the plot to a file
    filename = f"{title.replace(' ', '_').replace('/', '_').lower()}.png"
    plt.savefig(filename, dpi=150)
    plt.close() # Close the figure to free up memory
    print(f"Plot saved as {filename}\n")

# --- Generate Plots vs. Temperature ---
print("Generating plots as a function of Temperature...")
# 1. Density of liquid vs. Temperature
rho_liq_T = [PropsSI('D', 'T', T, 'Q', 0, fluid) for T in T_range]
plot_and_fit(T_range, rho_liq_T, 'Liquid Hydrogen Density vs. Temperature', 'Temperature (K)', 'Density (kg/m$^3$)', poly_deg=3)

# 2. Density of vapor vs. Temperature
rho_vap_T = [PropsSI('D', 'T', T, 'Q', 1, fluid) for T in T_range]
plot_and_fit(T_range, rho_vap_T, 'Vapor Hydrogen Density vs. Temperature', 'Temperature (K)', 'Density (kg/m$^3$)', poly_deg=4)

# 3. Viscosity of liquid vs. Temperature
mu_liq_T = [PropsSI('V', 'T', T, 'Q', 0, fluid) for T in T_range]
plot_and_fit(T_range, mu_liq_T, 'Liquid Hydrogen Viscosity vs. Temperature', 'Temperature (K)', 'Viscosity (Pa·s)', poly_deg=3)

# 4. Viscosity of vapor vs. Temperature
mu_vap_T = [PropsSI('V', 'T', T, 'Q', 1, fluid) for T in T_range]
plot_and_fit(T_range, mu_vap_T, 'Vapor Hydrogen Viscosity vs. Temperature', 'Temperature (K)', 'Viscosity (Pa·s)', poly_deg=4)

# 5. Specific heat of liquid vs. Temperature
cp_liq_T = [PropsSI('C', 'T', T, 'Q', 0, fluid) for T in T_range]
plot_and_fit(T_range, cp_liq_T, 'Liquid Hydrogen Specific Heat vs. Temperature', 'Temperature (K)', 'Specific Heat (J/kg·K)', poly_deg=5)

# 6. Specific heat of vapor vs. Temperature
cp_vap_T = [PropsSI('C', 'T', T, 'Q', 1, fluid) for T in T_range]
plot_and_fit(T_range, cp_vap_T, 'Vapor Hydrogen Specific Heat vs. Temperature', 'Temperature (K)', 'Specific Heat (J/kg·K)', poly_deg=5)

# 7. Thermal conductivity of liquid vs. Temperature
k_liq_T = [PropsSI('L', 'T', T, 'Q', 0, fluid) for T in T_range]
plot_and_fit(T_range, k_liq_T, 'Liquid Hydrogen Thermal Conductivity vs. Temperature', 'Temperature (K)', 'Thermal Conductivity (W/m·K)', poly_deg=3)

# 8. Thermal conductivity of vapor vs. Temperature
k_vap_T = [PropsSI('L', 'T', T, 'Q', 1, fluid) for T in T_range]
plot_and_fit(T_range, k_vap_T, 'Vapor Hydrogen Thermal Conductivity vs. Temperature', 'Temperature (K)', 'Thermal Conductivity (W/m·K)', poly_deg=4)

# 9. Saturation pressure vs. Temperature
P_sat_T = [PropsSI('P', 'T', T, 'Q', 0, fluid) for T in T_range]
plot_and_fit(T_range, P_sat_T, 'Hydrogen Saturation Pressure vs. Temperature', 'Temperature (K)', 'Pressure (Pa)', poly_deg=5)


# --- Generate Plots vs. Pressure ---
print("Generating plots as a function of Pressure...")
# 10. Saturation temperature vs. Pressure
T_sat_P = [PropsSI('T', 'P', P, 'Q', 0, fluid) for P in P_range]
plot_and_fit(P_range, T_sat_P, 'Hydrogen Saturation Temperature vs. Pressure', 'Pressure (Pa)', 'Temperature (K)', fit_type='log_poly')

# 11. Latent vaporization heat vs. Pressure
h_vap = np.array([PropsSI('H', 'P', P, 'Q', 1, fluid) for P in P_range])
h_liq = np.array([PropsSI('H', 'P', P, 'Q', 0, fluid) for P in P_range])
L_vap = h_vap - h_liq
plot_and_fit(P_range, L_vap, 'Hydrogen Latent Heat of Vaporization vs. Pressure', 'Pressure (Pa)', 'Latent Heat (J/kg)', poly_deg=5)

# 12. Density of liquid at saturation vs. Pressure
rho_liq_P = [PropsSI('D', 'P', P, 'Q', 0, fluid) for P in P_range]
plot_and_fit(P_range, rho_liq_P, 'Liquid Hydrogen Density vs. Saturation Pressure', 'Pressure (Pa)', 'Density (kg/m$^3$)', poly_deg=4)

# 13. Density of vapor at saturation vs. Pressure
rho_vap_P = [PropsSI('D', 'P', P, 'Q', 1, fluid) for P in P_range]
plot_and_fit(P_range, rho_vap_P, 'Vapor Hydrogen Density vs. Saturation Pressure', 'Pressure (Pa)', 'Density (kg/m$^3$)', poly_deg=2)

print("\nAll plots and fits have been generated and saved as PNG files.")

