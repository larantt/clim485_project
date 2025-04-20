import numpy as np
import matplotlib.pyplot as plt


# Optical depths
tau_total = 0.027357 + ((0.000023 * 14) * 3)
tau_cld = tau_total + (0.01376614 * 3)

# Altitude from 0 to 20 km
z = np.linspace(0, 20, 500)

# Split altitude into two parts
z_upper = z[z >= 10]   # 10 km to 20 km
z_lower = z[z < 10]    # 0 km to 10 km

# Upper layer: tau from 0 to halfway based on tau_total
tau_upper = tau_total * (1 - z_upper / 20)

# Lower layer: tau increases from tau_upper[-1] to tau_cld linearly over 10 km
tau_lower = np.linspace(tau_cld, tau_upper[0], len(z_lower))

tau_z = tau_total * (1 - z / 20)
# Combine full profile
tau_cloud = np.concatenate([tau_lower, tau_upper])
z_full = np.concatenate([z_lower, z_upper])



t_z = np.exp(-tau_z)
t_cloud = np.exp(-tau_cloud)

# Weighting function: negative derivative of transmittance w.r.t. altitude
WF_z = np.gradient(t_z, z)
WF_cloud = np.gradient(t_cloud, z)

def tau():

    # Print total optical depth
    print(f"Total Optical Depth: {tau_total:.6f}")

    # Plot it
    plt.figure(figsize=(5, 8))
    plt.plot(tau_z, z, label='Clear τ (tau_total)', color='orange')
    plt.plot(tau_cloud, z_full, label='Cloudy τ (Cloud starts at 10 km)', color='blue', linewidth=2)
    plt.xlabel('Optical Depth (τ)')
    plt.ylabel('Altitude (km)')
    plt.title('Hybrid Optical Depth Profile\nCloud Impact Below 10 km')
    plt.grid(True)
    plt.legend()
    plt.show()  

def t_cld():

    # Print total optical depth
    print(f"Total Transmission")

    # Plot
    plt.figure(figsize=(5, 8))
    plt.plot(t_z, z, color='darkorange', label='Clear')
    plt.plot(t_cloud, z, color='green', label='Cloud')
    plt.axhline(0, color='gray', linestyle='--', linewidth=1)
    plt.xlim(0.8,1)
    plt.xlabel('Optical Depth (τ)')
    plt.ylabel('Normalized Altitude')
    plt.title('Optical Depth Profile\nHomogeneous Atmospheric Layer')
    plt.grid(True)
    plt.legend()
    plt.show() 

def Weight_F():

    # Plot with zoomed x-axis
    plt.figure(figsize=(6, 8))
    plt.plot(WF_z, z, color='crimson', label='Clear')
    plt.plot(WF_cloud, z, color='blue', label='Cloud')
    plt.xlabel('Weighting Function Value')
    plt.ylabel('Altitude (km)')
    plt.title('Weighting Function Profile\nHomogeneous Atmospheric Layer (μ = 1)')
    plt.grid(True)
    plt.legend()
    plt.xlim(0, 0.1)
    # Zoom in based on actual max value
    plt.show()


