"""
Wavenumber triplet generator for Free-Stream-Turbulence.

Distributes (omega, gamma, beta) triplets on spherical shells using
dodecahedron vertices (randomly rotated and mirrored).

Ported from MATLAB: example/slot_FST/preprocessFST/wavenumber.m
"""

import numpy as np


# Golden ratio constant
PHI = (1.0 + np.sqrt(5.0)) / 2.0
RADIUS0 = np.sqrt(3.0)

# Dodecahedron vertices (3, 20) - copied exactly from wavenumber.m lines 15-95
# Columns correspond to the original vec_c(1:3, 1:20) in MATLAB
DODECAHEDRON = np.array([
    # Vertices 1-8: cube corners (±1, ±1, ±1)
    [1.0,    1.0,    1.0,    1.0,   -1.0,   -1.0,   -1.0,   -1.0,
    # Vertices 9-12: (0, ±1/phi, ±phi)
     0.0,    0.0,    0.0,    0.0,
    # Vertices 13-16: (±phi, 0, ±1/phi)
     PHI,    PHI,   -PHI,   -PHI,
    # Vertices 17-20: (±1/phi, ±phi, 0)
     1.0/PHI,  1.0/PHI, -1.0/PHI, -1.0/PHI],
    # Row 2 (Gamma = y)
    [1.0,    1.0,   -1.0,   -1.0,    1.0,    1.0,   -1.0,   -1.0,
     1.0/PHI, -1.0/PHI, 1.0/PHI, -1.0/PHI,
     0.0,    0.0,    0.0,    0.0,
     PHI,   -PHI,   -PHI,    PHI],
    # Row 3 (Beta = z)
    [1.0,   -1.0,    1.0,   -1.0,    1.0,   -1.0,    1.0,   -1.0,
     PHI,    PHI,   -PHI,   -PHI,
     1.0/PHI, -1.0/PHI, 1.0/PHI, -1.0/PHI,
     0.0,    0.0,    0.0,    0.0]
], dtype=np.float64)


def cart2sph(x, y, z):
    """
    Convert Cartesian to spherical coordinates (MATLAB convention).
    
    Returns:
        az : azimuth angle (atan2(y, x))
        elev : elevation angle (atan2(z, hypot(x, y)))
        r : radius (sqrt(x^2 + y^2 + z^2))
    """
    az = np.arctan2(y, x)
    elev = np.arctan2(z, np.hypot(x, y))
    r = np.sqrt(x**2 + y**2 + z**2)
    return az, elev, r


def sph2cart(az, elev, r):
    """
    Convert spherical to Cartesian coordinates (MATLAB convention).
    
    Returns:
        x : r*cos(elev)*cos(az)
        y : r*cos(elev)*sin(az)
        z : r*sin(elev)
    """
    cos_elev = np.cos(elev)
    x = r * cos_elev * np.cos(az)
    y = r * cos_elev * np.sin(az)
    z = r * np.sin(elev)
    return x, y, z


def cart2pol(x, y):
    """
    Convert Cartesian to polar coordinates in 2D.
    
    Returns:
        theta : atan2(y, x)
        rho : hypot(x, y)
    """
    theta = np.arctan2(y, x)
    rho = np.hypot(x, y)
    return theta, rho


def pol2cart(theta, rho):
    """
    Convert polar to Cartesian coordinates in 2D.
    
    Returns:
        x : rho*cos(theta)
        y : rho*sin(theta)
    """
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    return x, y


def generate_wavenumbers(numk, kkini, kkfin, seed=None):
    """
    Generate wavenumber triplets (omega, gamma, beta) distributed on spherical shells.
    
    Uses rejection sampling to ensure exactly 10 valid triplets per shell
    (those with omega > 0 AND gamma > 0).
    
    Parameters:
        numk : int
            Number of wavenumber shells
        kkini : float
            Minimum wavenumber radius
        kkfin : float
            Maximum wavenumber radius
        seed : int or None
            Random seed for reproducibility (passed to np.random.default_rng)
    
    Returns:
        wavenumbers : ndarray of shape (numk*10, 3)
            Each row is (omega, gamma, beta) triplet.
            Triplets are ordered by shell: shell 0's 10, then shell 1's 10, etc.
    """
    rng = np.random.default_rng(seed)
    
    kk = np.linspace(kkini, kkfin, numk)
    wavenumbers = []
    
    for i in range(numk):
        ok = False
        
        while not ok:
            # Draw 2 random uniforms in [0, 1)
            random_vals = rng.uniform(0.0, 1.0, size=2)
            
            # Convert dodecahedron vertices to spherical coordinates
            az, elev, r = cart2sph(DODECAHEDRON[0, :], 
                                    DODECAHEDRON[1, :], 
                                    DODECAHEDRON[2, :])
            
            # Scale radius to kk(i)
            r_scaled = (r / RADIUS0) * kk[i]
            
            # Convert back to Cartesian
            x, y, z = sph2cart(az, elev, r_scaled)
            
            # Rotation 1: about z-axis (beta direction)
            # Apply rotation using polar coords in (omega, gamma) plane
            theta_xy, rho_xy = cart2pol(x, y)
            theta_xy += 2.0 * np.pi * random_vals[0]
            x, y = pol2cart(theta_xy, rho_xy)
            
            # Rotation 2: about x-axis (omega direction)
            # Apply rotation using polar coords in (gamma, beta) plane
            theta_yz, rho_yz = cart2pol(y, z)
            theta_yz += np.pi * random_vals[1]
            y, z = pol2cart(theta_yz, rho_yz)
            
            # Mirror about gamma: append vertices with gamma -> -gamma
            x_mirrored = np.concatenate([x, x])
            y_mirrored = np.concatenate([y, -y])
            z_mirrored = np.concatenate([z, z])
            
            # Count points with omega > 0 AND gamma > 0
            valid_mask = (x_mirrored > 0) & (y_mirrored > 0)
            valid_count = np.sum(valid_mask)
            
            # Accept if exactly 10 valid points (i.e., 40/4)
            if valid_count == 10:
                ok = True
                # Extract the valid triplets
                valid_x = x_mirrored[valid_mask]
                valid_y = y_mirrored[valid_mask]
                valid_z = z_mirrored[valid_mask]
                
                # Stack into (omega, gamma, beta) columns
                shell_triplets = np.column_stack([valid_x, valid_y, valid_z])
                wavenumbers.append(shell_triplets)
    
    # Concatenate all shells
    result = np.vstack(wavenumbers)
    return result


if __name__ == "__main__":
    # Self-test
    print("=" * 70)
    print("Testing wavenumber.py port from MATLAB")
    print("=" * 70)
    
    # Compile check
    print("\n[1/5] Checking py_compile...")
    import py_compile
    try:
        py_compile.compile(__file__, doraise=True)
        print("      PASS: File compiles without syntax errors")
    except py_compile.PyCompileError as e:
        print(f"      FAIL: {e}")
        exit(1)
    
    # Generate wavenumbers
    print("\n[2/5] Generating wavenumbers (numk=5, kkini=0.5, kkfin=5.0, seed=42)...")
    w = generate_wavenumbers(numk=5, kkini=0.5, kkfin=5.0, seed=42)
    
    # Check shape
    print(f"\n[3/5] Verifying shape...")
    print(f"      Shape: {w.shape}")
    if w.shape != (50, 3):
        print(f"      FAIL: Expected (50, 3), got {w.shape}")
        exit(1)
    print(f"      PASS: Shape is (50, 3)")
    
    # Check positivity
    print(f"\n[4/5] Verifying omega > 0 and gamma > 0...")
    omega = w[:, 0]
    gamma = w[:, 1]
    beta = w[:, 2]
    
    min_omega = np.min(omega)
    min_gamma = np.min(gamma)
    
    print(f"      min(omega) = {min_omega:.15e}")
    print(f"      min(gamma) = {min_gamma:.15e}")
    
    if min_omega <= 0 or min_gamma <= 0:
        print(f"      FAIL: Found non-positive values")
        exit(1)
    print(f"      PASS: All omega > 0 and all gamma > 0")
    
    # Check per-shell radii
    print(f"\n[5/5] Verifying shell radii...")
    kk = np.linspace(0.5, 5.0, 5)
    print(f"      Shell targets: {kk}")
    
    max_dev_overall = 0.0
    for j in range(5):
        shell_points = w[j*10:(j+1)*10, :]
        radii = np.sqrt(shell_points[:, 0]**2 + shell_points[:, 1]**2 + shell_points[:, 2]**2)
        max_dev = np.max(np.abs(radii - kk[j]))
        print(f"      Shell {j}: kk={kk[j]:.1f}, max|r - kk| = {max_dev:.15e}")
        max_dev_overall = max(max_dev_overall, max_dev)
    
    if max_dev_overall >= 1e-9:
        print(f"      FAIL: Radius deviation exceeds 1e-9")
        exit(1)
    print(f"      PASS: All radius deviations < 1e-9")
    
    print("\n" + "=" * 70)
    print("ALL TESTS PASSED")
    print("=" * 70)
