# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline, make_interp_spline

# Read original mesh points from file, excluding the last point
data = np.loadtxt("original_mesh.txt")[:-1]

# Calculate original spacings
original_spacing = np.diff(data)
x_spacing = (data[1:] + data[:-1]) / 2  # Midpoints where spacing is measured

# Find the equidistant region around x=0
zero_idx = np.where(data == 0)[0][0]
spacing_around_zero = original_spacing[zero_idx - 1 : zero_idx + 1]
mean_spacing_zero = np.mean(spacing_around_zero)
tol = 1e-6  # Tolerance for identifying equidistant points

# Find extent of equidistant region
left_idx = zero_idx
right_idx = zero_idx
while left_idx > 0 and abs(original_spacing[left_idx - 1] - mean_spacing_zero) < tol:
    left_idx -= 1
while right_idx < len(original_spacing) and abs(original_spacing[right_idx] - mean_spacing_zero) < tol:
    right_idx += 1

print(f"\nEquidistant region around x=0:")
print(f"Range: [{data[left_idx]:.3f}, {data[right_idx]:.3f}]")
print(f"Spacing: {mean_spacing_zero:.3f}")
print(f"Number of points: {right_idx - left_idx + 1}")

# Create points array
x_pos = np.arange(len(data))

# Generate smooth curves for plotting
mesh_x_smooth = np.linspace(min(data), max(data), 1000)
x_smooth_spacing = np.linspace(min(x_spacing), max(x_spacing), 1000)


# Different spline orders for mesh distribution with preservation of equidistant region
def make_preserving_spline(x, y, k):
    # Split into three regions: before, equidistant, and after
    x_before = x[: left_idx + 1]
    y_before = y[: left_idx + 1]
    x_after = x[right_idx:]
    y_after = y[right_idx:]

    # Fit splines to outer regions
    spline_before = make_interp_spline(x_before, y_before, k=k) if len(x_before) > k else None
    spline_after = make_interp_spline(x_after, y_after, k=k) if len(x_after) > k else None

    # Return a function that combines all regions
    def combined_spline(x_new):
        y_new = np.zeros_like(x_new)
        for i, x_val in enumerate(x_new):
            if x_val <= data[left_idx]:
                y_new[i] = spline_before(x_val) if spline_before else y_before[0]
            elif x_val >= data[right_idx]:
                y_new[i] = spline_after(x_val) if spline_after else y_after[-1]
            else:
                # Linear interpolation in equidistant region
                idx = np.searchsorted(data[left_idx : right_idx + 1], x_val)
                y_new[i] = y[left_idx : right_idx + 1][idx - 1] + (y[left_idx : right_idx + 1][idx] - y[left_idx : right_idx + 1][idx - 1]) * (x_val - data[left_idx : right_idx + 1][idx - 1]) / (data[left_idx : right_idx + 1][idx] - data[left_idx : right_idx + 1][idx - 1])
        return y_new

    return combined_spline


# Create splines that preserve equidistant region
spline_mesh_2 = make_preserving_spline(data, x_pos, k=2)  # Quadratic
spline_mesh_3 = make_preserving_spline(data, x_pos, k=3)  # Cubic
spline_mesh_5 = make_preserving_spline(data, x_pos, k=5)  # 5th order

mesh_fit_spline2 = spline_mesh_2(mesh_x_smooth)
mesh_fit_spline3 = spline_mesh_3(mesh_x_smooth)
mesh_fit_spline5 = spline_mesh_5(mesh_x_smooth)

# Fits for spacing distribution with different smoothing factors
spline_spacing_2 = UnivariateSpline(x_spacing, original_spacing, k=2, s=0.05)  # Light smoothing
spline_spacing_3 = UnivariateSpline(x_spacing, original_spacing, k=3, s=0.1)  # Medium smoothing
spline_spacing_5 = UnivariateSpline(x_spacing, original_spacing, k=5, s=0.2)  # More smoothing

spacing_fit_spline2 = spline_spacing_2(x_smooth_spacing)
spacing_fit_spline3 = spline_spacing_3(x_smooth_spacing)
spacing_fit_spline5 = spline_spacing_5(x_smooth_spacing)

# Calculate errors at original points
# Use relative error and add small epsilon to avoid division by zero
eps = 1e-10


def calc_relative_error(predicted, actual):
    return np.abs(predicted - actual) / (np.abs(actual) + eps)


# Calculate relative errors and smooth them slightly for visualization
window = 5  # Small window for moving average


def smooth_error(error):
    return np.convolve(error, np.ones(window) / window, mode="valid")


# Calculate smooth x points for error plot (removing edges due to smoothing)
x_error = data[window // 2 : -(window // 2)]

mesh_error_spline2 = smooth_error(calc_relative_error(spline_mesh_2(data), x_pos))
mesh_error_spline3 = smooth_error(calc_relative_error(spline_mesh_3(data), x_pos))
mesh_error_spline5 = smooth_error(calc_relative_error(spline_mesh_5(data), x_pos))

# Create figure with mesh, spacing and errors
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(15, 15), height_ratios=[2, 1, 1])

# First plot: Mesh points (with swapped axes)
ax1.plot(data, x_pos, "o", label="Original Points", alpha=0.6, color="black")
ax1.plot(mesh_x_smooth, mesh_fit_spline2, "-", label="Quadratic Spline", alpha=0.8)
ax1.plot(mesh_x_smooth, mesh_fit_spline3, "--", label="Cubic Spline", alpha=0.8)
ax1.plot(mesh_x_smooth, mesh_fit_spline5, ":", label="5th Order Spline", alpha=0.8)
ax1.axvline(x=0, color="red", linestyle="--", alpha=0.3)
ax1.set_ylabel("Point Index")
ax1.set_xlabel("Mesh Point Location (D)")
ax1.set_title("Original Mesh Distribution")
ax1.legend()
ax1.grid(True)

# Second plot: Spacing vs mesh location
ax2.plot(x_spacing, original_spacing, "o", label="Original Spacing", alpha=0.6, color="black")
ax2.plot(x_smooth_spacing, spacing_fit_spline2, "-", label="Quadratic Spline (s=0.05)", alpha=0.8)
ax2.plot(x_smooth_spacing, spacing_fit_spline3, "--", label="Cubic Spline (s=0.1)", alpha=0.8)
ax2.plot(x_smooth_spacing, spacing_fit_spline5, ":", label="5th Order Spline (s=0.2)", alpha=0.8)
ax2.axvline(x=0, color="red", linestyle="--", alpha=0.3)
ax2.set_xlabel("Mesh Point Location (D)")
ax2.set_ylabel("Grid Spacing Δ")
ax2.set_title("Grid Spacing vs Location")
ax2.legend()
ax2.grid(True)

# Third plot: Logarithmic errors
ax3.semilogy(x_error, mesh_error_spline2, "-", label="Quadratic Spline", alpha=0.8)
ax3.semilogy(x_error, mesh_error_spline3, "--", label="Cubic Spline", alpha=0.8)
ax3.semilogy(x_error, mesh_error_spline5, ":", label="5th Order Spline", alpha=0.8)
ax3.axvline(x=0, color="red", linestyle="--", alpha=0.3)
ax3.set_xlabel("Mesh Point Location (D)")
ax3.set_ylabel("Relative Error (log scale)")
ax3.set_title("Smoothed Relative Error in Distribution Fits")
ax3.legend()
ax3.grid(True)

plt.tight_layout()

# Save the figure using the original filename with png extension
fig_filename = "original_mesh.png"
plt.savefig(fig_filename, dpi=600)

# Print some statistics about the original mesh
print(f"Original domain: [{data[0]:.3f}, {data[-1]:.3f}]")
print(f"Number of points: {len(data)}")
print(f"Spacing range: [{np.min(original_spacing):.3f}, {np.max(original_spacing):.3f}]")
print(f"Mean spacing: {np.mean(original_spacing):.3f}")
print(f"Median spacing: {np.median(original_spacing):.3f}")
print(f"\nSpacing measured at x-range: [{min(x_spacing):.3f}, {max(x_spacing):.3f}]")

# Print method information and error statistics
print("\nSpline Methods Used:")
print("1. Quadratic Spline (k=2)")
print("2. Cubic Spline (k=3)")
print("3. 5th Order Spline (k=5)")

print("\nSpline Smoothing Parameters:")
print("- Quadratic: s=0.05 (light smoothing)")
print("- Cubic: s=0.1 (medium smoothing)")
print("- 5th Order: s=0.2 (more smoothing)")

print("\nMaximum Relative Errors in Distribution Fits:")
print(f"Quadratic Spline: {np.max(mesh_error_spline2):.2e}")
print(f"Cubic Spline: {np.max(mesh_error_spline3):.2e}")
print(f"5th Order Spline: {np.max(mesh_error_spline5):.2e}")
