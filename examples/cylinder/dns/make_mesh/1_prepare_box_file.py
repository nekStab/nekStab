#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np


def fmt_coords(numbers, format_str="{:.4f}"):
    line = ""
    lines = []
    for num in numbers:
        num_str = format_str.format(num)
        if len(line + " " + num_str) > 120:  # Fortran line length limit
            lines.append(line.strip())
            line = num_str
        else:
            line = (line + " " + num_str) if line else num_str
    if line:
        lines.append(line.strip())
    return "\n".join(lines)


def str_to_nums(seq, dz=None):
    nums = [float(x) for x in seq.strip().split() if x]
    return [x * dz if dz is not None else x for x in nums]


def positions_to_deltas(positions):
    """Convert a sequence of positions to deltas between consecutive points."""
    return np.diff(positions)


def deltas_to_positions(deltas, start_pos):
    """Convert a sequence of deltas back to positions starting from start_pos."""
    positions = np.zeros(len(deltas) + 1)
    positions[0] = start_pos
    positions[1:] = start_pos + np.cumsum(deltas)
    return positions


# Define each direction string independently for future customization
pos_in_str = "0.3 0.6 0.9 1.2 1.5 1.722 2.3225 2.5951 2.8676 3.1676 3.4685 3.7702 4.0710 4.4470 5.0794 6.0839 7.4582 9.1422 11.0125 12.9670 14.9802 17.0373 19.0940 21.1072 23.0634 25.0000"

# Allow user to customize these independently in the future
pos_up_str = pos_in_str  # For now, same as pos_in_str
pos_down_str = pos_in_str  # For now, same as pos_in_str

pos_out_str = "0.3 0.6 0.9 1.2 1.5 1.722 2.1106 2.7435 3.4167 4.1278 4.8749 5.6558 6.4684 7.3106 8.1803 9.0753 9.9936 10.9331 11.8916 12.8670 13.8573 14.8603 15.8739 16.8960 17.9245 18.9572 19.9922 21.0272 22.0611 23.1072 24.1885 25.3115 26.4738 27.6769 28.9220 30.2108 31.5446 32.9251 34.3539 35.8327 37.3633 38.9475 40.5871 42.2842 44.0406 45.8585 47.7400 49.6873 51.7028 53.7889 55.9479 58.1826 60.4954 62.8892 65.3667 67.9310 70.5850 73.3320 76.1751 79.1177 82.1632 85.3154 88.5779 96.9862 100.0"

# Parse each string independently
pos_in = str_to_nums(pos_in_str)
pos_up = np.concatenate(([0.0], str_to_nums(pos_up_str)))
pos_down = np.concatenate(([0.0], str_to_nums(pos_down_str)))
pos_out = str_to_nums(pos_out_str)

dz = 0.3

deltas_in = positions_to_deltas(pos_in) / dz
print("deltas_in = [", end="")
print(", ".join(f"{x:.4f}" for x in deltas_in), end="")
print("]")

# User parameters for geometric progressions
n_initial = 6  # Number of initial unit values
dz = 0.3  # Base delta

# Inner region parameters (center)
n_in_points = 17  # Total points desired
target_in = 25  # Target final position (not delta)

# Up (positive y) region parameters
n_up_points = 17  # Default: same as center, can be customized
target_up = target_in  # Default: same as center, can be customized

# Down (negative y) region parameters
n_down_points = 17  # Default: same as center, can be customized
target_down = target_in  # Default: same as center, can be customized

# Outer region parameters
n_out_points = 51  # Total points desired
target_out = 100  # Target final position (not delta)

# X-direction (left/right) region parameters
n_x_left_points = n_in_points  # Default: same as center, can be customized
target_x_left = target_in  # Default: same as center, can be customized

n_x_right_points = n_out_points  # Default: same as outer, can be customized
target_x_right = target_out  # Default: same as outer, can be customized


def compute_position_progression(n_points, target_pos, dz, n_initial=6):
    """Compute geometric progression for positions.
    Args:
        n_points: Total number of points desired
        target_pos: Target final position to reach
        dz: Base delta for initial points
        n_initial: Number of initial unit-delta points
    Returns:
        List of positions following geometric progression
    """
    positions = np.zeros(n_points)

    # First n_initial points must have exactly dz spacing
    for i in range(n_initial):
        positions[i] = i * dz

    # Calculate growth rate to reach target from last constant point
    n_steps = n_points - n_initial  # Steps after initial points
    start_pos = positions[n_initial - 1]  # Last constant-spaced point

    # For smooth transition, first delta after constant zone should be dz
    # Then grow geometrically to reach target
    # We want: start_pos + dz * (growth_rate^0 + growth_rate^1 + ... + growth_rate^(n_steps-1)) = target_pos
    # This is a geometric series sum: a * (1-r^n)/(1-r) where a=dz and n=n_steps
    # Solve numerically for growth_rate
    def geo_sum(r, n):
        return (1 - pow(r, n)) / (1 - r) if r != 1 else n

    def target_func(r):
        return start_pos + dz * geo_sum(r, n_steps) - target_pos

    # Binary search for growth rate (faster and more stable than root finding)
    r_min, r_max = 1.001, 2.0  # Growth rate must be > 1
    while r_max - r_min > 1e-6:
        r_mid = (r_min + r_max) / 2
        if target_func(r_mid) > 0:
            r_max = r_mid
        else:
            r_min = r_mid

    growth_rate = r_min

    # Generate remaining positions using geometric progression of deltas
    current_pos = start_pos
    current_delta = dz
    for i in range(n_initial, n_points - 1):  # Leave last point for exact integer
        current_pos += current_delta
        current_delta *= growth_rate
        positions[i] = round(current_pos, 4)

    # Force last point to be integer
    positions[-1] = round(target_pos)

    return positions, growth_rate


# Compute inner progression
positions_in, growth_rate_in = compute_position_progression(n_in_points, target_in, dz, n_initial)
deltas_in = np.diff(positions_in)

print("\nInner region progression:")
print(f"Points: {len(positions_in)}, Target position: {target_in}")
print(f"First {n_initial} points (constant dz={dz}):", end=" ")
print(", ".join(f"{x:.4f}" for x in positions_in[:n_initial]))
print(f"Growth rate: {growth_rate_in:.4f}")
print("Deltas after constant zone:", end=" ")
current = dz
for _ in range(5):
    print(f"{current:.4f}", end=" ")
    current *= growth_rate_in
print("\npositions_in = [", end="")
print(", ".join(f"{x:.4f}" for x in positions_in), end="")
print("]")

# Compute outer progression
positions_out, growth_rate_out = compute_position_progression(n_out_points, target_out, dz, n_initial)
deltas_out = np.diff(positions_out)

print("\nOuter region progression:")
print(f"Points: {len(positions_out)}, Target position: {target_out}")
print(f"First {n_initial} points (constant dz={dz}):", end=" ")
print(", ".join(f"{x:.4f}" for x in positions_out[:n_initial]))
print(f"Growth rate: {growth_rate_out:.4f}")
print("Deltas after constant zone:", end=" ")
current = dz
for _ in range(5):
    print(f"{current:.4f}", end=" ")
    current *= growth_rate_out
print("\npositions_out = [", end="")
print(", ".join(f"{x:.4f}" for x in positions_out), end="")
print("]")

# Compute up progression
positions_up, growth_rate_up = compute_position_progression(n_up_points, target_up, dz, n_initial)

# Compute down progression
positions_down, growth_rate_down = compute_position_progression(n_down_points, target_down, dz, n_initial)

# Create the constant region with reversed part
x_coords = np.concatenate([-positions_in[::-1][:-1], positions_out])  # x
print(f"\nx_coords: {fmt_coords(x_coords)}")
y_coords = np.concatenate([-positions_down[::-1][:-1], positions_up])  # y
print(f"\ny_coords: {fmt_coords(y_coords)}")


nx = len(x_coords) - 1  # number of elements is number of points minus 1
ny = len(y_coords) - 1  # number of elements is number of points minus 1

output_file = "mybox.box"
box_content = """import.rea
2    spatial dimension
1    number of fields
Box
{nx} {ny}                        nelx,nely,nelz for Box)
{x} 
{y}
v  ,O  ,   ,   ,          bc's""".format(nx=nx, ny=ny, x=fmt_coords(x_coords), y=fmt_coords(y_coords))

with open(output_file, "w") as f:
    f.write(box_content)
print(f"\nSaved to {output_file} with {nx * ny} elements")

# Plotting
fig, (ax1, ax3) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

# Plot for x-coordinates
ax1.set_ylabel("x coordinates", color="b")
ax1.tick_params(axis="y", labelcolor="b")
ax2 = ax1.twinx()

x_indices = np.arange(len(x_coords))
ax1.plot(x_indices, x_coords, "b-", label="x coordinates", markersize=4)
ax2.plot(x_indices[1:], np.diff(x_coords), "r-o", mfc="none", label="x differences", markersize=4)
ax1.set_ylim(min(x_coords), max(x_coords))
ax2.set_ylabel("x differences", color="r")
ax2.tick_params(axis="y", labelcolor="r")
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left")
ax1.set_title("X and Y Coordinates and Differences")

# Plot for y-coordinates
ax3.set_xlabel("Index")
ax3.set_ylabel("y coordinates", color="g")
ax3.tick_params(axis="y", labelcolor="g")
ax4 = ax3.twinx()

y_indices = np.arange(len(y_coords))
ax3.plot(y_indices, y_coords, "g-", label="y coordinates", markersize=4)
ax4.plot(y_indices[1:], np.diff(y_coords), "m-o", mfc="none", label="y differences", markersize=4)
ax3.set_ylim(min(y_coords), max(y_coords))
ax4.set_ylabel("y differences", color="m")
ax4.tick_params(axis="y", labelcolor="m")
lines3, labels3 = ax3.get_legend_handles_labels()
lines4, labels4 = ax4.get_legend_handles_labels()
ax3.legend(lines3 + lines4, labels3 + labels4, loc="upper left")

fig.tight_layout()
plt.savefig("x_coords_analysis.png")
plt.close()
