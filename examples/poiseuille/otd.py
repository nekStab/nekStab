import matplotlib.pyplot as plt

# Configuration
log_files = {
    "Reflog": {"linestyle": ":", "linewidth": 2.5},
    "logfile": {"linestyle": "-", "linewidth": 1.5},
}

# Adjusted number of y-columns based on log file structure
n_lrr = 2  # For 'Lr | Re': columns 7 and 8
n_ftle = 2  # For 'FTLE PRD': columns 8 and 9 (adjust based on actual data)

# Initialize data storage
lrr_data = {}
ftle_data = {}

for lf in log_files.keys():
    lrr_data[lf] = {'t': [], 'y': [[] for _ in range(n_lrr)]}
    ftle_data[lf] = {'t': [], 'y': [[] for _ in range(n_ftle)]}
    
    try:
        with open(lf, 'r') as file:
            for line in file:
                if 'Lr | Re' in line:
                    parts = line.strip().split()
                    t = float(parts[2])  # time
                    lrr_data[lf]['t'].append(t)
                    for i in range(n_lrr):
                        y_val = float(parts[6 + i])  # Columns 7 and 8
                        lrr_data[lf]['y'][i].append(y_val)
                elif 'FTLE PRD' in line:
                    parts = line.strip().split()
                    t = float(parts[6])  # time 
                    ftle_data[lf]['t'].append(t)
                    for i in range(n_ftle):
                        y_val = float(parts[7 + i])  # Columns 8 and 9
                        ftle_data[lf]['y'][i].append(y_val)
    except FileNotFoundError:
        print(f"File {lf} not found. Skipping.")
        continue

# Define colors for each eigenvalue
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

# Plotting 'Lr | Re'
plt.figure(figsize=(10, 6))
for lf, data in lrr_data.items():
    t = data['t']
    for idx, y in enumerate(data['y']):
        plt.plot(t, y, label=f"Lr_Re_{idx+1} - {lf}", linestyle=log_files[lf]["linestyle"], linewidth=log_files[lf]["linewidth"], color=colors[idx % len(colors)])
plt.xlabel("$t$")
plt.ylabel("Re(Eigenvalues)")
plt.title("Re Part of Eigenvalues")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plot_Lr_Re.png", dpi=600)
plt.close()
print("Saved plot as 'plot_Lr_Re.png'")

# Plotting 'FTLE PRD'
plt.figure(figsize=(10, 6))
for lf, data in ftle_data.items():
    t = data['t']
    for idx, y in enumerate(data['y']):
        plt.plot(t, y, label=f"FTLE_PRD_{idx+1} - {lf}", linestyle=log_files[lf]["linestyle"], linewidth=log_files[lf]["linewidth"], color=colors[idx % len(colors)])
plt.xlabel("$t$")
plt.ylabel("Lambda Values")
plt.title("Lambda Values")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plot_FTLE_PRD.png", dpi=600)
plt.close()
print("Saved plot as 'plot_FTLE_PRD.png'")