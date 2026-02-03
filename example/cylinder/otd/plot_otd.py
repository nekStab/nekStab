import matplotlib.pyplot as plt

# Configuration
log_files = {
    "otd_Le.dat": {"linestyle": "-", "linewidth": 1.5},
    "otd_Lr.dat": {"linestyle": "--", "linewidth": 1.5},
    "otd_Ls.dat": {"linestyle": "-.", "linewidth": 1.0},
}

# Initialize data storage
data = {}

for lf in log_files.keys():
    data[lf] = {'t': [], 'y': []}
    
    try:
        with open(lf, 'r') as file:
            first_line = True
            for line in file:
                parts = line.strip().split()
                if not parts:
                    continue
                
                t = float(parts[0])
                
                if first_line:
                    n_cols = len(parts) - 1
                    data[lf]['y'] = [[] for _ in range(n_cols)]
                    first_line = False
                
                data[lf]['t'].append(t)
                for i in range(n_cols):
                    data[lf]['y'][i].append(float(parts[i+1]))
                        
    except FileNotFoundError:
        print(f"File {lf} not found. Skipping.")
        continue

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

plt.figure(figsize=(10, 6))
for lf, d in data.items():
    t = d['t']
    for idx, y in enumerate(d['y']):
        if lf == "otd_Ls.dat" and idx == 0:

            label = rf"$\sigma_{idx+1}$ - {lf}"
            linewidth = log_files[lf]["linewidth"]
            if idx == 0:
                linewidth = 3.0
            plt.plot(t, y, label=label, 
                     linestyle=log_files[lf]["linestyle"], 
                     linewidth=linewidth, 
                     color=colors[idx % len(colors)])
            
        if lf == "otd_Lr.dat":
            label = rf"$\Re(\lambda_{idx+1})$ - {lf}"
            plt.plot(t, y, label=label, 
                     linestyle=log_files[lf]["linestyle"], 
                     linewidth=log_files[lf]["linewidth"], 
                     color=colors[idx % len(colors)])
        
        if lf == "otd_Le.dat":
            label = rf"$\lambda_{idx+1}$ - {lf}"
            plt.plot(t, y, label=label, 
                     linestyle=log_files[lf]["linestyle"], 
                     linewidth=log_files[lf]["linewidth"], 
                     color=colors[idx % len(colors)])

# plt.axhline(y=-0.2928032E-01, color='k', linestyle='--', linewidth=0.5, label=r"$\sigma$ at $Re=40$")
plt.xlabel("$t$")
plt.ylabel(r"$\lambda$, $\Re$, $\sigma$")
plt.yscale("symlog")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("otd.png", dpi=600)
plt.close()