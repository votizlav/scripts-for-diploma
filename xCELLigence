import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
from numpy import trapz

df = pd.read_csv('/home/bmo/Downloads/xcell_data.csv', delimiter=',', decimal=',')
df.columns = ['Time', 'Time_str', 'Ctrl_Y', 'Ctrl_SD', 'C0.5_Y', 'C0.5_SD',
              'C1_Y', 'C1_SD', 'C6_Y', 'C6_SD', 'C10_Y', 'C10_SD']
colors = {
    'Control': 'red',
    'Cis_0.5': 'green',
    'Cis_1': 'blue',
    'Cis_6': 'magenta',
    'Cis_10': 'cyan'
}

# Подписи для легенды
labels = {
    'Control': '0 µg/ml (Control)',
    'Cis_0.5': '0.5 µg/ml',
    'Cis_1': '1 µg/ml',
    'Cis_6': '6 µg/ml',
    'Cis_10': '10 µg/ml'
}

time = df['Time']
groups = {
    'Control': df['Ctrl_Y'],
    'Cis_0.5': df['C0.5_Y'],
    'Cis_1': df['C1_Y'],
    'Cis_6': df['C6_Y'],
    'Cis_10': df['C10_Y']
}
plt.figure(figsize=(14, 8))
for key, y in groups.items():
    y_smooth = savgol_filter(y, window_length=11, polyorder=3)
    plt.plot(time, y, alpha=0.5, linestyle='--', color=colors[key], label=f"{labels[key]} raw")
    plt.plot(time, y_smooth, color=colors[key], label=f"{labels[key]} smoothed")
plt.rcParams.update({'font.size': 18})
plt.title('Comparison of the effect of cisplatin at different dosages on the CI depending on time.')
plt.xlabel("Time (hours)", fontsize = 18)
plt.ylabel("Normalized Cell Index", fontsize = 18)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('/run/media/bmo/пам-ПАМ/биофак/diploma/xCellingence_stat/pyplot.png', dpi=500)
plt.show()
groups_1 = {
    'Control': ('Ctrl_Y', 'Ctrl_SD'),
    'Cis_0.5': ('C0.5_Y', 'C0.5_SD'),
    'Cis_1': ('C1_Y', 'C1_SD'),
    'Cis_6': ('C6_Y', 'C6_SD'),
    'Cis_10': ('C10_Y', 'C10_SD')
}

auc_results = {}
auc_std = {}

for group, (y_col, sd_col) in groups.items():
    y = df[y_col]
    sd = df[sd_col]
    time = df['Time']

    auc = trapz(y, time)
    auc_results[group] = auc

    # Расчёт SD AUC через ошибку интегрирования:
    auc_sd = np.sqrt(np.sum((sd * np.gradient(time))**2))
    auc_std[group] = auc_sd

# === Построение графика ===
group_labels = list(auc_results.keys())
auc_vals = list(auc_results.values())
auc_errors = [auc_std[g] for g in group_labels]

plt.figure(figsize=(10, 7))
bars = plt.bar(group_labels, auc_vals, yerr=auc_errors, capsize=10, color='skyblue', edgecolor='black')

# Добавляем значения AUC ± SD на столбцы
for bar, val, err in zip(bars, auc_vals, auc_errors):
    height = bar.get_height()
    label = f'{val:.2f} ± {err:.2f}'
    plt.text(bar.get_x() + bar.get_width()/2., 
             height + 0.005 * max(auc_vals),     
             label,
             ha='center',
             va='bottom',
             fontsize=14)

plt.rcParams.update({'font.size': 16})
plt.ylabel('AUC (Total Proliferation)')
plt.xlabel("Groups")
plt.title('Proliferation Suppression by Cisplatin')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('/run/media/bmo/пам-ПАМ/биофак/diploma/xCellingence_stat/barplot_AUC_SD.png', dpi=500)
plt.show()
plt.figure(figsize=(12, 7))
for key, y in groups.items():
    y_smooth = savgol_filter(y, window_length=11, polyorder=3)
    dy_dt = np.gradient(y_smooth, time)
    plt.plot(time, dy_dt, label=f"{labels[key]} d(CI)/dt", color=colors[key])
plt.rcParams.update({'font.size': 16})
plt.title("Growth rate (d(CI)/dt) over time")
plt.xlabel("Time (hours)")
plt.ylabel("Growth rate")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('/run/media/bmo/пам-ПАМ/биофак/diploma/xCellingence_stat/pyplot_kinan.png', dpi=500)
plt.show()
def kinetic_analysis(time, ci):
    # Сглаживание данных
    ci_smooth = savgol_filter(ci, window_length=11, polyorder=3)

    # Производная — скорость роста
    dci_dt = np.gradient(ci_smooth, time)

    # Найдём max slope и когда она была
    max_slope = np.max(dci_dt)
    t_max_slope = time[np.argmax(dci_dt)]

    # Максимальное значение CI
    ci_max = np.max(ci)
    t_ci_max = time[np.argmax(ci)]

    # Найдём лаг-фазу (время, когда рост CI стал выше 10% от максимума)
    threshold = 0.1 * ci_max
    lag_indices = np.where(ci >= threshold)[0]
    lag_time = time[lag_indices[0]] if len(lag_indices) > 0 else np.nan

    return {
        'LagTime': lag_time,
        'MaxSlope': max_slope,
        'T_MaxSlope': t_max_slope,
        'CI_Max': ci_max,
        'T_CI_Max': t_ci_max
    }
results = {}
results['Control'] = kinetic_analysis(df['Time'], df['Ctrl_Y'])
results['Cis_0.5'] = kinetic_analysis(df['Time'], df['C0.5_Y'])
results['Cis_1'] = kinetic_analysis(df['Time'], df['C1_Y'])
results['Cis_6'] = kinetic_analysis(df['Time'], df['C6_Y'])
results['Cis_10'] = kinetic_analysis(df['Time'], df['C10_Y'])

for group, res in results.items():
    print(f"{group}:")
    for k, v in res.items():
        print(f"  {k}: {v:.2f}")
kinetic_data = []
for key, y in groups.items():
    y_smooth = savgol_filter(y, window_length=11, polyorder=3)
    dy_dt = np.gradient(y_smooth, time)
    
    # Tmax — время максимальной скорости роста
    max_slope_idx = np.argmax(dy_dt)
    Tmax = time[max_slope_idx]
    µmax = dy_dt[max_slope_idx]

    # Lag phase — первое значение, где скорость > 5% от µmax
    lag_threshold = 0.05 * µmax
    lag_indices = np.where(dy_dt > lag_threshold)[0]
    Lag = time[lag_indices[0]] if len(lag_indices) > 0 else np.nan

    kinetic_data.append({
        'Group': labels[key],
        'Tmax (h)': round(Tmax, 2),
        'µmax': round(µmax, 4),
        'Lag phase (h)': round(Lag, 2) if not np.isnan(Lag) else np.nan
    })

kinetic_df = pd.DataFrame(kinetic_data)
print(kinetic_df)
kinetic_df.to_csv("/run/media/bmo/пам-ПАМ/биофак/diploma/xCellingence_stat/kinetic_parameters.csv", index=False)

# Heatmap
plt.figure(figsize=(10, 7))
sns.heatmap(
    kinetic_df.set_index('Group'),
    annot=True, cmap='YlGnBu', fmt=".2f"
)
plt.rcParams.update({'font.size': 18})
plt.title("Kinetic parameters heatmap")
plt.tight_layout()
plt.savefig('/run/media/bmo/пам-ПАМ/биофак/diploma/xCellingence_stat/kinetic_heatmap.png', dpi=500)
plt.show()
