import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from statannotations.Annotator import Annotator
from pathlib import Path

# Загрузка данных
file_path = input('Input path to your output file from qPCR: ')
output_path = input('Input path where to save results: ')
output_dir = Path(f"{output_dir}/gene_plots")
output_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(file_path, encoding='cp1251', sep=';')

# Предобработка
df['Cq'] = df['Cq'].str.replace(',', '.', regex=False)
df['Cq'] = pd.to_numeric(df['Cq'], errors='coerce')
df.dropna(subset=['Cq'], inplace=True)
df['Sample'] = df['Sample'].astype(str).str.strip()
df['Target'] = df['Target'].astype(str).str.strip()

# Параметры
reference_gene = 'GAPDH'
control_sample = 'K-'
genes = [g for g in df['Target'].unique() if g != reference_gene]

# Расчёт ΔCt
merged = []
for sample in df['Sample'].unique():
    for gene in genes:
        cq_gene = df[(df['Sample'] == sample) & (df['Target'] == gene)]['Cq'].values
        cq_ref = df[(df['Sample'] == sample) & (df['Target'] == reference_gene)]['Cq'].values
        if len(cq_gene) > 0 and len(cq_ref) > 0:
            for g, r in zip(cq_gene, cq_ref):
                merged.append({'Sample': sample, 'Gene': gene, 'DeltaCt': g - r})

delta_ct_df = pd.DataFrame(merged)

# Расчёт ΔΔCt и Fold Change
control_mean = delta_ct_df[delta_ct_df['Sample'] == control_sample].groupby('Gene')['DeltaCt'].mean()
delta_ct_df = delta_ct_df.merge(control_mean.rename('DeltaCt_control'), on='Gene')
delta_ct_df['DeltaDeltaCt'] = delta_ct_df['DeltaCt'] - delta_ct_df['DeltaCt_control']
delta_ct_df['FoldChange'] = 2 ** (-delta_ct_df['DeltaDeltaCt'])

# Построение графика с аннотацией
# Увеличенные шрифты и стиль
sns.set(style='whitegrid')
sns.set_context("poster", font_scale=1.8)
plt.figure(figsize=(14, 11))
ax = sns.barplot(
    data=delta_ct_df,
    x='Gene',
    y='FoldChange',
    hue='Sample',
    ci='sd',
    palette='Set2'
)
plt.yscale('log')
plt.ylabel('Fold Change (log scale)', fontsize=24)
plt.xlabel("")
plt.title(
    'Gene expression changes of EMT-associated genes after treatment\n'
    'of A549 cells with 6 µg/mL cisplatin during the first 24 hours',
    fontsize=25
)
plt.xticks(fontsize=20)
plt.yticks(fontsize=24)
# Формируем сравнения
samples = sorted(delta_ct_df['Sample'].unique())
comparisons = []
for gene in genes:
    for sample in samples:
        if sample != control_sample:
            comparisons.append(((gene, control_sample), (gene, sample)))

# Добавляем аннотацию с помощью statannotations
annotator = Annotator(
    pairs=comparisons,
    data=delta_ct_df,
    x='Gene',
    y='FoldChange',
    hue='Sample',
    ax=ax
)
annotator.configure(
    test='t-test_ind',
    text_format='star',
    loc='inside',
    verbose=0,
    fontsize=14,  # уменьшенный размер звёздочек
    pvalue_thresholds=[
        (0.1, "*"),
        (0.05, "**"),
        (0.01, "***"),
        (0.001, "****")
    ]
)

annotator.apply_and_annotate()

plt.legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=18, title_fontsize=20)
plt.tight_layout()
plt.savefig(output_path)
plt.show()
plt.close()

# Индивидуальные графики по генам
for gene in genes:
    df_gene = delta_ct_df[delta_ct_df['Gene'] == gene]
    samples_gene = df_gene['Sample'].unique()
    comparisons = [(control_sample, sample) for sample in samples_gene if sample != control_sample]

    plt.figure(figsize=(8, 5))
    ax = sns.barplot(
        data=df_gene,
        x='Sample',
        y='FoldChange',
        ci='sd',
        palette='coolwarm'
    )
    plt.yscale('log')
    plt.ylabel('Fold Change (log scale)', fontsize=18)
    plt.xlabel('', fontsize=16)
    plt.title(f'Fold Change for {gene} (ΔΔCt method)', fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    annotator = Annotator(
        ax=ax,
        pairs=comparisons,
        data=df_gene,
        x='Sample',
        y='FoldChange'
    )
    annotator.configure(
        test='t-test_ind',
        text_format='star',
        loc='inside',
        verbose=0,
        fontsize=10,  # уменьшенный размер звёздочек
        pvalue_thresholds=[
            (0.1, "*"),
            (0.05, "**"),
            (0.01, "***"),
            (0.001, "****")
        ]
    )
    annotator.apply_and_annotate()

    plt.tight_layout()
    plt.savefig(output_dir / f"{gene}_foldchange_annotated.png")
    plt.close()

# Сводная таблица
summary = (
    delta_ct_df.groupby(['Gene', 'Sample'])['FoldChange']
    .agg(['mean', 'std', 'count'])
    .rename(columns={'mean': 'FoldChange_Mean', 'std': 'FoldChange_SD', 'count': 'N'})
    .reset_index()
)

# Сохраняем в CSV
summary_output = f"{output_path} /fold_change_summary.csv"
summary.to_csv(summary_output, index=False, sep=';')
