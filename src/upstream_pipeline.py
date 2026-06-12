import os
import scanpy as sc
import numpy as np
import pandas as pd
from src.api_clients import fetch_pubmed_count, fetch_chembl_max_phase
from src.scoring_engine import literature_saturation_score, chembl_phase_score, entropy_weights

# --- GENEL AYARLAR VE KLASÖR YÖNETİMİ ---
sc.settings.verbosity = 3
RESULTS_DIR = './scrna_results'
os.makedirs(RESULTS_DIR, exist_ok=True)

def run_target_prioritization_pipeline(top_n_markers=15):
    print("--- ADIM 1 & 2: Veri Yükleme ve QC Metrikleri ---")
    adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()
    
    # İnsan genomu kurallarına göre mitokondri ve ribozom tespiti
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
    
    print("--- ADIM 3: Filtreleme Uygulanıyor ---")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :].copy()
    
    print("--- ADIM 4: Normalizasyon ve Log-Transform (DÜZELTİLDİ 🛡️) ---")
    adata.layers['counts'] = adata.X.copy()
    
    # KORUMA KALKANI: Önce toplam transkriptleri eşitle, matrisi GÜNCELLE, sonra log1p al!
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # Downstream analizler için log-normalize matrisi raw alanına kilitle
    
    # Biyoistatistiksel Doğrulama Kontrolü
    row_sums = adata.raw.X.sum(axis=1)
    if isinstance(row_sums, np.matrix): row_sums = row_sums.A1
    print(f"[BAŞARILI] Normalizasyon Sonrası Hücre Kütüphane Varyansı (0.0 olmalı): {np.var(row_sums):.6f}")
    
    print("--- ADIM 5 - 9: HVG, Ölçekleme ve Boyut Azaltma ---")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    
    print("--- ADIM 10 - 12: Kümeleme (Leiden) ve Marker Analizi ---")
    # Faz 2 ölçeklenebilirliği için igraph altyapısı zorunlu kılınmıştır
    try:
        sc.tl.leiden(adata, resolution=0.5, flavor='igraph', n_iterations=2)
    except Exception:
        sc.tl.leiden(adata, resolution=0.5) # Fallback mekanizması
        
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    # --- ADIM 13: AGENTIC AI KNOWLEDGE-HARMONIZATION LAYER ---
    print("--- ADIM 13: Agentic Bilgi Harmanlama ve Hedef Önceliklendirme ---")
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    candidate_rows = []
    for group in groups:
        names = result['names'][group][:top_n_markers]
        lfc = result['logfoldchanges'][group][:top_n_markers]
        pvals_adj = result['pvals_adj'][group][:top_n_markers]
        for gene, fc, p in zip(names, lfc, pvals_adj):
            candidate_rows.append({
                "cluster": group, "gene": gene, "logFC": float(fc), "pval_adj": float(p)
            })
     
    master_df = pd.DataFrame(candidate_rows)
    unique_genes = master_df['gene'].unique().tolist()
    
    # Önbellekli (Cache) Canlı API Harvesting
    pubmed_cache = {}
    chembl_cache = {}
    for i, gene in enumerate(unique_genes):
        print(f"  [{i + 1}/{len(unique_genes)}] API Sorgulanıyor: {gene}")
        pubmed_cache[gene] = fetch_pubmed_count(gene)
        chembl_cache[gene] = fetch_chembl_max_phase(gene)
     
    master_df['pubmed_count'] = master_df['gene'].map(pubmed_cache)
    master_df['chembl_max_phase'] = master_df['gene'].map(chembl_cache)
    
    # Eksik verileri sıfırlama sigortası
    master_df['pubmed_count'] = master_df['pubmed_count'].fillna(0)
    master_df['chembl_max_phase'] = master_df['chembl_max_phase'].fillna(0)
     
    # Non-Lineer Skorlama Dönüşümleri
    master_df['literature_score'] = literature_saturation_score(master_df['pubmed_count'].values)
    master_df['chembl_score'] = chembl_phase_score(master_df['chembl_max_phase'].values)
     
    # LogFC Global Min-Max Normalizasyonu
    lfc = master_df['logFC'].values.astype(float)
    master_df['logfc_norm'] = (lfc - lfc.min()) / (lfc.max() - lfc.min() + 1e-12)
     
    # Entropi Tabanlı Çok Kriterli Karar Verme (MCDA) Ağırlıkları
    criteria_matrix = master_df[['logfc_norm', 'literature_score', 'chembl_score']].values
    weights = entropy_weights(criteria_matrix, min_logfc_weight=0.5)
     
    print("\n📊 Dinamik Hesaplanan MCDA Ağırlıkları:")
    print(f"  - LogFC (Biyolojik Dominans):     {weights[0]:.3f}")
    print(f"  - PubMed Literatür Doygunluğu:    {weights[1]:.3f}")
    print(f"  - ChEMBL Klinik Faz Güvencesi:    {weights[2]:.3f}")
     
    # Nihai Skorlama ve Sıralama
    master_df['priority_score'] = criteria_matrix @ weights
    master_df = master_df.sort_values(['cluster', 'priority_score'], ascending=[True, False]).reset_index(drop=True)
     
    output_csv = f'{RESULTS_DIR}/target_prioritization_master.csv'
    master_df.to_csv(output_csv, index=False)
    print(f"\n[BAŞARILI] Önceliklendirme tablosu diske mühürlendi: {output_csv}")
    return master_df

if __name__ == "__main__":
    run_target_prioritization_pipeline()
