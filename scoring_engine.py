import numpy as np

def literature_saturation_score(counts):
    """
    PubMed sayım dizisini [0, 1] aralığında bir 'literatür doygunluk skoruna' dönüştürür.
    Sigmoidin orta noktası log1p(sayı) dağılımının MEDYANINA, eğimi ise IQR'nin tersine bağlanır.
    Böylece eğri, taranan her gen kümesine kendiliğinden dinamik olarak uyum sağlar.
    """
    counts = np.asarray(counts, dtype=float)
    x = np.log1p(np.nan_to_num(counts, nan=0.0))
 
    median = np.median(x)
    q75, q25 = np.percentile(x, 75), np.percentile(x, 25)
    iqr = q75 - q25
    iqr = iqr if iqr > 1e-6 else 1.0
 
    k = 4.0 / iqr  # Eğimi dağılımın yayılımıyla ters orantılı kurgulayarak ayırt ediciliği korur
    score = 1.0 / (1.0 + np.exp(-k * (x - median)))
    return score

def chembl_phase_score(max_phases):
    """
    ChEMBL max_phase (0-4) değerlerini [0, 1] aralığında klinik doğrulama skoruna dönüştürür.
    Lojistik eğri Faz ~2.2 civarında merkezlenerek Faz 3 -> Faz 4 geçişindeki niteliksel
    sıçramayı (onaylı ilaç) güçlü şekilde ödüllendirir.
    """
    p = np.asarray(max_phases, dtype=float)
    p = np.nan_to_num(p, nan=0.0)
 
    raw = 1.0 / (1.0 + np.exp(-1.8 * (p - 2.2)))
    lo = 1.0 / (1.0 + np.exp(-1.8 * (0.0 - 2.2)))
    hi = 1.0 / (1.0 + np.exp(-1.8 * (4.0 - 2.2)))
    score = (raw - lo) / (hi - lo)
    return np.clip(score, 0.0, 1.0)

def entropy_weights(matrix, min_logfc_weight=0.5):
    """
    Shannon Entropisi kullanarak kriterlerin ayırt edicilik gücüne göre objektif ağırlık hesaplar.
    Sütun 0 = LogFC (Normalize), Sütun 1..k = Bilgi katmanı skorları ([0, 1]).
    LogFC'nin transkriptomik dominansını korumak için min_logfc_weight floor (taban) bariyeri uygulanır.
    """
    eps = 1e-12
    m = np.asarray(matrix, dtype=float)
 
    # Her sütunu [eps, 1] aralığına min-max normalize et
    col_min = m.min(axis=0)
    col_max = m.max(axis=0)
    denom = np.where((col_max - col_min) < eps, 1.0, col_max - col_min)
    norm = (m - col_min) / denom
    norm = np.clip(norm, eps, 1.0)
 
    # Bilgi Entropisi (Shannon Entropy) Hesaplaması
    P = norm / norm.sum(axis=0, keepdims=True)
    n = m.shape[0]
    k = 1.0 / np.log(n) if n > 1 else 1.0
    entropy = -k * (P * np.log(P)).sum(axis=0)
 
    diversity = np.clip(1.0 - entropy, eps, None)  # Ayırt edicilik derecesi
    raw_weights = diversity / diversity.sum()
 
    # Biyolojik sadakat kalkanı: LogFC taban ağırlığı kontrolü
    w_logfc = max(raw_weights[0], min_logfc_weight)
    remaining = 1.0 - w_logfc
 
    other_raw = raw_weights[1:]
    other_sum = other_raw.sum()
    if other_sum < eps:
        other_weights = np.full_like(other_raw, remaining / len(other_raw))
    else:
        other_weights = remaining * (other_raw / other_sum)
 
    return np.concatenate([[w_logfc], other_weights])
