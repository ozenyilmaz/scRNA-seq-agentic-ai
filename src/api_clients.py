import time
import requests
import numpy as np
from requests.adapters import HTTPAdapter, Retry

def get_requests_session():
    """Geri-deneme (retry) mantığı içeren paylaşılan bir HTTP oturumu döndürür."""
    session = requests.Session()
    retries = Retry(total=3, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504])
    session.mount('https://', HTTPAdapter(max_retries=retries))
    return session

SESSION = get_requests_session()
PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
CHEMBL_TARGET = "https://www.ebi.ac.uk/chembl/api/data/target.json"
CHEMBL_MECHANISM = "https://www.ebi.ac.uk/chembl/api/data/mechanism.json"
CHEMBL_MOLECULE = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
DISEASE_CONTEXT = "heart failure"

def fetch_pubmed_count(gene_symbol, disease_context=DISEASE_CONTEXT, session=SESSION, sleep=0.34):
    """Gen sembolü ile hastalık bağlamının PubMed ortak makale sayısını döndürür."""
    term = f"{gene_symbol}[Title/Abstract] AND {disease_context}[Title/Abstract]"
    params = {"db": "pubmed", "term": term, "retmode": "json", "retmax": 0}
    try:
        r = session.get(PUBMED_ESEARCH, params=params, timeout=15)
        r.raise_for_status()
        count = int(r.json().get("esearchresult", {}).get("count", 0))
    except Exception:
        count = np.nan
    time.sleep(sleep)
    return count

def fetch_chembl_max_phase(gene_symbol, session=SESSION, sleep=0.2):
    """ChEMBL ilişkisel şemasını tip-güvenlikli olarak sorgular ve max_phase döndürür."""
    try:
        # Adım A: İnsan hedef protein kimliği yakalama
        params = {"target_synonym__iexact": gene_symbol, "organism__iexact": "Homo sapiens", "format": "json"}
        r = session.get(CHEMBL_TARGET, params=params, timeout=15)
        r.raise_for_status()
        targets = r.json().get("targets", [])
        time.sleep(sleep)
        if not targets:
            return 0
            
        target_ids = [t["target_chembl_id"] for t in targets if "target_chembl_id" in t]
        
        # Adım B: Hedefe bağlı molekül kimliklerini toplama
        molecule_ids = set()
        for tid in target_ids:
            r = session.get(CHEMBL_MECHANISM, params={"target_chembl_id": tid, "format": "json"}, timeout=15)
            r.raise_for_status()
            for m in r.json().get("mechanisms", []):
                mol_id = m.get("molecule_chembl_id")
                if mol_id:
                    molecule_ids.add(mol_id)
            time.sleep(sleep)
            
        if not molecule_ids:
            return 0
            
        # Adım C: Gerçek max_phase değerini tip-güvenlikli okuma
        max_phase = 0
        molecule_ids = list(molecule_ids)
        CHUNK = 25
        for i in range(0, len(molecule_ids), CHUNK):
            chunk = molecule_ids[i:i + CHUNK]
            r = session.get(CHEMBL_MOLECULE, params={"molecule_chembl_id__in": ",".join(chunk), "format": "json"}, timeout=15)
            r.raise_for_status()
            for mol in r.json().get("molecules", []):
                raw_phase = mol.get("max_phase")
                if raw_phase is not None:
                    try:
                        # Tip Güvenlik Kalkanı: '4.0' string'ini float üzerinden int'e güvenle döküyoruz
                        phase = int(float(raw_phase))
                        if phase > max_phase:
                            max_phase = phase
                    except (ValueError, TypeError):
                        continue
            time.sleep(sleep)
            
        return max_phase
    except Exception:
        return np.nan
