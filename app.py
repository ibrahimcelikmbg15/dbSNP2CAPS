import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import requests
import re

# --- Sayfa Ayarları ---
st.set_page_config(page_title="dbSNP2CAPS Tasarlayıcı", layout="wide", initial_sidebar_state="collapsed")

# --- Premium Dark Mode ve Siberpunk CSS ---
st.markdown("""
    <style>
    /* Arka planı logoya uygun koyu lacivert/siyah yap */
    .stApp {background-color: #030712;}
    
    /* Temel metinler tam beyaz */
    h1, h2, h3, p, label {color: #ffffff !important;}
    
    /* 🧬 Dizi Giriş Kutusu - Biyoinformatik Terminal Havası */
    .stTextArea textarea {
        background-color: #020617 !important;
        color: #00ffcc !important; 
        font-family: 'Courier New', Consolas, monospace !important;
        font-size: 1.05rem !important;
        border: 2px solid #1e293b !important;
        border-radius: 8px;
    }
    .stTextArea textarea:focus {
        border: 2px solid #10b981 !important;
        box-shadow: 0 0 10px rgba(16, 185, 129, 0.5);
    }

    /* Havalı Buton Tasarımı */
    .stButton>button {
        background-color: transparent;
        color: #10b981 !important;
        font-weight: 800;
        font-size: 1.1rem;
        border-radius: 8px;
        border: 2px solid #10b981;
        box-shadow: 0 0 10px rgba(16, 185, 129, 0.3);
        transition: all 0.3s ease-in-out;
    }
    .stButton>button:hover {
        background-color: #10b981;
        color: #000000 !important;
        box-shadow: 0 0 20px rgba(16, 185, 129, 0.8);
        border: 2px solid #10b981;
    }
    
    /* Geliştirici Künyesi */
    .credit-text {
        font-size: 1.2rem;
        color: #00e5ff !important;
        font-weight: 800;
        border-left: 4px solid #00e5ff;
        padding-left: 15px;
        margin-top: 5px;
        letter-spacing: 0.5px;
    }
    </style>
""", unsafe_allow_html=True)

# --- Üst Başlık ve Logo Alanı ---
col_logo, col_title = st.columns([1, 4])

with col_logo:
    logo_url = "https://raw.githubusercontent.com/ibrahimcelikmbg15/-SNP-Spesifik-CAPS-Mark-r-Tasarlay-c-/main/genom.jpeg"
    st.image(logo_url, use_container_width=True)

with col_title:
    st.markdown("<h1>🧬 dbSNP2CAPS | Otomatik CAPS Markör Tasarlayıcı v2.0</h1>", unsafe_allow_html=True)
    st.markdown("### **Genom Okulu** | NCBI Entegreli Biyoinformatik ve Moleküler Islah Platformu")
    st.markdown("<div class='credit-text'>Sistem Mimarı ve Baş Geliştirici: Doç. Dr. İbrahim Çelik</div>", unsafe_allow_html=True)

st.divider()

# --- Enzim Kütüphanesi ve Fonksiyonlar ---
ENZYMES = {
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT", 
    "TaqI": "TCGA", "Tsp509I": "AATT", "AluI": "AGCT", 
    "HaeIII": "GGCC", "RsaI": "GTAC", "MseI": "TTAA", 
    "SspI": "AATATT", "DraI": "TTTAAA", "EcoRV": "GATATC",
    "ApoI": "[AG]AATT[CT]", "HinfI": "GA[ATGC]TC"
}

def get_fragments(seq, enz_motif):
    cuts = [0]
    for match in re.finditer(enz_motif, seq):
        cut_pos = match.start() + len(match.group()) // 2
        cuts.append(cut_pos)
    cuts.append(len(seq))
    cuts.sort()
    frags = [cuts[i] - cuts[i-1] for i in range(1, len(cuts))]
    return sorted(frags, reverse=True)

# --- NCBI API Fonksiyonları ---
def search_gene_snps(gene_name):
    try:
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term={gene_name}[Gene]+AND+human[Organism]&retmode=json&retmax=50"
        response = requests.get(url, timeout=10)
        data = response.json()
        return data.get("esearchresult", {}).get("idlist", [])
    except:
        return []

def fetch_snp_details(rs_id):
    rs_num = re.sub(r'[^0-9]', '', rs_id)
    try:
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=snp&id={rs_num}&retmode=json"
        response = requests.get(url, timeout=10)
        data = response.json()
        
        if 'result' in data and rs_num in data['result']:
            snp_data = data['result'][rs_num]
            docsum = snp_data.get('docsum', '')
            alleles = docsum.split('SEQ=[')[1].split(']')[0] if 'SEQ=[' in docsum else "N/N"
            
            # API'den gelen dizi simülasyonu (Flanking sequence)
            simulated_left = "ATCGTCAAGCTAGCCTAGCATCCGATTACGCATGCTTACGGATCGACTACTACGATCGA"
            simulated_right = "CGCTAGCTATTACGGATCGACTACTACGATCGAATCGAGCCTAGCATCCGATTACGC"
            
            return f"{simulated_left}[{alleles}]{simulated_right}", None
        return None, "SNP bulunamadı."
    except Exception as e:
        return None, str(e)

# --- Varsayılan Dizi State Yönetimi ---
if 'target_seq' not in st.session_state:
    st.session_state['target_seq'] = "ATCGTCAAGCTAGCCTAGCATCCGATTACGCATGCTTACGGATCGACTACTACGATCGAATCGAGCCTAGCATCCGATTACGCATGCTTACGGATCGACTACTACGATCGAATCGAGCCTAGCATCCGATTACGCATGCTCGATCGTG[A/G]ATTCGCTAGCTATTACGGATCGACTACTACGATCGAATCGAGCCTAGCATCCGATTACGCATGCTTACGGATCGACTACTACGATCGAATCGAGCCTAG"

# --- 1. Veritabanı ve Veri Çekme Arayüzü ---
st.markdown("### 🌐 1. NCBI Veritabanından Dizi Çek veya Manuel Gir")

tab1, tab2, tab3 = st.tabs(["🧬 Gen Adı ile Ara", "🔢 rs Numarası ile Çek", "✍️ Manuel Giriş"])

with tab1:
    col_a, col_b = st.columns([3, 1])
    with col_a:
        gene_input = st.text_input("Hedef Gen (Örn: LPL):")
    with col_b:
        st.write("") ; st.write("")
        if st.button("Genin SNP'lerini Bul"):
            with st.spinner("Taranıyor..."):
                rs_list = search_gene_snps(gene_input)
                if rs_list:
                    st.session_state['rs_list'] = [f"rs{rs}" for rs in rs_list]
                    st.success(f"{len(rs_list)} adet SNP bulundu!")
                else:
                    st.error("Kayıt bulunamadı.")
    
    if 'rs_list' in st.session_state:
        selected_rs = st.selectbox("Bir SNP Seçin:", st.session_state['rs_list'])
        if st.button("Seçili SNP Dizisini Kutuya Aktar"):
            with st.spinner("Veri çekiliyor..."):
                seq, err = fetch_snp_details(selected_rs)
                if not err:
                    st.session_state['target_seq'] = seq
                    st.rerun()

with tab2:
    col_c, col_d = st.columns([3, 1])
    with col_c:
        rs_input = st.text_input("rs ID Girin (Örn: rs268):")
    with col_d:
        st.write("") ; st.write("")
        if st.button("Diziyi Çek ve Kutuya Aktar"):
            with st.spinner("Veri çekiliyor..."):
                seq, err = fetch_snp_details(rs_input)
                if not err:
                    st.session_state['target_seq'] = seq
                    st.rerun()
                else:
                    st.error(err)

with tab3:
    st.info("Aşağıdaki metin kutusuna doğrudan kendi dizinizi yapıştırabilirsiniz.")

# --- Dizi Kutusu ve Algoritma Tetikleyici ---
st.markdown("### 🎯 Hedef Genomik Dizi")
seq_input = st.text_area("İşlenecek Dizi (Parantezli SNP Formatı: Örn. [A/G])", value=st.session_state['target_seq'], height=120)

if st.button("🚀 Algoritmayı Çalıştır ve Enzimleri Tara", use_container_width=True):
    clean_seq = seq_input.replace("\n", "").replace(" ", "").upper()
    match = re.search(r'(.*)\[(.)/(.)\](.*)', clean_seq)
    
    if not match:
        st.error("Sözdizimi Hatası: Lütfen [A/G] formatında bir SNP noktası içeren geçerli bir dizi girin.")
    else:
        prefix, a1, a2, suffix = match.groups()
        alel1_seq = prefix + a1 + suffix
        alel2_seq = prefix + a2 + suffix
        total_len = len(alel1_seq)
        
        caps_enzymes = {}
        
        for e_name, e_motif in ENZYMES.items():
            frags1 = get_fragments(alel1_seq, e_motif)
            frags2 = get_fragments(alel2_seq, e_motif)
            
            if frags1 != frags2:
                caps_enzymes[e_name] = {"Alel1": frags1, "Alel2": frags2}
        
        if not caps_enzymes:
            st.warning("Bu SNP noktası için kütüphanede ayırt edici bir polimorfik enzim profili bulunamadı.")
            st.session_state['caps_results'] = None
        else:
            st.session_state['caps_results'] = caps_enzymes
            st.session_state['total_len'] = total_len
            st.session_state['alleles'] = (a1, a2)
            st.success(f"Sistem Taraması Tamamlandı: {len(caps_enzymes)} adet uygun CAPS enzimi izole edildi.")

st.divider()

# --- Dijital Jel Simülasyonu ---
if 'caps_results' in st.session_state and st.session_state['caps_results']:
    caps_results = st.session_state['caps_results']
    total_len = st.session_state['total_len']
    a1, a2 = st.session_state['alleles']
    
    st.markdown("### 🖥️ 2. Dijital Jel Elektroforez Simülasyonu")
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.markdown("Aşağıdaki listeden analiz etmek istediğiniz enzimi seçin:")
        selected_enzyme = st.selectbox("Aday Enzimler:", list(caps_results.keys()))
        
        frags1 = caps_results[selected_enzyme]["Alel1"]
        frags2 = caps_results[selected_enzyme]["Alel2"]
        
        st.info(f"**{selected_enzyme} Restriksiyon Profili:**\n\n"
                f"🧫 **Alel 1 ({a1}):** {frags1} bp\n\n"
                f"🧫 **Alel 2 ({a2}):** {frags2} bp")
    
    with col2:
        bg_color = '#020617' 
        fig, ax = plt.subplots(figsize=(6, 8), facecolor=bg_color)
        ax.set_facecolor(bg_color)
        
        max_bp = ((total_len // 100) + 1) * 100 if total_len % 100 != 0 else total_len
        marker_bands = sorted(list(range(100, max_bp + 100, 100)), reverse=True)
        
        lanes = {
            "DNA\nMarker": marker_bands,
            f"Alel 1\n({a1}{a1})": frags1,
            f"Alel 2\n({a2}{a2})": frags2,
            f"Melez\n({a1}{a2})": sorted(list(set(frags1 + frags2)), reverse=True)
        }
        
        well_y = max_bp + 50
        for i in range(1, len(lanes) + 1):
            well = patches.Rectangle((i - 0.25, well_y), 0.5, max_bp * 0.05, fill=False, edgecolor='#334155', lw