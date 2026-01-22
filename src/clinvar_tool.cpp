// clinvar_tool.cpp
// Analiza wariantów z VCF z użyciem ClinVar + gnomAD + COSMIC + TP53db + OncoKB/CIViC
// Tryby: "sclc" (TP53/RB1, somatic P/LP, lung/SCLC) oraz "global"
//
// Kompilacja:
//  - C++17
//  - zlib (dla .gz)

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <zlib.h> // gzopen, gzgets, gzclose

// ==================== NARZĘDZIA POMOCNICZE ====================

static std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, '\t')) out.push_back(item);
    return out;
}

static std::string trim(const std::string& s) {
    size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

static std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c) { return (unsigned char)std::tolower(c); });
    return s;
}

static std::string make_coord_key(const std::string& chr,
    const std::string& pos,
    const std::string& ref,
    const std::string& alt) {
    return chr + ":" + pos + ":" + ref + ":" + alt;
}

static std::string join_semicolon(const std::vector<std::string>& vec) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i) oss << ";";
        oss << vec[i];
    }
    return oss.str();
}

// ==================== LineReader: obsługa .txt i .gz ====================

class LineReader {
public:
    explicit LineReader(const std::string& path)
        : is_gz_(ends_with(path, ".gz")), gz_(nullptr) {
        if (is_gz_) {
            gz_ = gzopen(path.c_str(), "rb");
            if (!gz_) throw std::runtime_error("Cannot open gz file: " + path);
        }
        else {
            in_.open(path);
            if (!in_) throw std::runtime_error("Cannot open file: " + path);
        }
    }

    ~LineReader() {
        if (is_gz_ && gz_) gzclose(gz_);
    }

    bool getline(std::string& out) {
        out.clear();
        if (is_gz_) {
            const int BUF = 1 << 16;
            static thread_local std::vector<char> buf(BUF);
            if (!gzgets(gz_, buf.data(), BUF)) return false;
            out.assign(buf.data());
            while (!out.empty() && (out.back() == '\n' || out.back() == '\r')) out.pop_back();
            return true;
        }
        return static_cast<bool>(std::getline(in_, out));
    }

    bool good() const {
        return is_gz_ ? (gz_ != nullptr) : static_cast<bool>(in_);
    }

private:
    bool is_gz_;
    std::ifstream in_;
    gzFile gz_;

    static bool ends_with(const std::string& s, const std::string& suf) {
        if (s.size() < suf.size()) return false;
        return s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
    }
};

// ==================== STRUKTURY DANYCH ====================

struct ClinVarRecord {
    std::string variation_id;
    std::string allele_id; // mapping z variation_allele.txt
    std::string gene_symbol;
    std::string clinical_significance;
    std::string origin;
    std::string phenotype_list;
    std::string chrom;
    std::string pos;
    std::string ref;
    std::string alt;
};

struct HgvsInfo {
    std::vector<std::string> nuc;  // HGVS nucleotide
    std::vector<std::string> prot; // HGVS protein
};

struct GnomADInfo { std::string af; };

struct CosmicInfo {
    std::string ids;
    std::string primary_site;
    std::string histology;
    std::string count;
};

struct TP53FuncInfo {
    std::string classification;
    std::string func_score;
    std::string notes;
};

struct TherapyInfo {
    std::string evidence_level;
    std::string drug;
    std::string disease;
    std::string source;
};

// ==================== Filtr fenotypu SCLC / lung ====================

static bool is_sclc_or_lung(const std::string& phenotype_list) {
    std::string low = to_lower(phenotype_list);
    if (low.find("small cell lung") != std::string::npos) return true;
    if (low.find("sclc") != std::string::npos) return true;
    if (low.find("lung") != std::string::npos) return true;
    if (low.find("pulmonary") != std::string::npos) return true;
    return false;
}

// ==================== variation_allele: VariationID -> AlleleID ====================

static std::unordered_map<std::string, std::string>
load_variation_allele(const std::string& path) {
    std::unordered_map<std::string, std::string> out;
    if (path.empty() || path == "-") return out;

    LineReader lr(path);
    if (!lr.good()) return out;

    std::string header;
    if (!lr.getline(header)) return out;

    auto h = split_tab(header);
    std::unordered_map<std::string, size_t> idx;
    for (size_t i = 0; i < h.size(); ++i) idx[h[i]] = i;

    int iVar = idx.count("VariationID") ? (int)idx["VariationID"] : -1;
    int iAll = idx.count("AlleleID") ? (int)idx["AlleleID"] : -1;
    if (iVar < 0 || iAll < 0) return out;

    std::string line;
    while (lr.getline(line)) {
        if (line.empty()) continue;
        auto c = split_tab(line);
        if ((int)c.size() <= std::max(iVar, iAll)) continue;
        out[c[iVar]] = c[iAll];
    }
    return out;
}

// ==================== ClinVar: variant_summary (filtrowanie) ====================

static std::unordered_map<std::string, ClinVarRecord>
load_clinvar_filtered(const std::string& variant_summary_path,
    const std::unordered_map<std::string, std::string>& varid_to_alleleid,
    bool sclc_mode)
{
    std::unordered_map<std::string, ClinVarRecord> result;

    LineReader lr(variant_summary_path);
    if (!lr.good()) {
        std::cerr << "ERROR: cannot open variant_summary\n";
        return result;
    }

    std::string header_line;
    if (!lr.getline(header_line)) {
        std::cerr << "ERROR: empty variant_summary\n";
        return result;
    }

    auto header_cols = split_tab(header_line);
    std::unordered_map<std::string, size_t> col_idx;
    for (size_t i = 0; i < header_cols.size(); ++i) col_idx[header_cols[i]] = i;

    auto get_idx = [&](const std::string& name) -> int {
        auto it = col_idx.find(name);
        if (it == col_idx.end()) return -1;
        return (int)it->second;
        };

    int idx_VariationID = get_idx("VariationID");
    int idx_GeneSymbol = get_idx("GeneSymbol");
    int idx_ClinicalSignificance = get_idx("ClinicalSignificance");
    int idx_Origin = get_idx("Origin");
    int idx_Chromosome = get_idx("Chromosome");
    int idx_Start = get_idx("Start");
    int idx_RefAllele = get_idx("ReferenceAllele");
    int idx_AltAllele = get_idx("AlternateAllele");
    int idx_Phenotype = get_idx("PhenotypeList");
    if (idx_Phenotype < 0) idx_Phenotype = get_idx("Disease/Phenotype");

    if (idx_VariationID < 0 || idx_GeneSymbol < 0 ||
        idx_ClinicalSignificance < 0 || idx_Origin < 0 ||
        idx_Chromosome < 0 || idx_Start < 0 ||
        idx_RefAllele < 0 || idx_AltAllele < 0) {
        std::cerr << "ERROR: required columns not found in variant_summary header\n";
        return result;
    }

    std::string line;
    while (lr.getline(line)) {
        if (line.empty()) continue;
        auto cols = split_tab(line);
        if ((int)cols.size() <= idx_AltAllele) continue;

        std::string variation_id = cols[idx_VariationID];
        std::string gene_symbol = trim(cols[idx_GeneSymbol]);
        std::string clinsig = trim(cols[idx_ClinicalSignificance]);
        std::string origin = trim(cols[idx_Origin]);
        std::string chrom = cols[idx_Chromosome];
        std::string pos = cols[idx_Start];
        std::string ref = cols[idx_RefAllele];
        std::string alt = cols[idx_AltAllele];
        std::string phen_list = (idx_Phenotype >= 0 && (int)cols.size() > idx_Phenotype)
            ? cols[idx_Phenotype] : "";

        if (sclc_mode) {
            if (gene_symbol != "TP53" && gene_symbol != "RB1") continue;
            if (to_lower(origin).find("somatic") == std::string::npos) continue;
            if (to_lower(clinsig).find("pathogenic") == std::string::npos) continue;
            if (idx_Phenotype >= 0 && !is_sclc_or_lung(phen_list)) continue;
        }

        ClinVarRecord rec;
        rec.variation_id = variation_id;
        rec.allele_id = "";
        auto itA = varid_to_alleleid.find(variation_id);
        if (itA != varid_to_alleleid.end()) rec.allele_id = itA->second;
        rec.gene_symbol = gene_symbol;
        rec.clinical_significance = clinsig;
        rec.origin = origin;
        rec.phenotype_list = phen_list;
        rec.chrom = chrom;
        rec.pos = pos;
        rec.ref = ref;
        rec.alt = alt;

        result[make_coord_key(chrom, pos, ref, alt)] = rec;
    }

    return result;
}

// ==================== HGVS: hgvs4variation ====================

static std::unordered_map<std::string, HgvsInfo>
load_hgvs_for_variants(const std::string& hgvs4variation_path,
    const std::unordered_set<std::string>& wanted_variation_ids)
{
    std::unordered_map<std::string, HgvsInfo> out;

    LineReader lr(hgvs4variation_path);
    std::string header_line;
    if (!lr.getline(header_line)) {
        std::cerr << "ERROR: empty hgvs4variation\n";
        return out;
    }

    auto header_cols = split_tab(header_line);
    std::unordered_map<std::string, size_t> col_idx;
    for (size_t i = 0; i < header_cols.size(); ++i) col_idx[header_cols[i]] = i;

    auto get_idx = [&](const std::string& name) -> int {
        auto it = col_idx.find(name);
        if (it == col_idx.end()) return -1;
        return (int)it->second;
        };

    int idx_VariationID = get_idx("VariationID");
    int idx_HGVS = get_idx("HGVS");
    int idx_Type = get_idx("Type");

    if (idx_VariationID < 0 || idx_HGVS < 0 || idx_Type < 0) {
        std::cerr << "ERROR: required columns not found in hgvs4variation\n";
        return out;
    }

    std::string line;
    while (lr.getline(line)) {
        if (line.empty()) continue;
        auto cols = split_tab(line);
        if ((int)cols.size() <= idx_Type) continue;

        std::string variation_id = cols[idx_VariationID];
        if (!wanted_variation_ids.empty() &&
            wanted_variation_ids.find(variation_id) == wanted_variation_ids.end()) {
            continue;
        }

        std::string hgvs = cols[idx_HGVS];
        std::string type = cols[idx_Type];

        auto& info = out[variation_id];
        if (type == "nucleotide") info.nuc.push_back(hgvs);
        else if (type == "protein") info.prot.push_back(hgvs);
    }

    return out;
}

// ==================== gnomAD TSV ====================

static std::unordered_map<std::string, GnomADInfo>
load_gnomad(const std::string& path) {
    std::unordered_map<std::string, GnomADInfo> out;
    if (path == "-" || path.empty()) return out;

    LineReader lr(path);
    if (!lr.good()) {
        std::cerr << "WARN: cannot open gnomAD file: " << path << "\n";
        return out;
    }

    std::string header;
    if (!lr.getline(header)) {
        std::cerr << "WARN: empty gnomAD file\n";
        return out;
    }

    auto header_cols = split_tab(header);
    std::unordered_map<std::string, size_t> col_idx;
    for (size_t i = 0; i < header_cols.size(); ++i) col_idx[header_cols[i]] = i;

    auto get_idx = [&](const std::string& name) -> int {
        auto it = col_idx.find(name);
        if (it == col_idx.end()) return -1;
        return (int)it->second;
        };

    int idx_CHROM = get_idx("CHROM");
    int idx_POS = get_idx("POS");
    int idx_REF = get_idx("REF");
    int idx_ALT = get_idx("ALT");
    int idx_AF = get_idx("AF");

    if (idx_CHROM < 0 || idx_POS < 0 || idx_REF < 0 || idx_ALT < 0 || idx_AF < 0) {
        std::cerr << "WARN: missing columns in gnomAD TSV\n";
        return out;
    }

    std::string line;
    while (lr.getline(line)) {
        if (line.empty()) continue;
        auto cols = split_tab(line);
        if ((int)cols.size() <= idx_AF) continue;

        std::string chr = cols[idx_CHROM];
        std::string pos = cols[idx_POS];
        std::string ref = cols[idx_REF];
        std::string alt = cols[idx_ALT];
        std::string af = cols[idx_AF];

        out[make_coord_key(chr, pos, ref, alt)] = GnomADInfo{ af };
    }

    return out;
}

// ==================== COSMIC TSV ====================

static std::unordered_map<std::string, CosmicInfo>
load_cosmic(const std::string& path) {
    std::unordered_map<std::string, CosmicInfo> out;
    if (path == "-" || path.empty()) return out;

    LineReader lr(path);
    if (!lr.good()) {
        std::cerr << "WARN: cannot open COSMIC file: " << path << "\n";
        return out;
    }

    std::string header;
    if (!lr.getline(header)) {
        std::cerr << "WARN: empty COSMIC file\n";
        return out;
    }

    auto header_cols = split_tab(header);
    std::unordered_map<std::string, size_t> col_idx;
    for (size_t i = 0; i < header_cols.size(); ++i) col_idx[header_cols[i]] = i;

    auto get_idx = [&](const std::string& name) -> int {
        auto it = col_idx.find(name);
        if (it == col_idx.end()) return -1;
        return (int)it->second;
        };

    int idx_CHROM = get_idx("CHROM");
    int idx_POS = get_idx("POS");
    int idx_REF = get_idx("REF");
    int idx_ALT = get_idx("ALT");
    int idx_IDS = get_idx("COSMIC_IDS");
    int idx_SITE = get_idx("PRIMARY_SITE");
    int idx_HIST = get_idx("HISTOLOGY");
    int idx_COUNT = get_idx("COUNT");

    if (idx_CHROM < 0 || idx_POS < 0 || idx_REF < 0 || idx_ALT < 0 ||
        idx_IDS < 0 || idx_SITE < 0 || idx_HIST < 0 || idx_COUNT < 0) {
        std::cerr << "WARN: missing columns in COSMIC TSV\n";
        return out;
    }

    std::string line;
    while (lr.getline(line)) {
        if (line.empty()) continue;
        auto cols = split_tab(line);
        if ((int)cols.size() <= idx_COUNT) continue;

        std::string chr = cols[idx_CHROM];
        std::string pos = cols[idx_POS];
        std::string ref = cols[idx_REF];
        std::string alt = cols[idx_ALT];

        CosmicInfo ci;
        ci.ids = cols[idx_IDS];
        ci.primary_site = cols[idx_SITE];
        ci.histology = cols[idx_HIST];
        ci.count = cols[idx_COUNT];

        out[make_coord_key(chr, pos, ref, alt)] = ci;
    }

    return out;
}

// ==================== TP53db TSV ====================

static std::unordered_map<std::string, TP53FuncInfo>
load_tp53db(const std::string& path) {
    std::unordered_map<std::string, TP53FuncInfo> out;
    if (path == "-" || path.empty()) return out;

    LineReader lr(path);
    if (!lr.good()) {
        std::cerr << "WARN: cannot open TP53db file: " << path << "\n";
        return out;
    }

    std::string header;
    if (!lr.getline(header)) {
        std::cerr << "WARN: empty TP53db file\n";
        return out;
    }

    auto header_cols = split_tab(header);
    std::unordered_map<std::string, size_t> col_idx;
    for (size_t i = 0; i < header_cols.size(); ++i) col_idx[header_cols[i]] = i;

    auto get_idx = [&](const std::string& name) -> int {
        auto it = col_idx.find(name);
        if (it == col_idx.end()) return -1;
        return (int)it->second;
        };

    int idx_HGVS_P = get_idx("HGVS_P");
    int idx_CLASS = get_idx("CLASS");
    int idx_SCORE = get_idx("FUNC_SCORE");
    int idx_NOTES = get_idx("NOTES");

    if (idx_HGVS_P < 0) {
        std::cerr << "WARN: TP53db must contain column HGVS_P\n";
        return out;
    }

    std::string line;
    while (lr.getline(line)) {
        if (line.empty()) continue;
        auto cols = split_tab(line);
        if ((int)cols.size() <= idx_HGVS_P) continue;

        std::string hgvsp = cols[idx_HGVS_P];

        TP53FuncInfo info;
        info.classification = (idx_CLASS >= 0 && (int)cols.size() > idx_CLASS) ? cols[idx_CLASS] : "";
        info.func_score = (idx_SCORE >= 0 && (int)cols.size() > idx_SCORE) ? cols[idx_SCORE] : "";
        info.notes = (idx_NOTES >= 0 && (int)cols.size() > idx_NOTES) ? cols[idx_NOTES] : "";

        out[hgvsp] = info;
    }

    return out;
}

// ==================== Therapy TSV (OncoKB/CIViC) ====================

static std::unordered_map<std::string, TherapyInfo>
load_therapy(const std::string& path) {
    std::unordered_map<std::string, TherapyInfo> out;
    if (path == "-" || path.empty()) return out;

    LineReader lr(path);
    if (!lr.good()) {
        std::cerr << "WARN: cannot open therapy file: " << path << "\n";
        return out;
    }

    std::string header;
    if (!lr.getline(header)) {
        std::cerr << "WARN: empty therapy file\n";
        return out;
    }

    auto header_cols = split_tab(header);
    std::unordered_map<std::string, size_t> col_idx;
    for (size_t i = 0; i < header_cols.size(); ++i) col_idx[header_cols[i]] = i;

    auto get_idx = [&](const std::string& name) -> int {
        auto it = col_idx.find(name);
        if (it == col_idx.end()) return -1;
        return (int)it->second;
        };

    int idx_HGVS_P = get_idx("HGVS_P");
    int idx_EVID = get_idx("EVIDENCE_LEVEL");
    int idx_DRUG = get_idx("DRUG");
    int idx_DIS = get_idx("DISEASE");
    int idx_SRC = get_idx("SOURCE");

    if (idx_HGVS_P < 0) {
        std::cerr << "WARN: therapy TSV must contain HGVS_P\n";
        return out;
    }

    std::string line;
    while (lr.getline(line)) {
        if (line.empty()) continue;
        auto cols = split_tab(line);
        if ((int)cols.size() <= idx_HGVS_P) continue;

        std::string hgvsp = cols[idx_HGVS_P];

        TherapyInfo info;
        info.evidence_level = (idx_EVID >= 0 && (int)cols.size() > idx_EVID) ? cols[idx_EVID] : "";
        info.drug = (idx_DRUG >= 0 && (int)cols.size() > idx_DRUG) ? cols[idx_DRUG] : "";
        info.disease = (idx_DIS >= 0 && (int)cols.size() > idx_DIS) ? cols[idx_DIS] : "";
        info.source = (idx_SRC >= 0 && (int)cols.size() > idx_SRC) ? cols[idx_SRC] : "";

        out[hgvsp] = info;
    }

    return out;
}

// ==================== ANALIZA VCF (SCLC) ====================

static void analyze_vcf_sclc(
    const std::string& vcf_path,
    const std::unordered_map<std::string, ClinVarRecord>& clinvar_by_coord,
    const std::unordered_map<std::string, HgvsInfo>& hgvs_by_varid,
    const std::unordered_map<std::string, GnomADInfo>& gnomad_by_coord,
    const std::unordered_map<std::string, CosmicInfo>& cosmic_by_coord,
    const std::unordered_map<std::string, TP53FuncInfo>& tp53_by_hgvsp,
    const std::unordered_map<std::string, TherapyInfo>& therapy_by_hgvsp
) {
    LineReader lr(vcf_path);

    std::cout
        << "CHROM\tPOS\tREF\tALT\tVariationID\tAlleleID\tGeneSymbol\tClinicalSignificance\tOrigin\tPhenotype\t"
        << "HGVS_NUC\tHGVS_PROT\tAF\tCOSMIC_IDS\tPRIMARY_SITE\tHISTOLOGY\tCOUNT\t"
        << "TP53_CLASS\tTP53_SCORE\tTP53_NOTES\t"
        << "THERAPY_LEVEL\tDRUG\tDISEASE\tSOURCE\n";

    std::string line;
    while (lr.getline(line)) {
        if (line.empty() || line[0] == '#') continue;

        auto cols = split_tab(line);
        if ((int)cols.size() < 5) continue;

        std::string chrom = cols[0];
        std::string pos = cols[1];
        std::string ref = cols[3];
        std::string alt = cols[4];

        std::string key = make_coord_key(chrom, pos, ref, alt);

        auto it = clinvar_by_coord.find(key);
        if (it == clinvar_by_coord.end()) continue;

        const ClinVarRecord& rec = it->second;

        std::string hgvs_nuc, hgvs_prot, primary_hgvsp;

        auto hg_it = hgvs_by_varid.find(rec.variation_id);
        if (hg_it != hgvs_by_varid.end()) {
            hgvs_nuc = join_semicolon(hg_it->second.nuc);
            hgvs_prot = join_semicolon(hg_it->second.prot);
            if (!hg_it->second.prot.empty()) primary_hgvsp = hg_it->second.prot.front();
        }

        std::string gnomad_af;
        auto g_it = gnomad_by_coord.find(key);
        if (g_it != gnomad_by_coord.end()) gnomad_af = g_it->second.af;

        std::string cosmic_ids, cosmic_site, cosmic_hist, cosmic_count;
        auto c_it = cosmic_by_coord.find(key);
        if (c_it != cosmic_by_coord.end()) {
            cosmic_ids = c_it->second.ids;
            cosmic_site = c_it->second.primary_site;
            cosmic_hist = c_it->second.histology;
            cosmic_count = c_it->second.count;
        }

        std::string tp53_class, tp53_score, tp53_notes;
        if (!primary_hgvsp.empty()) {
            auto t_it = tp53_by_hgvsp.find(primary_hgvsp);
            if (t_it != tp53_by_hgvsp.end()) {
                tp53_class = t_it->second.classification;
                tp53_score = t_it->second.func_score;
                tp53_notes = t_it->second.notes;
            }
        }

        std::string th_level, th_drug, th_disease, th_source;
        if (!primary_hgvsp.empty()) {
            auto th_it = therapy_by_hgvsp.find(primary_hgvsp);
            if (th_it != therapy_by_hgvsp.end()) {
                th_level = th_it->second.evidence_level;
                th_drug = th_it->second.drug;
                th_disease = th_it->second.disease;
                th_source = th_it->second.source;
            }
        }

        std::cout
            << chrom << "\t"
            << pos << "\t"
            << ref << "\t"
            << alt << "\t"
            << rec.variation_id << "\t"
            << rec.allele_id << "\t"
            << rec.gene_symbol << "\t"
            << rec.clinical_significance << "\t"
            << rec.origin << "\t"
            << rec.phenotype_list << "\t"
            << hgvs_nuc << "\t"
            << hgvs_prot << "\t"
            << gnomad_af << "\t"
            << cosmic_ids << "\t"
            << cosmic_site << "\t"
            << cosmic_hist << "\t"
            << cosmic_count << "\t"
            << tp53_class << "\t"
            << tp53_score << "\t"
            << tp53_notes << "\t"
            << th_level << "\t"
            << th_drug << "\t"
            << th_disease << "\t"
            << th_source
            << "\n";
    }
}

// ==================== TRYB GLOBAL ====================

static void dump_global(const std::unordered_map<std::string, ClinVarRecord>& clinvar_by_coord,
    const std::unordered_map<std::string, HgvsInfo>& hgvs_by_varid)
{
    std::cout
        << "VariationID\tAlleleID\tGeneSymbol\tClinicalSignificance\tOrigin\tPhenotype\tChrom\tPos\tRef\tAlt\tHGVS_NUC\tHGVS_PROT\n";

    for (const auto& kv : clinvar_by_coord) {
        const ClinVarRecord& rec = kv.second;

        auto hg_it = hgvs_by_varid.find(rec.variation_id);
        std::string hgvs_nuc, hgvs_prot;
        if (hg_it != hgvs_by_varid.end()) {
            hgvs_nuc = join_semicolon(hg_it->second.nuc);
            hgvs_prot = join_semicolon(hg_it->second.prot);
        }

        std::cout
            << rec.variation_id << "\t"
            << rec.allele_id << "\t"
            << rec.gene_symbol << "\t"
            << rec.clinical_significance << "\t"
            << rec.origin << "\t"
            << rec.phenotype_list << "\t"
            << rec.chrom << "\t"
            << rec.pos << "\t"
            << rec.ref << "\t"
            << rec.alt << "\t"
            << hgvs_nuc << "\t"
            << hgvs_prot << "\n";
    }
}

// ==================== MAIN ====================

int main(int argc, char* argv[]) {
    if (argc < 6 || argc > 10) {
        std::cerr
            << "Usage:\n"
            << "  " << argv[0] << " sclc   sample.vcf[.gz] variant_summary.txt[.gz] hgvs4variation.txt[.gz] variation_allele.txt[.gz] [gnomad.tsv] [cosmic.tsv] [tp53db.tsv] [therapy.tsv]\n"
            << "  " << argv[0] << " global dummy.vcf      variant_summary.txt[.gz] hgvs4variation.txt[.gz] variation_allele.txt[.gz]\n";
        return 1;
    }

    try {
        std::string mode = argv[1];
        bool sclc_mode = (mode == "sclc");

        std::string vcf_path = argv[2];
        std::string variant_summary_path = argv[3];
        std::string hgvs4variation_path = argv[4];
        std::string variation_allele_path = argv[5];

        std::string gnomad_path = (argc > 6) ? argv[6] : "-";
        std::string cosmic_path = (argc > 7) ? argv[7] : "-";
        std::string tp53db_path = (argc > 8) ? argv[8] : "-";
        std::string therapy_path = (argc > 9) ? argv[9] : "-";

        auto varid_to_alleleid = load_variation_allele(variation_allele_path);

        auto clinvar_by_coord = load_clinvar_filtered(
            variant_summary_path,
            varid_to_alleleid,
            sclc_mode
        );

        std::unordered_set<std::string> wanted_variation_ids;
        wanted_variation_ids.reserve(clinvar_by_coord.size());
        for (const auto& kv : clinvar_by_coord) wanted_variation_ids.insert(kv.second.variation_id);

        auto hgvs_by_varid = load_hgvs_for_variants(hgvs4variation_path, wanted_variation_ids);

        auto gnomad_by_coord = load_gnomad(gnomad_path);
        auto cosmic_by_coord = load_cosmic(cosmic_path);
        auto tp53_by_hgvsp = load_tp53db(tp53db_path);
        auto therapy_by_hgvsp = load_therapy(therapy_path);

        if (sclc_mode) {
            analyze_vcf_sclc(
                vcf_path,
                clinvar_by_coord,
                hgvs_by_varid,
                gnomad_by_coord,
                cosmic_by_coord,
                tp53_by_hgvsp,
                therapy_by_hgvsp
            );
        }
        else {
            dump_global(clinvar_by_coord, hgvs_by_varid);
        }

    }
    catch (const std::exception& ex) {
        std::cerr << "EXCEPTION: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
