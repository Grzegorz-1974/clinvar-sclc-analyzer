# clinvar-sclc-analyzer
# clinvar-tool


C++17 tool for ClinVar-driven variant analysis (VCF/VCF.gz) with HGVS mapping and optional gnomAD, COSMIC, TP53 functional and therapy annotations.

The tool supports a dedicated **SCLC mode** focused on somatic pathogenic variants in **TP53/RB1**, as well as a **global mode** for exporting ClinVar + HGVS tables.

---

## Features

- Parses VCF and VCF.gz files
- Integrates ClinVar:
  - `variant_summary`
  - `hgvs4variation`
  - `variation_allele`
- HGVS nucleotide and protein mapping
- Optional annotations:
  - gnomAD allele frequency
  - COSMIC cancer occurrence data
  - TP53 functional database
  - Therapy evidence (OncoKB / CIViC–style TSV)
- Two operating modes:
  - `sclc` – TP53/RB1, somatic, pathogenic(/likely pathogenic), lung/SCLC phenotype
  - `global` – ClinVar + HGVS export without phenotype filtering
- Plain TSV output suitable for downstream analysis

---

## Requirements

Tested on Linux (Fedora).

### System dependencies
```bash
sudo dnf groupinstall "Development Tools" -y
sudo dnf install -y cmake ninja-build gdb zlib-devel
