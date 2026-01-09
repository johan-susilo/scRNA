#!/bin/bash
###############################################################################
# Variant Calling Visualization Script
#
# Purpose: Generate comprehensive visualizations for trio variant analysis
#
# Visualizations:
#   1. Variant quality metrics (QC plots)
#   2. De novo variant distribution
#   3. Ti/Tv ratio analysis
#   4. Variant type distribution
#   5. Impact distribution (SnpEff)
#   6. Autism gene analysis
#   7. Coverage statistics
#   8. Summary dashboard
#
# All steps skip if output already exists
###############################################################################

set -euo pipefail

# Configuration
OUTPUT_DIR="${HOME}/output/autism"
RESULTS_DIR="${OUTPUT_DIR}/variants"
ANNOTATION_DIR="${OUTPUT_DIR}/comprehensive_analysis"
VIZ_DIR="${OUTPUT_DIR}/visualizations"
SCRIPTS_DIR="${VIZ_DIR}/scripts"

# Input files
COHORT_VCF="${RESULTS_DIR}/trio_joint_called.vcf.gz"
DENOVO_VCF="${RESULTS_DIR}/denovo_high_confidence.vcf.gz"
ANNOTATED_VCF="${ANNOTATION_DIR}/denovo_annotated_clean.vcf.gz"
LOF_FILE="${ANNOTATION_DIR}/lof_variants.txt"
AUTISM_VARS="${ANNOTATION_DIR}/autism_gene_details.txt"

echo "=========================================="
echo "Variant Analysis Visualization"
echo "=========================================="
echo "Date: $(date)"
echo "Output: ${VIZ_DIR}"
echo ""

# Activate conda environment
source "${HOME}/tool/anaconda3/etc/profile.d/conda.sh"
conda activate varcall

# Create directories
mkdir -p "${VIZ_DIR}"
mkdir -p "${SCRIPTS_DIR}"
cd "${VIZ_DIR}"

###############################################################################
# Check if visualization tools are available
###############################################################################

echo "=== Checking visualization tools ==="

# Check for R (optional but recommended)
if command -v R &> /dev/null; then
    echo "✓ R found"
    HAS_R=true
else
    echo "⚠ R not found - some plots will be text-based only"
    HAS_R=false
fi

# Check for required tools
for tool in bcftools python3; do
    if command -v $tool &> /dev/null; then
        echo "✓ $tool found"
    else
        echo "ERROR: $tool not found - required for visualization"
        exit 1
    fi
done

echo ""

###############################################################################
# STEP 1: Extract Variant Statistics
###############################################################################

echo "=== STEP 1: Extracting variant statistics ==="

STATS_FILE="${VIZ_DIR}/variant_stats.txt"

if [ -f "${STATS_FILE}" ]; then
    echo "✓ Variant statistics already exist - skipping"
else
    echo "Calculating statistics from VCFs..."

    {
        echo "=========================================="
        echo "VARIANT STATISTICS SUMMARY"
        echo "=========================================="
        echo ""
        echo "Generated: $(date)"
        echo ""

        # Total variants
        if [ -f "${COHORT_VCF}" ]; then
            TOTAL_VARS=$(bcftools view -H "${COHORT_VCF}" | wc -l)
            echo "Total variants called: ${TOTAL_VARS}"

            # SNPs vs INDELs
            SNPS=$(bcftools view -v snps -H "${COHORT_VCF}" | wc -l)
            INDELS=$(bcftools view -v indels -H "${COHORT_VCF}" | wc -l)
            echo "  SNPs: ${SNPS}"
            echo "  INDELs: ${INDELS}"
        fi

        echo ""

        # De novo variants
        if [ -f "${DENOVO_VCF}" ]; then
            DENOVO_COUNT=$(bcftools view -H "${DENOVO_VCF}" | wc -l)
            echo "De novo variants: ${DENOVO_COUNT}"

            DENOVO_SNPS=$(bcftools view -v snps -H "${DENOVO_VCF}" | wc -l)
            DENOVO_INDELS=$(bcftools view -v indels -H "${DENOVO_VCF}" | wc -l)
            echo "  SNPs: ${DENOVO_SNPS}"
            echo "  INDELs: ${DENOVO_INDELS}"
        fi

        echo ""

        # SnpEff annotation stats
        if [ -f "${ANNOTATED_VCF}" ]; then
            echo "SnpEff Annotation Summary:"

            HIGH=$(bcftools view -H "${ANNOTATED_VCF}" | grep -c "|HIGH|" || echo "0")
            MODERATE=$(bcftools view -H "${ANNOTATED_VCF}" | grep -c "|MODERATE|" || echo "0")
            LOW=$(bcftools view -H "${ANNOTATED_VCF}" | grep -c "|LOW|" || echo "0")
            MODIFIER=$(bcftools view -H "${ANNOTATED_VCF}" | grep -c "|MODIFIER|" || echo "0")

            echo "  HIGH impact: ${HIGH}"
            echo "  MODERATE impact: ${MODERATE}"
            echo "  LOW impact: ${LOW}"
            echo "  MODIFIER: ${MODIFIER}"
        fi

        echo ""

    } > "${STATS_FILE}"

    echo "✓ Statistics saved: ${STATS_FILE}"
fi

echo ""

###############################################################################
# STEP 2: Chromosome Distribution
###############################################################################

echo "=== STEP 2: Chromosome distribution ==="

CHR_DIST="${VIZ_DIR}/chromosome_distribution.txt"

if [ -f "${CHR_DIST}" ]; then
    echo "✓ Chromosome distribution already exists - skipping"
else
    echo "Calculating chromosome distribution..."

    {
        echo "=========================================="
        echo "DE NOVO VARIANTS BY CHROMOSOME"
        echo "=========================================="
        echo ""

        if [ -f "${DENOVO_VCF}" ]; then
            echo "Chromosome  Count  Percentage"
            echo "----------  -----  ----------"

            TOTAL=$(bcftools view -H "${DENOVO_VCF}" | wc -l)

            bcftools view -H "${DENOVO_VCF}" | cut -f1 | sort | uniq -c | \
            awk -v total="$TOTAL" '{
                chr=$2; count=$1
                pct=count*100/total
                printf "%-10s  %5d  %6.2f%%\n", chr, count, pct
            }' | sort -k1V
        fi

        echo ""
    } > "${CHR_DIST}"

    echo "✓ Chromosome distribution saved: ${CHR_DIST}"
fi

cat "${CHR_DIST}"
echo ""

###############################################################################
# STEP 3: Quality Metrics Distribution
###############################################################################

echo "=== STEP 3: Quality metrics distribution ==="

QUAL_DIST="${VIZ_DIR}/quality_distribution.txt"

if [ -f "${QUAL_DIST}" ]; then
    echo "✓ Quality distribution already exists - skipping"
else
    echo "Extracting quality metrics..."

    {
        echo "=========================================="
        echo "VARIANT QUALITY DISTRIBUTION"
        echo "=========================================="
        echo ""

        if [ -f "${DENOVO_VCF}" ]; then
            # Extract QUAL values
            bcftools query -f '%QUAL\n' "${DENOVO_VCF}" | \
            awk '
            BEGIN {
                q0_50=0; q50_100=0; q100_200=0; q200_500=0; q500_1000=0; q1000plus=0
            }
            {
                if ($1 < 50) q0_50++
                else if ($1 < 100) q50_100++
                else if ($1 < 200) q100_200++
                else if ($1 < 500) q200_500++
                else if ($1 < 1000) q500_1000++
                else q1000plus++
                sum += $1
                count++
                if (min == "" || $1 < min) min = $1
                if ($1 > max) max = $1
            }
            END {
                print "Quality Range    Count  Percentage"
                print "-------------    -----  ----------"
                printf "0-50            %6d  %6.2f%%\n", q0_50, q0_50*100/count
                printf "50-100          %6d  %6.2f%%\n", q50_100, q50_100*100/count
                printf "100-200         %6d  %6.2f%%\n", q100_200, q100_200*100/count
                printf "200-500         %6d  %6.2f%%\n", q200_500, q200_500*100/count
                printf "500-1000        %6d  %6.2f%%\n", q500_1000, q500_1000*100/count
                printf "1000+           %6d  %6.2f%%\n", q1000plus, q1000plus*100/count
                print ""
                printf "Min quality: %.2f\n", min
                printf "Max quality: %.2f\n", max
                printf "Mean quality: %.2f\n", sum/count
            }'
        fi

        echo ""
    } > "${QUAL_DIST}"

    echo "✓ Quality distribution saved: ${QUAL_DIST}"
fi

cat "${QUAL_DIST}"
echo ""

###############################################################################
# STEP 4: Ti/Tv Ratio Calculation
###############################################################################

echo "=== STEP 4: Ti/Tv ratio analysis ==="

TITV_FILE="${VIZ_DIR}/titv_ratio.txt"

if [ -f "${TITV_FILE}" ]; then
    echo "✓ Ti/Tv ratio already calculated - skipping"
else
    echo "Calculating Ti/Tv ratio..."

    {
        echo "=========================================="
        echo "TRANSITION / TRANSVERSION RATIO"
        echo "=========================================="
        echo ""

        if [ -f "${DENOVO_VCF}" ]; then
            bcftools view -v snps -H "${DENOVO_VCF}" | \
            awk '{
                ref=$4; alt=$5
                change=ref">"alt

                # Transitions (A<->G, C<->T)
                if (change=="A>G" || change=="G>A" || change=="C>T" || change=="T>C")
                    ti++
                # Transversions (all others)
                else
                    tv++
            }
            END {
                print "Transitions (Ti): " ti
                print "Transversions (Tv): " tv
                printf "Ti/Tv ratio: %.3f\n", ti/tv
                print ""
                print "Expected Ti/Tv ratio:"
                print "  Whole genome: ~2.0-2.1"
                print "  Exome: ~3.0-3.5"
                print ""
                if (ti/tv >= 1.5 && ti/tv <= 2.5)
                    print "✓ Ti/Tv ratio is within expected range (good quality)"
                else
                    print "⚠ Ti/Tv ratio is outside expected range (check quality)"
            }'
        fi

    } > "${TITV_FILE}"

    echo "✓ Ti/Tv ratio saved: ${TITV_FILE}"
fi

cat "${TITV_FILE}"
echo ""

###############################################################################
# STEP 5: Impact Distribution (SnpEff)
###############################################################################

echo "=== STEP 5: SnpEff impact distribution ==="

IMPACT_DIST="${VIZ_DIR}/impact_distribution.txt"

if [ -f "${IMPACT_DIST}" ]; then
    echo "✓ Impact distribution already exists - skipping"
else
    echo "Analyzing impact distribution..."

    if [ -f "${ANNOTATED_VCF}" ]; then
        {
            echo "=========================================="
            echo "SNPEFF IMPACT DISTRIBUTION"
            echo "=========================================="
            echo ""

            # Count impacts using grep
            TOTAL=$(bcftools view -H "${ANNOTATED_VCF}" | wc -l)
            HIGH=$(bcftools view -H "${ANNOTATED_VCF}" | grep -c "|HIGH|" || echo "0")
            MODERATE=$(bcftools view -H "${ANNOTATED_VCF}" | grep -c "|MODERATE|" || echo "0")
            LOW=$(bcftools view -H "${ANNOTATED_VCF}" | grep -c "|LOW|" || echo "0")
            MODIFIER=$(bcftools view -H "${ANNOTATED_VCF}" | grep -c "|MODIFIER|" || echo "0")

            echo "Impact Level    Count  Percentage"
            echo "------------    -----  ----------"
            printf "HIGH           %6d  %6.2f%%\n" $HIGH $(echo "scale=2; $HIGH*100/$TOTAL" | bc)
            printf "MODERATE       %6d  %6.2f%%\n" $MODERATE $(echo "scale=2; $MODERATE*100/$TOTAL" | bc)
            printf "LOW            %6d  %6.2f%%\n" $LOW $(echo "scale=2; $LOW*100/$TOTAL" | bc)
            printf "MODIFIER       %6d  %6.2f%%\n" $MODIFIER $(echo "scale=2; $MODIFIER*100/$TOTAL" | bc)
            echo ""
            printf "Total variants: %d\n" $TOTAL
            echo ""
            echo "Impact Categories:"
            echo "  HIGH     - Loss of function (stop-gain, frameshift, splice)"
            echo "  MODERATE - Missense variants (amino acid changes)"
            echo "  LOW      - Synonymous variants (no amino acid change)"
            echo "  MODIFIER - Non-coding variants"
            echo ""

        } > "${IMPACT_DIST}"

        echo "✓ Impact distribution saved: ${IMPACT_DIST}"
        cat "${IMPACT_DIST}"
    else
        echo "⚠ Annotated VCF not found - skipping impact distribution"
    fi
fi

echo ""

###############################################################################
# STEP 6: Effect Type Distribution
###############################################################################

echo "=== STEP 6: Variant effect types ==="

EFFECT_DIST="${VIZ_DIR}/effect_distribution.txt"

if [ -f "${EFFECT_DIST}" ]; then
    echo "✓ Effect distribution already exists - skipping"
else
    echo "Analyzing effect types..."

    if [ -f "${ANNOTATED_VCF}" ]; then
        {
            echo "=========================================="
            echo "TOP VARIANT EFFECT TYPES"
            echo "=========================================="
            echo ""

            bcftools view -H "${ANNOTATED_VCF}" | \
            awk '{
                for(i=1;i<=NF;i++){
                    if($i ~ /^ANN=/){
                        ann=$i
                        sub(/.*;ANN=/,"",ann)
                        split(ann, effects, ",")
                        for(j in effects){
                            split(effects[j], fields, "|")
                            effect=fields[2]
                            gsub(/&/, "\n", effect)
                            print effect
                        }
                        break
                    }
                }
            }' | sort | uniq -c | sort -rn | head -20 | \
            awk '{printf "%-40s %6d\n", $2, $1}'

            echo ""

        } > "${EFFECT_DIST}"

        echo "✓ Effect distribution saved: ${EFFECT_DIST}"
        cat "${EFFECT_DIST}"
    else
        echo "⚠ Annotated VCF not found - skipping effect distribution"
    fi
fi

echo ""

###############################################################################
# STEP 7: Autism Gene Analysis
###############################################################################

echo "=== STEP 7: Autism gene analysis ==="

AUTISM_SUMMARY="${VIZ_DIR}/autism_genes_summary.txt"

if [ -f "${AUTISM_SUMMARY}" ]; then
    echo "✓ Autism gene summary already exists - skipping"
else
    echo "Summarizing autism gene findings..."

    {
        echo "=========================================="
        echo "AUTISM GENE VARIANT SUMMARY"
        echo "=========================================="
        echo ""

        if [ -f "${AUTISM_VARS}" ] && [ $(wc -l < "${AUTISM_VARS}") -gt 1 ]; then
            # Count variants per gene
            echo "Variants per autism gene:"
            echo ""
            echo "Gene          Count  HIGH  MODERATE  LOW  MODIFIER"
            echo "----          -----  ----  --------  ---  --------"

            tail -n +2 "${AUTISM_VARS}" | \
            awk '{gene=$1; impact=$7; count[gene]++; impacts[gene,impact]++}
            END {
                for(gene in count){
                    printf "%-12s  %5d  %4d  %8d  %3d  %8d\n",
                        gene, count[gene],
                        impacts[gene,"HIGH"]+0,
                        impacts[gene,"MODERATE"]+0,
                        impacts[gene,"LOW"]+0,
                        impacts[gene,"MODIFIER"]+0
                }
            }' | sort -k2 -rn | head -20

            echo ""
            echo "(Showing top 20 genes by variant count)"
            echo ""

            # Summary
            TOTAL_GENES=$(tail -n +2 "${AUTISM_VARS}" | cut -f1 | sort -u | wc -l)
            TOTAL_VARS=$(tail -n +2 "${AUTISM_VARS}" | wc -l)

            echo "Total autism genes with variants: ${TOTAL_GENES}"
            echo "Total variant annotations in autism genes: ${TOTAL_VARS}"

        else
            echo "No autism gene variants found or file missing"
        fi

        echo ""

    } > "${AUTISM_SUMMARY}"

    echo "✓ Autism gene summary saved: ${AUTISM_SUMMARY}"
fi

cat "${AUTISM_SUMMARY}"
echo ""

###############################################################################
# STEP 8: Loss-of-Function Summary
###############################################################################

echo "=== STEP 8: Loss-of-function summary ==="

LOF_SUMMARY="${VIZ_DIR}/lof_summary.txt"

if [ -f "${LOF_SUMMARY}" ]; then
    echo "✓ LoF summary already exists - skipping"
else
    echo "Summarizing loss-of-function variants..."

    {
        echo "=========================================="
        echo "LOSS-OF-FUNCTION VARIANTS SUMMARY"
        echo "=========================================="
        echo ""

        if [ -f "${LOF_FILE}" ] && [ $(wc -l < "${LOF_FILE}") -gt 1 ]; then
            TOTAL_LOF=$(tail -n +2 "${LOF_FILE}" | wc -l)
            echo "Total LoF variants: ${TOTAL_LOF}"
            echo ""

            # Count by effect type
            echo "LoF by effect type:"
            echo ""
            tail -n +2 "${LOF_FILE}" | cut -f7 | \
            awk '{
                gsub(/&/, "\n", $0)
                print
            }' | sort | uniq -c | sort -rn | head -10 | \
            awk '{printf "  %-40s %5d\n", $2, $1}'

            echo ""

            # Top genes with LoF
            echo "Top 10 genes with LoF variants:"
            echo ""
            tail -n +2 "${LOF_FILE}" | cut -f1 | sort | uniq -c | sort -rn | head -10 | \
            awk '{printf "  %-20s %5d variants\n", $2, $1}'

        else
            echo "No LoF variants found or file missing"
        fi

        echo ""

    } > "${LOF_SUMMARY}"

    echo "✓ LoF summary saved: ${LOF_SUMMARY}"
fi

cat "${LOF_SUMMARY}"
echo ""

###############################################################################
# STEP 9: Create Python Visualization Script
###############################################################################

echo "=== STEP 9: Creating Python visualization script ==="

PYTHON_VIZ="${SCRIPTS_DIR}/plot_variants.py"

if [ -f "${PYTHON_VIZ}" ]; then
    echo "✓ Python script already exists - skipping"
else
    echo "Creating Python plotting script..."

cat > "${PYTHON_VIZ}" << 'EOFPYTHON'
#!/usr/bin/env python3
"""
Variant Analysis Visualization Script
Generates plots from variant calling results
"""

import sys
import os

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    print("⚠ matplotlib not available - skipping plots")
    HAS_MATPLOTLIB = False
    sys.exit(0)

# Set style
plt.style.use('seaborn-v0_8-darkgrid' if 'seaborn-v0_8-darkgrid' in plt.style.available else 'default')

def plot_impact_distribution(viz_dir):
    """Plot SnpEff impact distribution"""
    impact_file = os.path.join(viz_dir, 'impact_distribution.txt')

    if not os.path.exists(impact_file):
        return

    # Parse impact data
    impacts = []
    counts = []

    with open(impact_file) as f:
        for line in f:
            if any(imp in line for imp in ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']):
                parts = line.split()
                if len(parts) >= 2:
                    impacts.append(parts[0])
                    counts.append(int(parts[1]))

    if not counts:
        return

    # Create pie chart
    fig, ax = plt.subplots(figsize=(10, 8))
    colors = ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4']
    explode = (0.1, 0, 0, 0)  # Explode HIGH impact

    ax.pie(counts, labels=impacts, autopct='%1.1f%%', startangle=90,
           colors=colors, explode=explode, shadow=True)
    ax.set_title('SnpEff Impact Distribution', fontsize=16, fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(viz_dir, 'impact_distribution.png'), dpi=300)
    plt.close()

    print("  ✓ Created: impact_distribution.png")

def plot_chromosome_distribution(viz_dir):
    """Plot chromosome distribution"""
    chr_file = os.path.join(viz_dir, 'chromosome_distribution.txt')

    if not os.path.exists(chr_file):
        return

    # Parse chromosome data
    chromosomes = []
    counts = []

    with open(chr_file) as f:
        for line in f:
            if line.startswith('chr'):
                parts = line.split()
                if len(parts) >= 2:
                    chromosomes.append(parts[0])
                    counts.append(int(parts[1]))

    if not counts:
        return

    # Create bar plot
    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(chromosomes))
    bars = ax.bar(x, counts, color='steelblue', edgecolor='black', linewidth=0.5)

    ax.set_xlabel('Chromosome', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of De Novo Variants', fontsize=12, fontweight='bold')
    ax.set_title('De Novo Variant Distribution by Chromosome', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(chromosomes, rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3)

    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(viz_dir, 'chromosome_distribution.png'), dpi=300)
    plt.close()

    print("  ✓ Created: chromosome_distribution.png")

def plot_quality_distribution(viz_dir):
    """Plot quality score distribution"""
    qual_file = os.path.join(viz_dir, 'quality_distribution.txt')

    if not os.path.exists(qual_file):
        return

    # Parse quality data
    ranges = []
    counts = []

    with open(qual_file) as f:
        for line in f:
            if '-' in line or '+' in line:
                parts = line.split()
                if len(parts) >= 2 and parts[1].isdigit():
                    ranges.append(parts[0])
                    counts.append(int(parts[1]))

    if not counts:
        return

    # Create bar plot
    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(ranges))
    colors = ['#d62728' if i == 0 else '#ff7f0e' if i == 1 else '#2ca02c'
              for i in range(len(ranges))]
    bars = ax.bar(x, counts, color=colors, edgecolor='black', linewidth=0.5)

    ax.set_xlabel('Quality Score Range', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Variants', fontsize=12, fontweight='bold')
    ax.set_title('Variant Quality Score Distribution', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(ranges, rotation=0)
    ax.grid(axis='y', alpha=0.3)

    # Add value labels
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(viz_dir, 'quality_distribution.png'), dpi=300)
    plt.close()

    print("  ✓ Created: quality_distribution.png")

def main():
    viz_dir = sys.argv[1] if len(sys.argv) > 1 else '/home/johan/output/autism/visualizations'

    print("\nGenerating plots...")

    plot_impact_distribution(viz_dir)
    plot_chromosome_distribution(viz_dir)
    plot_quality_distribution(viz_dir)

    print("\n✓ All plots generated successfully")

if __name__ == '__main__':
    main()
EOFPYTHON

    chmod +x "${PYTHON_VIZ}"
    echo "✓ Python script created: ${PYTHON_VIZ}"
fi

echo ""

###############################################################################
# STEP 10: Generate Plots
###############################################################################

echo "=== STEP 10: Generating plots ==="

# Check if plots already exist
PLOT_EXISTS=false
for plot in impact_distribution.png chromosome_distribution.png quality_distribution.png; do
    if [ -f "${VIZ_DIR}/${plot}" ]; then
        PLOT_EXISTS=true
        break
    fi
done

if [ "${PLOT_EXISTS}" = true ]; then
    echo "✓ Plots already exist - skipping generation"
else
    # Check for matplotlib
    if python3 -c "import matplotlib" 2>/dev/null; then
        echo "Running Python visualization script..."
        python3 "${PYTHON_VIZ}" "${VIZ_DIR}"
    else
        echo "⚠ matplotlib not available - skipping plot generation"
        echo "  Install with: pip install matplotlib"
    fi
fi

echo ""

###############################################################################
# STEP 11: Create HTML Dashboard
###############################################################################

echo "=== STEP 11: Creating HTML dashboard ==="

DASHBOARD="${VIZ_DIR}/index.html"

if [ -f "${DASHBOARD}" ]; then
    echo "✓ Dashboard already exists - skipping"
else
    echo "Creating HTML dashboard..."

cat > "${DASHBOARD}" << 'EOFHTML'
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Autism Trio Variant Analysis - Dashboard</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
        }
        .container {
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.3);
            overflow: hidden;
        }
        header {
            background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }
        header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
        }
        header p {
            font-size: 1.1em;
            opacity: 0.9;
        }
        .nav-tabs {
            display: flex;
            background: #34495e;
            padding: 0 20px;
        }
        .nav-tab {
            padding: 15px 25px;
            color: white;
            cursor: pointer;
            border: none;
            background: none;
            font-size: 1em;
            transition: background 0.3s;
        }
        .nav-tab:hover {
            background: rgba(255,255,255,0.1);
        }
        .nav-tab.active {
            background: white;
            color: #2c3e50;
        }
        .content {
            padding: 30px;
        }
        .tab-pane {
            display: none;
        }
        .tab-pane.active {
            display: block;
        }
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }
        .stat-card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        .stat-card h3 {
            font-size: 0.9em;
            opacity: 0.9;
            margin-bottom: 10px;
        }
        .stat-card .value {
            font-size: 2.5em;
            font-weight: bold;
        }
        .plot-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
            gap: 30px;
            margin: 30px 0;
        }
        .plot-container {
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .plot-container h3 {
            color: #2c3e50;
            margin-bottom: 15px;
            padding-bottom: 10px;
            border-bottom: 3px solid #667eea;
        }
        .plot-container img {
            width: 100%;
            height: auto;
            border-radius: 4px;
        }
        pre {
            background: #2c3e50;
            color: #ecf0f1;
            padding: 20px;
            border-radius: 8px;
            overflow-x: auto;
            font-size: 0.9em;
            line-height: 1.5;
        }
        .file-links {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }
        .file-link {
            background: #ecf0f1;
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
            transition: transform 0.2s, box-shadow 0.2s;
        }
        .file-link:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        }
        .file-link a {
            color: #2c3e50;
            text-decoration: none;
            font-weight: 500;
        }
        .file-link p {
            font-size: 0.85em;
            color: #7f8c8d;
            margin-top: 5px;
        }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>🧬 Autism Trio Variant Analysis</h1>
            <p>Comprehensive Visualization Dashboard</p>
            <p style="font-size: 0.9em; margin-top: 10px;">Generated: <span id="datetime"></span></p>
        </header>

        <div class="nav-tabs">
            <button class="nav-tab active" onclick="showTab('overview')">Overview</button>
            <button class="nav-tab" onclick="showTab('plots')">Visualizations</button>
            <button class="nav-tab" onclick="showTab('details')">Detailed Reports</button>
            <button class="nav-tab" onclick="showTab('files')">Output Files</button>
        </div>

        <div class="content">
            <!-- Overview Tab -->
            <div id="overview" class="tab-pane active">
                <h2 style="color: #2c3e50; margin-bottom: 20px;">Analysis Summary</h2>

                <div class="stats-grid">
                    <div class="stat-card">
                        <h3>Total Variants</h3>
                        <div class="value" id="total-vars">-</div>
                    </div>
                    <div class="stat-card">
                        <h3>De Novo Variants</h3>
                        <div class="value" id="denovo-vars">-</div>
                    </div>
                    <div class="stat-card">
                        <h3>HIGH Impact</h3>
                        <div class="value" id="high-impact">-</div>
                    </div>
                    <div class="stat-card">
                        <h3>Ti/Tv Ratio</h3>
                        <div class="value" id="titv-ratio">-</div>
                    </div>
                </div>

                <h3 style="color: #2c3e50; margin: 30px 0 15px 0;">Key Statistics</h3>
                <pre id="stats-content">Loading...</pre>
            </div>

            <!-- Plots Tab -->
            <div id="plots" class="tab-pane">
                <h2 style="color: #2c3e50; margin-bottom: 20px;">Visual Analysis</h2>

                <div class="plot-grid">
                    <div class="plot-container">
                        <h3>📊 Impact Distribution</h3>
                        <img src="impact_distribution.png" alt="Impact Distribution" onerror="this.src='data:image/svg+xml,<svg xmlns=%22http://www.w3.org/2000/svg%22><text y=%2250%%22 x=%2250%%22 text-anchor=%22middle%22>Plot not generated</text></svg>'">
                    </div>
                    <div class="plot-container">
                        <h3>📊 Chromosome Distribution</h3>
                        <img src="chromosome_distribution.png" alt="Chromosome Distribution" onerror="this.src='data:image/svg+xml,<svg xmlns=%22http://www.w3.org/2000/svg%22><text y=%2250%%22 x=%2250%%22 text-anchor=%22middle%22>Plot not generated</text></svg>'">
                    </div>
                    <div class="plot-container">
                        <h3>📊 Quality Distribution</h3>
                        <img src="quality_distribution.png" alt="Quality Distribution" onerror="this.src='data:image/svg+xml,<svg xmlns=%22http://www.w3.org/2000/svg%22><text y=%2250%%22 x=%2250%%22 text-anchor=%22middle%22>Plot not generated</text></svg>'">
                    </div>
                </div>
            </div>

            <!-- Details Tab -->
            <div id="details" class="tab-pane">
                <h2 style="color: #2c3e50; margin-bottom: 20px;">Detailed Reports</h2>

                <h3 style="color: #2c3e50; margin-top: 20px;">Chromosome Distribution</h3>
                <pre id="chr-dist">Loading...</pre>

                <h3 style="color: #2c3e50; margin-top: 30px;">Quality Metrics</h3>
                <pre id="qual-dist">Loading...</pre>

                <h3 style="color: #2c3e50; margin-top: 30px;">Ti/Tv Ratio</h3>
                <pre id="titv">Loading...</pre>

                <h3 style="color: #2c3e50; margin-top: 30px;">Autism Genes</h3>
                <pre id="autism-summary">Loading...</pre>
            </div>

            <!-- Files Tab -->
            <div id="files" class="tab-pane">
                <h2 style="color: #2c3e50; margin-bottom: 20px;">Output Files</h2>

                <div class="file-links">
                    <div class="file-link">
                        <a href="../TRUE_TOP_CANDIDATES.md">⭐ TRUE_TOP_CANDIDATES.md</a>
                        <p>Prioritized variant candidates for review</p>
                    </div>
                    <div class="file-link">
                        <a href="../comprehensive_analysis/denovo_snpeff_stats.html">📊 SnpEff HTML Report</a>
                        <p>Interactive SnpEff visualization</p>
                    </div>
                    <div class="file-link">
                        <a href="../comprehensive_analysis/lof_variants.txt">🧬 Loss-of-Function Variants</a>
                        <p>HIGH impact variants list</p>
                    </div>
                    <div class="file-link">
                        <a href="../comprehensive_analysis/autism_gene_details.txt">🎯 Autism Gene Variants</a>
                        <p>Variants in known autism genes</p>
                    </div>
                    <div class="file-link">
                        <a href="../variants/denovo_high_confidence.vcf.gz">📁 De Novo VCF</a>
                        <p>High-confidence de novo variants</p>
                    </div>
                    <div class="file-link">
                        <a href="../comprehensive_analysis/denovo_annotated_clean.vcf.gz">📁 Annotated VCF</a>
                        <p>SnpEff annotated variants (clean)</p>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <script>
        // Set datetime
        document.getElementById('datetime').textContent = new Date().toLocaleString();

        // Tab switching
        function showTab(tabName) {
            // Hide all tabs
            document.querySelectorAll('.tab-pane').forEach(pane => {
                pane.classList.remove('active');
            });
            document.querySelectorAll('.nav-tab').forEach(tab => {
                tab.classList.remove('active');
            });

            // Show selected tab
            document.getElementById(tabName).classList.add('active');
            event.target.classList.add('active');
        }

        // Load text files
        function loadFile(filename, elementId) {
            fetch(filename)
                .then(response => response.text())
                .then(data => {
                    document.getElementById(elementId).textContent = data;

                    // Parse stats for cards
                    if (elementId === 'stats-content') {
                        parseStats(data);
                    }
                })
                .catch(err => {
                    document.getElementById(elementId).textContent = 'File not found or error loading';
                });
        }

        function parseStats(data) {
            // Extract key metrics
            const totalMatch = data.match(/Total variants called:\s*(\d+)/);
            const denovoMatch = data.match(/De novo variants:\s*(\d+)/);
            const highMatch = data.match(/HIGH impact:\s*(\d+)/);
            const titvMatch = data.match(/Ti\/Tv ratio:\s*([\d.]+)/);

            if (totalMatch) document.getElementById('total-vars').textContent = parseInt(totalMatch[1]).toLocaleString();
            if (denovoMatch) document.getElementById('denovo-vars').textContent = parseInt(denovoMatch[1]).toLocaleString();
            if (highMatch) document.getElementById('high-impact').textContent = parseInt(highMatch[1]).toLocaleString();
            if (titvMatch) document.getElementById('titv-ratio').textContent = parseFloat(titvMatch[1]).toFixed(2);
        }

        // Load all text files
        loadFile('variant_stats.txt', 'stats-content');
        loadFile('chromosome_distribution.txt', 'chr-dist');
        loadFile('quality_distribution.txt', 'qual-dist');
        loadFile('titv_ratio.txt', 'titv');
        loadFile('autism_genes_summary.txt', 'autism-summary');
    </script>
</body>
</html>
EOFHTML

    echo "✓ HTML dashboard created: ${DASHBOARD}"
fi

echo ""

###############################################################################
# STEP 12: Create Summary Report
###############################################################################

echo "=== STEP 12: Creating visualization summary ==="

VIZ_REPORT="${VIZ_DIR}/VISUALIZATION_REPORT.txt"

if [ -f "${VIZ_REPORT}" ]; then
    echo "✓ Visualization report already exists - skipping"
else
    {
        echo "=========================================="
        echo "VISUALIZATION SUMMARY REPORT"
        echo "=========================================="
        echo ""
        echo "Generated: $(date)"
        echo "Location: ${VIZ_DIR}"
        echo ""
        echo "=========================================="
        echo "TEXT REPORTS"
        echo "=========================================="
        echo ""
        echo "✓ variant_stats.txt - Overall statistics"
        echo "✓ chromosome_distribution.txt - Variants by chromosome"
        echo "✓ quality_distribution.txt - Quality score analysis"
        echo "✓ titv_ratio.txt - Transition/transversion ratio"
        echo "✓ impact_distribution.txt - SnpEff impact categories"
        echo "✓ effect_distribution.txt - Variant effect types"
        echo "✓ autism_genes_summary.txt - Autism gene analysis"
        echo "✓ lof_summary.txt - Loss-of-function summary"
        echo ""
        echo "=========================================="
        echo "VISUALIZATIONS"
        echo "=========================================="
        echo ""

        if [ -f "${VIZ_DIR}/impact_distribution.png" ]; then
            echo "✓ impact_distribution.png - Impact pie chart"
            echo "✓ chromosome_distribution.png - Chromosome bar plot"
            echo "✓ quality_distribution.png - Quality histogram"
        else
            echo "⚠ PNG plots not generated (matplotlib not available)"
        fi

        echo ""
        echo "=========================================="
        echo "INTERACTIVE DASHBOARD"
        echo "=========================================="
        echo ""
        echo "✓ index.html - Main HTML dashboard"
        echo ""
        echo "Open in browser:"
        echo "  firefox ${DASHBOARD}"
        echo ""
        echo "Or copy to local machine:"
        echo "  scp -r user@server:${VIZ_DIR} ."
        echo ""
        echo "=========================================="
        echo "NEXT STEPS"
        echo "=========================================="
        echo ""
        echo "1. Open HTML dashboard in web browser"
        echo "2. Review text reports for detailed statistics"
        echo "3. Check plots for visual patterns"
        echo "4. Focus on HIGH impact variants"
        echo "5. Review autism gene findings"
        echo ""

    } > "${VIZ_REPORT}"

    echo "✓ Report created: ${VIZ_REPORT}"
fi

cat "${VIZ_REPORT}"

###############################################################################
# Final Summary
###############################################################################

echo ""
echo "=========================================="
echo "VISUALIZATION COMPLETE!"
echo "=========================================="
echo ""
echo "📊 Output location: ${VIZ_DIR}"
echo ""
echo "⭐ KEY FILES:"
echo ""
echo "  1. HTML Dashboard (open this first!)"
echo "     ${DASHBOARD}"
echo ""
echo "  2. Text Reports"
echo "     ${VIZ_DIR}/*.txt"
echo ""
echo "  3. Plots (if generated)"
echo "     ${VIZ_DIR}/*.png"
echo ""
echo "To view dashboard:"
echo "  firefox ${DASHBOARD}"
echo ""
echo "All visualizations use smart file skipping!"
echo "Re-run anytime - existing files won't be regenerated."
echo ""
echo "=========================================="
echo "Completed: $(date)"
echo "=========================================="
echo ""
