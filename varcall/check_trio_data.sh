#!/bin/bash
#######################################################################
# Check Trio Data Availability
# This script helps you verify if you have all trio samples
#######################################################################

echo "=========================================="
echo "Trio Data Checker"
echo "=========================================="
echo ""

# Check proband
PROBAND_SRA="SRR8697636"
echo "Checking PROBAND: ${PROBAND_SRA}"
if [ -f "/home/johan/output/autism/${PROBAND_SRA}/${PROBAND_SRA}.sra" ]; then
    echo "  ✓ SRA file found"
    ls -lh "/home/johan/output/autism/${PROBAND_SRA}/${PROBAND_SRA}.sra"
else
    echo "  ✗ SRA file NOT found"
fi

echo ""
echo "Fetching metadata from NCBI for ${PROBAND_SRA}..."
echo "(This will show the BioProject and related samples)"
echo ""

# Try to get SRA metadata
if command -v esearch &> /dev/null; then
    esearch -db sra -query "${PROBAND_SRA}" | efetch -format runinfo | head -20
else
    echo "esearch not found. Install NCBI E-utilities to fetch metadata."
    echo ""
    echo "Manual check:"
    echo "  1. Visit: https://www.ncbi.nlm.nih.gov/sra/${PROBAND_SRA}"
    echo "  2. Look for 'BioProject' link"
    echo "  3. Find father and mother samples in the same project"
fi

echo ""
echo "=========================================="
echo "What you need:"
echo "=========================================="
echo "For de novo analysis, you need THREE samples:"
echo "  1. Proband (affected): ${PROBAND_SRA} [FOUND]"
echo "  2. Father (unaffected): SRR??????? [NEEDED]"
echo "  3. Mother (unaffected): SRR??????? [NEEDED]"
echo ""
echo "Next steps:"
echo "  1. Find parent SRA IDs from the BioProject"
echo "  2. Edit: /home/johan/pipeline/varcall/denovo_fam92_analysis.sh"
echo "  3. Update FATHER_SRA and MOTHER_SRA variables"
echo "  4. Run: ./denovo_fam92_analysis.sh"
echo ""
