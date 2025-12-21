# 🔧 Fix for SRA Download Error

## Error You're Seeing
```
prefetch.3.2.0 int: name not found while resolving tree within virtual file system module
cannot get remote location for 'CM000663.1'
```

**Cause**: This is a known issue with SRA Toolkit's prefetch. It's trying to resolve chromosome references instead of downloading the actual data.

---

## ✅ Solution: Use Robust Download Script

I've created a script that tries **4 different methods** to download your trio data:

```bash
cd /home/johan/pipeline/varcall

# Run the robust downloader
./download_trio_robust.sh 2>&1 | tee download.log
```

---

## 🔄 What the Script Does (4 Fallback Methods)

### **Method 1**: Direct fasterq-dump (RECOMMENDED)
- Bypasses prefetch entirely
- Downloads and converts to FASTQ in one step
- Usually works when prefetch fails

### **Method 2**: SRA Toolkit with config fix
- Reconfigures SRA Toolkit settings
- Uses local cache
- Retries with corrected paths

### **Method 3**: ENA Mirror (European Alternative)
- Downloads from ENA FTP instead of NCBI
- Often more reliable outside USA
- Gets pre-converted FASTQ files

### **Method 4**: AWS S3 (if AWS CLI available)
- Downloads from Amazon cloud
- Part of NCBI's cloud distribution
- Very fast if you have AWS access

---

## 🎯 Quick Fix (Alternative Manual Method)

If the script still fails, download manually from ENA:

### Proband (SRR8697636):
```bash
cd ~/output/autism

wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR869/006/SRR8697636/SRR8697636_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR869/006/SRR8697636/SRR8697636_2.fastq.gz

mv SRR8697636_1.fastq.gz proband_1.fastq.gz
mv SRR8697636_2.fastq.gz proband_2.fastq.gz
```

### Mother (SRR8697627):
```bash
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR869/007/SRR8697627/SRR8697627_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR869/007/SRR8697627/SRR8697627_2.fastq.gz

mv SRR8697627_1.fastq.gz mother_1.fastq.gz
mv SRR8697627_2.fastq.gz mother_2.fastq.gz
```

### Father (SRR8697645):
```bash
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR869/005/SRR8697645/SRR8697645_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR869/005/SRR8697645/SRR8697645_2.fastq.gz

mv SRR8697645_1.fastq.gz father_1.fastq.gz
mv SRR8697645_2.fastq.gz father_2.fastq.gz
```

---

## 🔍 Alternative: Use SRA Run Selector

1. Go to: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRR8697636
2. Click "Metadata" tab
3. Select all 3 samples (SRR8697636, SRR8697627, SRR8697645)
4. Click "Download" → Choose "FASTQ files from ENA"
5. This gives you direct download links

---

## 🛠️ Fix SRA Toolkit Configuration (if needed)

If you want to fix prefetch permanently:

```bash
# Reset SRA Toolkit config
vdb-config --restore-defaults

# Enable cloud access
vdb-config --set /repository/user/main/public/root="${HOME}/output/autism/sra_cache"

# Disable prefetch to current directory (sometimes helps)
vdb-config --prefetch-to-cwd

# Test with one sample
fasterq-dump --split-files -X 10000 SRR8697636
# (downloads only first 10k reads as test)
```

---

## 📊 Check What You Already Have

```bash
# Check if SRA files exist
ls -lh ~/output/autism/SRR*/

# Check if FASTQ files exist
ls -lh ~/output/autism/*.fastq.gz

# Check disk space
df -h ~
```

---

## ⚡ Skip Download If Proband Already Downloaded

Since you already have `SRR8697636.sra` (2.4GB), you can convert it:

```bash
cd ~/output/autism

# Convert the existing SRA to FASTQ
fasterq-dump \
    --split-files \
    --threads 60 \
    --progress \
    --outdir . \
    SRR8697636/SRR8697636.sra

# Rename and compress
mv SRR8697636_1.fastq proband_1.fastq
mv SRR8697636_2.fastq proband_2.fastq
pigz -p 60 proband_1.fastq &
pigz -p 60 proband_2.fastq &
wait

echo "✓ Proband FASTQ ready!"
```

Then only download the parents (mother & father) using the robust script.

---

## 🎯 Recommended Workflow

```bash
# Step 1: Try the robust downloader
cd /home/johan/pipeline/varcall
./download_trio_robust.sh 2>&1 | tee download.log

# Step 2: If that fails, use manual ENA download (see above)

# Step 3: Verify files
ls -lh ~/output/autism/*.fastq.gz

# Step 4: Proceed with preprocessing
./preprocess_trio_optimized.sh 2>&1 | tee preprocess.log
```

---

## 💡 Why This Happens

The SRA Toolkit error occurs because:
1. **NCBI changed their infrastructure** (moved to cloud)
2. **Old toolkit versions** don't handle new paths
3. **Network/firewall** blocking NCBI but not ENA
4. **VDB (Virtual Database)** corruption in cache

**Solution**: Use direct download or alternative mirrors (ENA is most reliable).

---

## ✅ Success Indicators

You'll know it's working when you see:
```
spots read      : 50,000,000
reads read      : 100,000,000
reads written   : 100,000,000
```

And files created:
```
-rw-rw-r-- proband_1.fastq.gz  (5-8 GB)
-rw-rw-r-- proband_2.fastq.gz  (5-8 GB)
-rw-rw-r-- mother_1.fastq.gz   (5-8 GB)
-rw-rw-r-- mother_2.fastq.gz   (5-8 GB)
-rw-rw-r-- father_1.fastq.gz   (5-8 GB)
-rw-rw-r-- father_2.fastq.gz   (5-8 GB)
```

---

**Try the robust script first - it handles all these issues automatically!**
