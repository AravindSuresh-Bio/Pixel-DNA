# Pixel-DNA: Lossless Genomic Storage & Diagnostic Engine

**Pixel-DNA** is a high-performance bio-compression tool that converts raw DNA sequencing data (FASTQ) into **Lossless Diagnostic Images (PNG)**.

Unlike traditional binary compression (GZIP) which creates opaque "black box" files, Pixel-DNA transforms genomic data into **visually inspectable assets**, allowing researchers to diagnose sequencing failures (flow cell issues, quality drops) instantly by looking at the file thumbnail.

## 🚀 Key Features

* **Hybrid Archival Format:**
    * **Visual Layer:** Stores DNA sequences and Quality Scores as image pixels (snake-scan layout).
    * **Metadata Sidecar:** Preserves original read headers in a lightweight `.hdr` file for bit-perfect restoration.
* **Compression Efficiency:**
    * Achieves **~3:1 Lossless Compression** (comparable to GZIP).
    * Reduces storage footprint by **~65%** while adding visualization capabilities.
* **Instant Diagnostics:**
    * **Bright/Vivid Image:** High-quality run.
    * **Faded/White Areas:** Sequencing quality drop-off.
    * **Black/Dark Areas:** No signal/N-bases.
* **Fail-Safe Integrity:**
    * In binary formats (GZIP), a single bit flip corrupts the entire file.
    * In Pixel-DNA, a bit flip is just a single "dead pixel"—the rest of the genome remains readable.

## 🛠️ Technical Stack

* **Core Logic:** Python 3.13 + NumPy (Vectorized for C-like speed).
* **Compression:** LZ77 (via PNG Deflate algorithm) optimized for 2D biological patterns.
* **GUI:** Tkinter (Native Windows Interface).
* **Mapping Logic:**
    * **Hue:** Base Identity (A=Red, C=Green, G=Blue, T=Yellow).
    * **Saturation:** Phred Quality Score (0-40).

## ⚡ Quick Start

### 1. Encode (FASTQ $\to$ PNG)
1.  Launch `PixelDNA_Pro.exe`.
2.  Click **"IMPORT FASTQ"**.
3.  Check **"Preserve Metadata"** for archival mode.
4.  *Result:* Generates a `.png` image and a `.hdr` sidecar file.

### 2. Restore (PNG $\to$ FASTQ)
1.  Click **"IMPORT PNG"**.
2.  Select the Pixel-DNA image.
3.  *Result:* The engine reads the image + sidecar and reconstructs the original FASTQ file byte-for-byte.

## 📊 Performance Benchmark

| Metric | Original (FASTQ) | Pixel-DNA (PNG+HDR) | Efficiency |
| :--- | :--- | :--- | :--- |
| **File Size** | 6.54 MB | 2.24 MB | **2.91x Compression** |
| **Integrity** | N/A | 100% Bit-Identical | Lossless |
| **Visual QC** | Impossible (Text) | **Instant** | < 1 Second |

---
*Developed by Aravind Suresh as a prototype for high-throughput genomic storage research.*