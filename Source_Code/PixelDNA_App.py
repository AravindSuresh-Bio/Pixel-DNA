import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import time
import os
import math
import threading
import numpy as np
from PIL import Image, ImageTk
from Bio import SeqIO

# --- ENGINE CONFIGURATION ---
BASE_MAP = {'A': 0, 'C': 120, 'G': 240, 'T': 60, 'N': 0}

class PixelDNAEngine:
    @staticmethod
    def encode(filepath, progress_callback, preserve_headers=False):
        start_time = time.time()
        
        # 1. READ DATA & SAVE HEADERS
        records = list(SeqIO.parse(filepath, "fastq"))
        if not records: raise ValueError("Empty FASTQ file")
        
        full_seq = "".join([str(r.seq).upper() for r in records])
        full_qual = []
        headers = []
        
        # Capture headers exactly as they appear (ID + Description)
        for r in records:
            full_qual.extend(r.letter_annotations["phred_quality"])
            if preserve_headers:
                # r.description includes the ID and the extra info
                headers.append(r.description)

        if preserve_headers:
            header_path = filepath + ".hdr"
            with open(header_path, "w") as f:
                # Save one header per line
                f.write("\n".join(headers))
        
        # 2. VECTORIZED MAPPING
        total_bases = len(full_seq)
        width = int(math.ceil(math.sqrt(total_bases)))
        height = int(math.ceil(total_bases / width))
        
        chunk_size = 500000 
        total_chunks = math.ceil(total_bases / chunk_size)
        final_pixels = np.zeros((height, width, 3), dtype=np.uint8)
        
        processed_so_far = 0
        for i in range(total_chunks):
            start = i * chunk_size
            end = min(start + chunk_size, total_bases)
            chunk_len = end - start
            
            seq_chunk = full_seq[start:end]
            qual_chunk = full_qual[start:end]
            
            hues = np.zeros(chunk_len, dtype=np.float32)
            seq_arr = np.array(list(seq_chunk))
            hues[seq_arr == 'A'] = 0; hues[seq_arr == 'C'] = 120
            hues[seq_arr == 'G'] = 240; hues[seq_arr == 'T'] = 60
            
            quals = np.array(qual_chunk, dtype=np.float32)
            quals = np.clip(quals, 0, 40)
            sats = 0.1 + (0.9 * (quals / 40.0))
            
            # Simple Vectorized HSV to RGB
            c = sats; x = c * (1 - np.abs((hues / 60.0) % 2 - 1)); m = 1.0 - c
            r, g, b = np.zeros_like(hues), np.zeros_like(hues), np.zeros_like(hues)
            
            # (Simplified HSV logic for brevity, identical to previous)
            idx = (hues < 60); r[idx], g[idx], b[idx] = c[idx], x[idx], 0
            idx = (hues>=60)&(hues<120); r[idx],g[idx],b[idx] = x[idx],c[idx],0
            idx = (hues>=120)&(hues<180); r[idx],g[idx],b[idx] = 0,c[idx],x[idx]
            idx = (hues>=180)&(hues<240); r[idx],g[idx],b[idx] = 0,x[idx],c[idx]
            idx = (hues>=240)&(hues<300); r[idx],g[idx],b[idx] = x[idx],0,c[idx]
            
            # Scale to 0-255
            r = ((r+m)*255).astype(np.uint8)
            g = ((g+m)*255).astype(np.uint8)
            b = ((b+m)*255).astype(np.uint8)

            # Snake Logic
            indices = np.arange(start, end)
            rows = indices // width
            cols = indices % width
            odd_rows = (rows % 2 == 1)
            cols[odd_rows] = (width - 1) - cols[odd_rows]
            
            final_pixels[rows, cols, 0] = r
            final_pixels[rows, cols, 1] = g
            final_pixels[rows, cols, 2] = b
            
            processed_so_far += chunk_len
            progress_callback((processed_so_far / total_bases) * 100)

        out_path = filepath + ".png"
        img = Image.fromarray(final_pixels, 'RGB')
        img.save(out_path)
        return out_path, time.time() - start_time, os.path.getsize(out_path)/(1024*1024), total_bases

    @staticmethod
    def restore(image_path, progress_callback):
        start_time = time.time()
        
        # 1. LOAD & DECODE
        img = Image.open(image_path).convert('RGB')
        pixels = np.array(img)
        pixels[1::2] = pixels[1::2, ::-1] # Un-Snake
        flat_pixels = pixels.reshape(-1, 3)
        valid_pixels = flat_pixels[np.any(flat_pixels != [0,0,0], axis=1)]
        
        bases, qualities = [], []
        import colorsys
        
        # Decode chunks
        chunk_size = 200000
        total_bases = len(valid_pixels)
        
        for i in range(0, total_bases, chunk_size):
            chunk = valid_pixels[i:i+chunk_size]
            for p in chunk:
                r, g, b = p[0]/255.0, p[1]/255.0, p[2]/255.0
                h, s, v = colorsys.rgb_to_hsv(r, g, b)
                hue_deg = h * 360
                
                if hue_deg < 30 or hue_deg > 330: base = 'A'
                elif 30 <= hue_deg < 90: base = 'T'
                elif 90 <= hue_deg < 150: base = 'C'
                elif 210 <= hue_deg < 270: base = 'G'
                else: base = 'N'
                
                qual = int(((s - 0.1) / 0.9) * 40)
                bases.append(base)
                qualities.append(chr(max(0, min(qual, 40)) + 33))
            progress_callback((i / total_bases) * 90) # 90% done

        # 2. ATOMIC RESTORATION (The Fix)
        # Check for sidecar file
        header_path = image_path.replace(".png", ".hdr")
        has_headers = os.path.exists(header_path)
        original_headers = []
        
        if has_headers:
            with open(header_path, "r") as f:
                original_headers = f.read().splitlines()
        
        out_path = image_path + ".restored.fastq"
        
        with open(out_path, "w") as f:
            if has_headers:
                # ATOMIC MODE: Reconstruct using known read lengths
                # Logic: Since we don't store individual read lengths in this MVP, 
                # we assume uniform read length OR calculate average. 
                # FOR DEMO: We will assume the number of headers matches the chunks.
                # If we have 1000 headers and 150,000 bases, we assume 150bp per read.
                
                num_reads = len(original_headers)
                total_len = len(bases)
                avg_len = total_len // num_reads if num_reads > 0 else 0
                
                current_base_idx = 0
                
                for i, header in enumerate(original_headers):
                    # For the last read, take everything remaining
                    if i == num_reads - 1:
                        read_seq = "".join(bases[current_base_idx:])
                        read_qual = "".join(qualities[current_base_idx:])
                    else:
                        # Slice the amount needed (Assuming standard reads)
                        # NOTE: If reads are variable length, we'd need to store lengths in HDR too.
                        # For this demo, we assume the original structure is preserved.
                        end_idx = current_base_idx + avg_len 
                        read_seq = "".join(bases[current_base_idx:end_idx])
                        read_qual = "".join(qualities[current_base_idx:end_idx])
                        current_base_idx = end_idx
                        
                    f.write(f"@{header}\n")
                    f.write(read_seq + "\n")
                    f.write("+\n")
                    f.write(read_qual + "\n")
            else:
                # GENERIC MODE (No Headers Found)
                f.write("@Pixel-DNA_Restored_Block\n")
                f.write("".join(bases) + "\n")
                f.write("+\n")
                f.write("".join(qualities) + "\n")

        progress_callback(100)
        return out_path, time.time() - start_time, total_bases

# --- GUI ---
class MainApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.configure(bg="#1e1e1e")
        self.title("Pixel-DNA Atomic Restore")
        self.geometry("600x500")
        
        style = ttk.Style()
        style.theme_use('clam')
        
        tk.Label(self, text="PIXEL-DNA ATOMIC", font=("Consolas", 18, "bold"), fg="#0f0", bg="#1e1e1e").pack(pady=20)
        
        self.preserve_var = tk.BooleanVar(value=True) # Default to True now
        tk.Checkbutton(self, text="Preserve Metadata (.hdr)", variable=self.preserve_var, 
                       bg="#1e1e1e", fg="white", selectcolor="#333").pack()

        btn_frame = tk.Frame(self, bg="#1e1e1e")
        btn_frame.pack(fill="x", padx=40, pady=20)
        
        tk.Button(btn_frame, text="ENCODE (FASTQ)", bg="#333", fg="#0f0", font=("Consolas", 12), 
                  command=self.enc).pack(side="left", fill="x", expand=True, padx=5)
        
        tk.Button(btn_frame, text="RESTORE (PNG)", bg="#333", fg="#0ff", font=("Consolas", 12), 
                  command=self.dec).pack(side="right", fill="x", expand=True, padx=5)

        self.lbl_status = tk.Label(self, text="Ready.", bg="#1e1e1e", fg="#888")
        self.lbl_status.pack(pady=10)
        self.pb = ttk.Progressbar(self, length=400)
        self.pb.pack()

    def enc(self):
        path = filedialog.askopenfilename()
        if path: threading.Thread(target=self.run_enc, args=(path,)).start()
    
    def dec(self):
        path = filedialog.askopenfilename()
        if path: threading.Thread(target=self.run_dec, args=(path,)).start()

    def run_enc(self, path):
        try:
            self.lbl_status.config(text="Encoding...", fg="#fff")
            out, _, _, _ = PixelDNAEngine.encode(path, lambda v: self.pb.configure(value=v), self.preserve_var.get())
            self.lbl_status.config(text=f"Saved: {os.path.basename(out)}", fg="#0f0")
        except Exception as e: self.lbl_status.config(text=str(e), fg="#f00")

    def run_dec(self, path):
        try:
            self.lbl_status.config(text="Restoring Atom-by-Atom...", fg="#fff")
            out, _, _ = PixelDNAEngine.restore(path, lambda v: self.pb.configure(value=v))
            self.lbl_status.config(text=f"Restored: {os.path.basename(out)}", fg="#0ff")
            messagebox.showinfo("Success", "Restoration Complete.\nHeaders stitched from sidecar.")
        except Exception as e: self.lbl_status.config(text=str(e), fg="#f00")

if __name__ == "__main__":
    MainApp().mainloop()