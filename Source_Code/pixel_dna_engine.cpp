#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <algorithm>

// --- CONFIGURATION ---
const int H_A = 0;   // Red
const int H_C = 120; // Green
const int H_G = 240; // Blue
const int H_T = 60;  // Yellow
const int H_N = 0;   // Black

struct Pixel { unsigned char r, g, b; };

Pixel hsv_to_rgb(float h, float s, float v) {
    float c = v * s;
    float x = c * (1 - std::abs(std::fmod(h / 60.0f, 2.0f) - 1));
    float m = v - c;
    float r_p, g_p, b_p;
    if (0 <= h && h < 60) { r_p = c; g_p = x; b_p = 0; }
    else if (60 <= h && h < 120) { r_p = x; g_p = c; b_p = 0; }
    else if (120 <= h && h < 180) { r_p = 0; g_p = c; b_p = x; }
    else if (180 <= h && h < 240) { r_p = 0; g_p = x; b_p = c; }
    else if (240 <= h && h < 300) { r_p = x; g_p = 0; b_p = c; }
    else { r_p = c; g_p = 0; b_p = x; }
    Pixel p;
    p.r = (unsigned char)((r_p + m) * 255);
    p.g = (unsigned char)((g_p + m) * 255);
    p.b = (unsigned char)((b_p + m) * 255);
    return p;
}

float qual_to_sat(char q_char) {
    int score = q_char - 33;
    if (score < 0) score = 0; if (score > 40) score = 40;
    return 0.1f + (0.9f * (score / 40.0f));
}

int get_hue(char base) {
    switch (base) {
        case 'A': return H_A; case 'C': return H_C;
        case 'G': return H_G; case 'T': return H_T;
        default: return H_N;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) { std::cout << "Usage: ./pixel_dna.exe input.fastq\n"; return 1; }

    std::string input_path = argv[1];
    std::string output_path = input_path + ".png";

    std::cout << "Loading: [" << input_path << "] ...\n";
    std::ifstream file(input_path);
    if (!file.is_open()) { std::cerr << "CRITICAL ERROR: Could not open input file.\n"; return 1; }

    std::string full_seq, full_qual;
    // Reserve memory to prevent reallocations
    full_seq.reserve(50000000); 
    full_qual.reserve(50000000);

    std::string line;
    int line_idx = 0;
    while (std::getline(file, line)) {
        if (line_idx % 4 == 1) full_seq += line;
        else if (line_idx % 4 == 3) full_qual += line;
        line_idx++;
    }
    file.close();

    size_t total_bases = full_seq.length();
    if (total_bases == 0) { std::cerr << "Error: File was empty.\n"; return 1; }

    int width = std::ceil(std::sqrt(total_bases));
    int height = std::ceil((double)total_bases / width);
    std::cout << "Grid: " << width << "x" << height << " (" << total_bases << " bases)\n";

    std::vector<unsigned char> img_data(width * height * 3, 0);

    for (size_t i = 0; i < total_bases; ++i) {
        size_t row = i / width;
        size_t col = i % width;
        if (row % 2 == 1) col = (width - 1) - col; // Snake Logic

        char base = full_seq[i];
        if (base >= 'a' && base <= 'z') base -= 32;

        int hue = get_hue(base);
        float sat = qual_to_sat(full_qual[i]);
        
        Pixel p = {0,0,0};
        if (hue != H_N || base == 'A') p = hsv_to_rgb((float)hue, sat, 1.0f);

        size_t idx = (row * width + col) * 3;
        img_data[idx] = p.r; img_data[idx+1] = p.g; img_data[idx+2] = p.b;
    }

    std::cout << "Writing PNG to [" << output_path << "] ...\n";
    int result = stbi_write_png(output_path.c_str(), width, height, 3, img_data.data(), width * 3);

    if (result == 0) {
        std::cerr << "CRITICAL ERROR: Failed to write PNG file.\n";
        std::cerr << "Check write permissions or disk space.\n";
        return 1;
    }
    std::cout << "SUCCESS: Image saved.\n";
    return 0;
}