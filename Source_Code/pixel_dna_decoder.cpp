#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>

char get_base_from_hue(float h) {
    if (h < 30 || h > 330) return 'A';
    if (30 <= h && h < 90) return 'T';
    if (90 <= h && h < 150) return 'C';
    if (210 <= h && h < 270) return 'G';
    return 'N';
}

char get_qual_from_sat(float s) {
    float val = (s - 0.1f) / 0.9f * 40.0f;
    int score = (int)(val + 0.5f);
    if (score < 0) score = 0; if (score > 40) score = 40;
    return (char)(score + 33);
}

void rgb_to_hsv(unsigned char r, unsigned char g, unsigned char b, float &h, float &s, float &v) {
    float rf = r/255.0f, gf = g/255.0f, bf = b/255.0f;
    float cmax = std::max({rf, gf, bf}), cmin = std::min({rf, gf, bf});
    float delta = cmax - cmin;
    if (delta == 0) h = 0;
    else if (cmax == rf) h = 60 * std::fmod(((gf-bf)/delta), 6.0f);
    else if (cmax == gf) h = 60 * (((bf-rf)/delta) + 2);
    else h = 60 * (((rf-gf)/delta) + 4);
    if (h < 0) h += 360;
    s = (cmax == 0) ? 0 : delta/cmax;
    v = cmax;
}

int main(int argc, char* argv[]) {
    if (argc < 2) { std::cout << "Usage: ./pixel_dna_decoder.exe input.png\n"; return 1; }

    std::string input_path = argv[1];
    std::string output_path = input_path + ".restored.fastq";
    
    int width, height, channels;
    unsigned char *img = stbi_load(input_path.c_str(), &width, &height, &channels, 3);
    
    if (img == NULL) {
        std::cerr << "CRITICAL ERROR: Could not load image.\n";
        std::cerr << "Reason: " << stbi_failure_reason() << "\n";
        std::cerr << "Attempted Path: [" << input_path << "]\n";
        return 1;
    }
    
    std::string seq_buffer, qual_buffer;
    long total = width * height;
    seq_buffer.reserve(total); qual_buffer.reserve(total);
    
    for (long i = 0; i < total; ++i) {
        int row = i / width;
        int col = i % width;
        if (row % 2 == 1) col = (width - 1) - col; // Snake Logic
        
        int idx = (row * width + col) * 3;
        unsigned char r = img[idx], g = img[idx+1], b = img[idx+2];
        if (r==0 && g==0 && b==0) continue; // Skip padding
        
        float h, s, v;
        rgb_to_hsv(r, g, b, h, s, v);
        seq_buffer += get_base_from_hue(h);
        qual_buffer += get_qual_from_sat(s);
    }
    stbi_image_free(img);
    
    std::ofstream outfile(output_path);
    outfile << "@Restored_Pixel_DNA\n" << seq_buffer << "\n+\n" << qual_buffer << "\n";
    outfile.close();
    
    std::cout << "SUCCESS: Restored to [" << output_path << "]\n";
    return 0;
}