#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <cstdint>
#include <limits>
#include <random>
#include <unordered_map>
#include <iomanip> 

uint32_t murmur_hash(const void * key, int len, uint32_t seed) {
    const uint8_t * data = (const uint8_t*)key;
    const int nblocks = len / 4;
    uint32_t h1 = seed;
    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;
    const uint32_t * blocks = (const uint32_t *)(data + nblocks*4);
    for(int i = -nblocks; i; i++) {
        uint32_t k1 = blocks[i]; k1 *= c1; k1 = (k1 << 15) | (k1 >> (32 - 15)); k1 *= c2;
        h1 ^= k1; h1 = (h1 << 13) | (h1 >> (32 - 13)); h1 = h1*5+0xe6546b64;
    }
    const uint8_t * tail = (const uint8_t*)(data + nblocks*4);
    uint32_t k1 = 0;
    switch(len & 3) {
        case 3: k1 ^= tail[2] << 16; case 2: k1 ^= tail[1] << 8; case 1: k1 ^= tail[0];
                k1 *= c1; k1 = (k1 << 15) | (k1 >> (32 - 15)); k1 *= c2; h1 ^= k1;
    };
    h1 ^= len; h1 ^= h1 >> 16; h1 *= 0x85ebca6b; h1 ^= h1 >> 13; h1 *= 0xc2b2ae35; h1 ^= h1 >> 16;
    return h1; 
}

class CountMinSketch {
private:
    int width;
    int depth;
    uint64_t total_count;
    std::vector<std::vector<uint32_t>> tables;
    std::vector<uint32_t> seeds;

    CountMinSketch(int w, int d, uint64_t tc) : width(w), depth(d), total_count(tc) {
        tables.resize(depth, std::vector<uint32_t>(width));
        std::mt19937 gen(1234);
        std::uniform_int_distribution<uint32_t> dist;
        for (int i = 0; i < depth; ++i) seeds.push_back(dist(gen));
    }

public:
    uint32_t query(const std::string& item) {
        uint32_t min_count = std::numeric_limits<uint32_t>::max();
        for (int i = 0; i < depth; ++i) {
            uint32_t index = murmur_hash(item.data(), item.length(), seeds[i]) % width;
            min_count = std::min(min_count, tables[i][index]);
        }
        return min_count;
    }
    static CountMinSketch load_from_file(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) { throw std::runtime_error("Error: No se pudo abrir el archivo .bin: " + filename); }
        int width, depth;
        uint64_t total_count;
        file.read(reinterpret_cast<char*>(&width), sizeof(width));
        file.read(reinterpret_cast<char*>(&depth), sizeof(depth));
        file.read(reinterpret_cast<char*>(&total_count), sizeof(total_count));
        CountMinSketch cms(width, depth, total_count);
        for (int i = 0; i < depth; ++i) {
            file.read(reinterpret_cast<char*>(cms.tables[i].data()), width * sizeof(uint32_t));
        }
        file.close();
        return cms;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Uso correcto: " << argv[0] << " <archivo_sketch.bin> <archivo_conteos.txt>" << std::endl;
        return 1;
    }

    std::string SKETCH_FILE = argv[1];
    std::string COUNTS_FILE = argv[2];

    std::cout << "Cargando sketch desde: " << SKETCH_FILE << std::endl;
    CountMinSketch cms = CountMinSketch::load_from_file(SKETCH_FILE);

    std::cout << "Cargando conteos exactos desde: " << COUNTS_FILE << std::endl;
    std::unordered_map<std::string, uint64_t> exact_counts;
    std::ifstream counts_in(COUNTS_FILE);
    if (!counts_in.is_open()) {
        std::cerr << "Error: No se pudo abrir el archivo de conteos." << std::endl;
        return 1;
    }
    std::string line;
    while (std::getline(counts_in, line)) {
        std::stringstream ss(line);
        std::string item;
        uint64_t count;
        ss >> item >> count;
        exact_counts[item] = count;
    }
    counts_in.close();
    
    uint64_t total_unique = exact_counts.size();
    std::cout << "Unicos cargados: " << total_unique << std::endl;
    std::cout << "Iniciando analisis global..." << std::endl;

    double total_relative_error = 0.0;
    uint64_t total_absolute_error = 0;
    uint64_t error_count = 0;

    for (const auto& pair : exact_counts) {
        const std::string& item = pair.first;
        uint64_t exact = pair.second;
        uint64_t estim = cms.query(item);
        uint64_t err = (estim > exact) ? (estim - exact) : 0;

        if (err > 0) {
            total_absolute_error += err;
            error_count++;
        }
        if (exact > 0) {
            total_relative_error += (double)err / exact;
        }
    }

    std::cout << "\n--- Resultados Globales ---" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Total Elementos: " << total_unique << std::endl;
    std::cout << "Elementos con Error: " << error_count << " (" << ((double)error_count/total_unique*100.0) << "%)" << std::endl;
    std::cout << "Piso de Ruido Promedio (Error Absoluto): " << ((double)total_absolute_error/total_unique) << std::endl;
    std::cout << "Error Relativo Global Promedio: " << ((total_relative_error/total_unique)*100.0) << "%" << std::endl;

    return 0;
}