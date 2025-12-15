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

uint32_t murmur_hash(const void * key, int len, uint32_t seed) {
    const uint8_t * data = (const uint8_t*)key;
    const int nblocks = len / 4;
    uint32_t h1 = seed;
    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;
    const uint32_t * blocks = (const uint32_t *)(data + nblocks*4);
    for(int i = -nblocks; i; i++) {
        uint32_t k1 = blocks[i];
        k1 *= c1; k1 = (k1 << 15) | (k1 >> (32 - 15)); k1 *= c2;
        h1 ^= k1; h1 = (h1 << 13) | (h1 >> (32 - 13)); h1 = h1*5+0xe6546b64;
    }
    const uint8_t * tail = (const uint8_t*)(data + nblocks*4);
    uint32_t k1 = 0;
    switch(len & 3) {
        case 3: k1 ^= tail[2] << 16;
        case 2: k1 ^= tail[1] << 8;
        case 1: k1 ^= tail[0];
                k1 *= c1; k1 = (k1 << 15) | (k1 >> (32 - 15)); k1 *= c2; h1 ^= k1;
    };
    h1 ^= len; h1 ^= h1 >> 16; h1 *= 0x85ebca6b; h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35; h1 ^= h1 >> 16;
    return h1;
}

class CountMinSketch {
private:
    int width;
    int depth;
    uint64_t total_count;
    std::vector<std::vector<uint32_t>> tables;
    std::vector<uint32_t> seeds;

    CountMinSketch(int w, int d, uint64_t tc)
        : width(w), depth(d), total_count(tc) {
        tables.resize(depth, std::vector<uint32_t>(width));
        
        std::mt19937 gen(1234); 
        std::uniform_int_distribution<uint32_t> dist;
        for (int i = 0; i < depth; ++i) {
            seeds.push_back(dist(gen));
        }
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

    uint64_t getTotalCount() const {
        return total_count;
    }

    static CountMinSketch load_from_file(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) {
            std::cerr << "Error: No se pudo abrir el archivo para cargar: " << filename << std::endl;
            throw std::runtime_error("No se pudo cargar el sketch.");
        }

        int width;
        int depth;
        uint64_t total_count;
        file.read(reinterpret_cast<char*>(&width), sizeof(width));
        file.read(reinterpret_cast<char*>(&depth), sizeof(depth));
        file.read(reinterpret_cast<char*>(&total_count), sizeof(total_count));

        CountMinSketch cms(width, depth, total_count);

        for (int i = 0; i < depth; ++i) {
            file.read(reinterpret_cast<char*>(cms.tables[i].data()), width * sizeof(uint32_t));
        }

        file.close();
        std::cout << "Sketch cargado exitosamente (" << total_count << " inserciones totales)." << std::endl;
        return cms;
    }
};

bool is_valid_kmer_char(char c) {
    return (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
}

std::string reverse_complement(const std::string& kmer) {
    static const std::map<char, char> complement = {
        {'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'}, {'N', 'N'}
    };
    std::string rev_comp;
    rev_comp.reserve(kmer.length());
    for (int i = kmer.length() - 1; i >= 0; --i) {
        if (complement.count(kmer[i])) {
            rev_comp += complement.at(kmer[i]);
        } else {
            rev_comp += 'N'; 
        }
    }
    return rev_comp;
}

std::string get_canonical_kmer(const std::string& kmer) {
    std::string rev_comp = reverse_complement(kmer);
    return std::min(kmer, rev_comp);
}

int main() {
    std::string SKETCH_FILE = "./sketch.bin";
    int K_MER_LENGTH = 21; 

    CountMinSketch cms = CountMinSketch::load_from_file(SKETCH_FILE);
    
    std::string input_kmer;
    std::cout << "\n--- Consultor de Sketch ---" << std::endl;
    std::cout << "Escribe un k-mer (k=" << K_MER_LENGTH << ") para consultar su frecuencia." << std::endl;
    std::cout << "Escribe 'exit' o 'quit' para salir." << std::endl;

    while (true) {
        std::cout << "\n> ";
        std::getline(std::cin, input_kmer);

        if (input_kmer == "exit" || input_kmer == "quit") {
            break;
        }

        std::string kmer_upper;
        bool is_valid = true;
        for (char c : input_kmer) {
            if (std::isspace(c)) continue; 
            char upper_c = std::toupper(static_cast<unsigned char>(c));
            if (!is_valid_kmer_char(upper_c)) {
                std::cout << "Error: El k-mer contiene caracteres invalidos ('" << upper_c << "')." << std::endl;
                is_valid = false;
                break;
            }
            kmer_upper += upper_c;
        }
        
        if (!is_valid) continue;

        if (kmer_upper.length() != K_MER_LENGTH) {
            std::cout << "Error: El k-mer debe tener longitud " << K_MER_LENGTH << "." << std::endl;
            continue;
        }

        std::string canonical = get_canonical_kmer(kmer_upper);
        uint32_t frequency = cms.query(canonical);

        std::cout << "  k-mer: " << kmer_upper << std::endl;
        std::cout << "  Canonico: " << canonical << std::endl;
        std::cout << "  Frecuencia Estimada (f_hat): " << frequency << std::endl;
    }

    std::cout << "Saliendo." << std::endl;
    return 0;
}