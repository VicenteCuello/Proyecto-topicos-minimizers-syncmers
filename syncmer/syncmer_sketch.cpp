#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <deque>
#include <filesystem>
#include <cstdint>
#include <limits>
#include <random>
#include <chrono>

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
public:
    CountMinSketch(int w, int d) : width(w), depth(d), total_count(0) {
        tables.resize(depth, std::vector<uint32_t>(width, 0));
        std::mt19937 gen(1234);
        std::uniform_int_distribution<uint32_t> dist;
        for (int i = 0; i < depth; ++i) {
            seeds.push_back(dist(gen));
        }
    }
    void update(const std::string& item, uint32_t count = 1) {
        total_count += count;
        for (int i = 0; i < depth; ++i) {
            uint32_t index = murmur_hash(item.data(), item.length(), seeds[i]) % width;
            tables[i][index] += count;
        }
    }
    uint64_t getTotalCount() const { return total_count; }
    
    bool save_to_file(const std::string& filename) const {
        std::ofstream file(filename, std::ios::binary);
        if (!file) return false;
        file.write(reinterpret_cast<const char*>(&width), sizeof(width));
        file.write(reinterpret_cast<const char*>(&depth), sizeof(depth));
        file.write(reinterpret_cast<const char*>(&total_count), sizeof(total_count));
        for (int i = 0; i < depth; ++i) {
            file.write(reinterpret_cast<const char*>(tables[i].data()), width * sizeof(uint32_t));
        }
        file.close();
        return true;
    }
};

long get_peak_memory_rss() {
    std::ifstream status_file("/proc/self/status");
    if (!status_file.is_open()) return -1;
    std::string line;
    while (std::getline(status_file, line)) {
        if (line.rfind("VmHWM:", 0) == 0) {
            std::stringstream ss(line);
            std::string key;
            long value;
            ss >> key >> value; 
            return value;
        }
    }
    return -1; 
}

bool is_valid_kmer(const std::string& kmer) {
    for (char c : kmer) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') return false;
    }
    return true;
}
std::string reverse_complement(const std::string& kmer) {
    static const std::map<char, char> complement = {
        {'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'}, {'N', 'N'}
    };
    std::string rev_comp;
    rev_comp.reserve(kmer.length());
    for (int i = kmer.length() - 1; i >= 0; --i) {
        rev_comp += complement.at(kmer[i]);
    }
    return rev_comp;
}
std::string get_canonical_kmer(const std::string& kmer) {
    std::string rev_comp = reverse_complement(kmer);
    return std::min(kmer, rev_comp);
}


void process_sequence_syncmer(const std::string& sequence, int k, int s, int t, CountMinSketch& cms) {
    if (sequence.length() < k) {
        return;
    }

    std::deque<std::pair<uint32_t, int>> mq;
    
    uint32_t hash_seed = 0;

    int syncmer_window_size = k - s + 1;

    for (int i = 0; i <= (int)sequence.length() - s; ++i) {
        
        std::string s_mer = sequence.substr(i, s);

        if (!is_valid_kmer(s_mer)) {
            mq.clear();
            continue;
        }

        std::string canonical_s = get_canonical_kmer(s_mer);
        uint32_t hash_val = murmur_hash(canonical_s.data(), canonical_s.length(), hash_seed);

        while (!mq.empty() && mq.back().first > hash_val) {
            mq.pop_back();
        }
        mq.push_back({hash_val, i});

        if (i >= syncmer_window_size - 1) {
            int window_start_pos = i - syncmer_window_size + 1; 

            while (!mq.empty() && mq.front().second < window_start_pos) {
                mq.pop_front();
            }

            int min_s_mer_pos = mq.front().second;
            int relative_pos = min_s_mer_pos - window_start_pos;

            if (relative_pos == t) {
                std::string kmer = sequence.substr(window_start_pos, k);
                std::string canonical_k = get_canonical_kmer(kmer);
                
                cms.update(canonical_k);
            }
        }
    }
}

void process_fasta_file(const std::filesystem::path& filepath, int k, int s, int t, CountMinSketch& cms) {
    std::ifstream file(filepath);
    if (!file.is_open()) {  return; }

    std::string line;
    std::stringstream seq_buffer;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            std::string sequence = seq_buffer.str();
            if (!sequence.empty()) {
                process_sequence_syncmer(sequence, k, s, t, cms);
            }
            seq_buffer.str(""); 
        } else {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            for(char c : line) {
                char upper_c = std::toupper(static_cast<unsigned char>(c));
                if (upper_c == 'A' || upper_c == 'C' || upper_c == 'G' || upper_c == 'T' || upper_c == 'N') {
                    seq_buffer << upper_c;
                }
            }
        }
    }
    std::string sequence = seq_buffer.str();
    if (!sequence.empty()) {
        process_sequence_syncmer(sequence, k, s, t, cms);
    }
}

int main() {
    std::string INPUT_DIR = "../Genomas_Helicobacter"; 
    std::string SKETCH_FILE = "./syncmer_sketch.bin";
    std::string METRICS_FILE = "./syncmer_metrics.txt";
    
    int K_MER = 21;  
    int S_MER = 11; 
    int T_OFFSET = 5; 

    int CMS_WIDTH = 2000000;
    int CMS_DEPTH = 5;

    std::cout << "Iniciando Pipeline C++ (Open Syncmers)..." << std::endl;
    std::cout << "Parametros: k=" << K_MER << ", s=" << S_MER << ", t=" << T_OFFSET << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();

    CountMinSketch cms(CMS_WIDTH, CMS_DEPTH);
    
    int processed_files = 0;
    int total_files = 0;
    
    if (!std::filesystem::exists(INPUT_DIR)) { return 1; }
    for (const auto& entry : std::filesystem::directory_iterator(INPUT_DIR)) {
        std::string ext = entry.path().extension().string();
        if (ext == ".fna" || ext == ".fasta") { total_files++; }
    }

    for (const auto& entry : std::filesystem::directory_iterator(INPUT_DIR)) {
        std::string ext = entry.path().extension().string();
        if (ext == ".fna" || ext == ".fasta") {
            processed_files++;
            std::cout << "[" << processed_files << "/" << total_files << "] Procesando: " 
                      << entry.path().filename() << "...\r" << std::flush;
            
            process_fasta_file(entry.path(), K_MER, S_MER, T_OFFSET, cms);
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    long peak_mem_kb = get_peak_memory_rss();
    double peak_mem_mb = (peak_mem_kb > 0) ? (peak_mem_kb / 1024.0) : 0.0;
    uint64_t total_insertions = cms.getTotalCount();

    std::cout << "\n\n--- Procesamiento de Syncmers Completado ---" << std::endl;
    std::cout << "Total de inserciones en el Sketch: " << total_insertions << std::endl;
    std::cout << "Tiempo total: " << duration.count() << " segundos." << std::endl;
    std::cout << "Peak de memoria: " << peak_mem_mb << " MB." << std::endl;

    cms.save_to_file(SKETCH_FILE);
    
    std::ofstream metrics_file(METRICS_FILE);
    if (metrics_file.is_open()) {
        metrics_file << "--- Métricas Syncmers ---" << std::endl;
        metrics_file << "Total Inserciones: " << total_insertions << std::endl;
        metrics_file << "Tiempo: " << duration.count() << std::endl;
        metrics_file << "Peak Memoria: " << peak_mem_mb << std::endl;
        metrics_file.close();
    }

    return 0;
}