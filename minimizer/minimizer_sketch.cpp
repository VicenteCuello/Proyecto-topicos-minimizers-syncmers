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


uint32_t murmur_hash(const void * key, int len, uint32_t seed) {
    const uint8_t * data = (const uint8_t*)key;
    const int nblocks = len / 4;
    uint32_t h1 = seed;
    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;

    const uint32_t * blocks = (const uint32_t *)(data + nblocks*4);
    for(int i = -nblocks; i; i++) {
        uint32_t k1 = blocks[i];
        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> (32 - 15));
        k1 *= c2;
        h1 ^= k1;
        h1 = (h1 << 13) | (h1 >> (32 - 13)); 
        h1 = h1*5+0xe6546b64;
    }

    const uint8_t * tail = (const uint8_t*)(data + nblocks*4);
    uint32_t k1 = 0;
    switch(len & 3) {
        case 3: k1 ^= tail[2] << 16;
        case 2: k1 ^= tail[1] << 8;
        case 1: k1 ^= tail[0];
                k1 *= c1; k1 = (k1 << 15) | (k1 >> (32 - 15)); k1 *= c2; h1 ^= k1;
    };

    h1 ^= len;
    h1 ^= h1 >> 16;
    h1 *= 0x85ebca6b;
    h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35;
    h1 ^= h1 >> 16;
    return h1;
}

class CountMinSketch {
private:
    int width;
    int depth;
    std::vector<std::vector<uint32_t>> tables;
    std::vector<uint32_t> seeds;
    uint64_t total_count;

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

    bool save_to_file(const std::string& filename) const {
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            std::cerr << "Error: No se pudo abrir el archivo para guardar: " << filename << std::endl;
            return false;
        }

        file.write(reinterpret_cast<const char*>(&width), sizeof(width));
        file.write(reinterpret_cast<const char*>(&depth), sizeof(depth));
        file.write(reinterpret_cast<const char*>(&total_count), sizeof(total_count));

        for (int i = 0; i < depth; ++i) {
            file.write(reinterpret_cast<const char*>(tables[i].data()), width * sizeof(uint32_t));
        }

        file.close();
        return true;
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
        return cms;
    }
};

long get_peak_memory_rss() {
    std::ifstream status_file("/proc/self/status");
    if (!status_file.is_open()) {
        return -1;
    }
    
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
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
            return false;
        }
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


void process_sequence(const std::string& sequence, int w, int k, CountMinSketch& cms) {
    if (sequence.length() < k + w - 1) {
        return;
    }

    std::deque<std::pair<uint32_t, int>> mq;
    std::deque<std::string> kmer_buffer;
    int last_minimizer_pos = -1;
    uint32_t hash_seed = 0;

    for (int i = 0; i <= (int)sequence.length() - k; ++i) {
        
        std::string kmer = sequence.substr(i, k);

        if (!is_valid_kmer(kmer)) {

            mq.clear();
            kmer_buffer.clear();
            continue; 
        }


        std::string canonical = get_canonical_kmer(kmer);
        uint32_t hash_val = murmur_hash(canonical.data(), canonical.length(), hash_seed);
        
        kmer_buffer.push_back(canonical);

        while (!mq.empty() && mq.back().first > hash_val) {
            mq.pop_back();
        }
        mq.push_back({hash_val, i});

        if (i >= w - 1) {
            int window_start_pos = i - w + 1;

            while (!mq.empty() && mq.front().second < window_start_pos) {
                mq.pop_front();
            }

            int minimizer_pos = mq.front().second;

            if (minimizer_pos != last_minimizer_pos) {
                int buffer_idx = minimizer_pos - window_start_pos;

                if (buffer_idx < kmer_buffer.size()) {
                    cms.update(kmer_buffer[buffer_idx]);
                    last_minimizer_pos = minimizer_pos;
                } else {
                    
                }
            }
            
           
            if (!kmer_buffer.empty()) {
                kmer_buffer.pop_front();
            }
        }
    }
}


void process_fasta_file(const std::filesystem::path& filepath, int w, int k, CountMinSketch& cms) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error abriendo archivo: " << filepath << std::endl;
        return;
    }

    std::string line;
    std::stringstream seq_buffer;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            std::string sequence = seq_buffer.str();
            if (!sequence.empty()) {
                process_sequence(sequence, w, k, cms);
            }
            seq_buffer.str(""); 
        } else {
           
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
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
        process_sequence(sequence, w, k, cms);
    }
}


int main() {
    std::string INPUT_DIR = "../Genomas_Helicobacter"; 
    std::string SKETCH_FILE = "./sketch.bin";
    std::string METRICS_FILE = "./processing_metrics.txt";
    int W_WINDOW = 10;
    int K_MER = 21;
    
    int CMS_WIDTH = 2000000;
    int CMS_DEPTH = 5;

    std::cout << "Iniciando Pipeline C++ (v2 Corregido)..." << std::endl;
    std::cout << "Parametros: w=" << W_WINDOW << ", k=" << K_MER << std::endl;
    std::cout << "Sketch: width=" << CMS_WIDTH << ", depth=" << CMS_DEPTH << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    CountMinSketch cms(CMS_WIDTH, CMS_DEPTH);
    
    int processed_files = 0;
    int total_files = 0;
    
    if (!std::filesystem::exists(INPUT_DIR)) {
        std::cerr << "Error: El directorio no existe: " << INPUT_DIR << std::endl;
        return 1;
    }
    for (const auto& entry : std::filesystem::directory_iterator(INPUT_DIR)) {
        std::string ext = entry.path().extension().string();
        if (ext == ".fna" || ext == ".fasta") {
            total_files++;
        }
    }

    if (total_files == 0) {
        std::cerr << "Error: No se encontraron archivos .fna o .fasta en " << INPUT_DIR << std::endl;
        return 1;
    }

    for (const auto& entry : std::filesystem::directory_iterator(INPUT_DIR)) {
        std::string ext = entry.path().extension().string();
        if (ext == ".fna" || ext == ".fasta") {
            processed_files++;
            std::cout << "[" << processed_files << "/" << total_files << "] Procesando: " 
                      << entry.path().filename() << "...\r" << std::flush;
            
            process_fasta_file(entry.path(), W_WINDOW, K_MER, cms);
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    double total_time_seconds = duration.count();

    long peak_mem_kb = get_peak_memory_rss();
    double peak_mem_mb = (peak_mem_kb > 0) ? (peak_mem_kb / 1024.0) : 0.0;
    uint64_t total_insertions = cms.getTotalCount();

    std::cout << "\n\n--- Procesamiento completado ---" << std::endl;
    std::cout << "Total de inserciones en el Sketch: " << cms.getTotalCount() << std::endl;
    std::cout << "Tiempo total de procesamiento: " << total_time_seconds << " segundos." << std::endl;
    std::cout << "Peak de memoria (RSS HWM): " << peak_mem_mb << " MB." << std::endl;

    std::cout << "Guardando sketch en " << SKETCH_FILE << "..." << std::endl;
    if (cms.save_to_file(SKETCH_FILE)) {
        std::cout << "Sketch guardado exitosamente." << std::endl;
    } else {
        std::cout << "Error al guardar el sketch." << std::endl;
    }

    std::cout << "Guardando métricas en " << METRICS_FILE << "..." << std::endl;
    std::ofstream metrics_file(METRICS_FILE);
    if (metrics_file.is_open()) {
        metrics_file << "--- Métricas de Procesamiento ---" << std::endl;
        metrics_file << "Total Archivos FNA: " << total_files << std::endl;
        metrics_file << "Parametros (w, k): (" << W_WINDOW << ", " << K_MER << ")" << std::endl;
        metrics_file << "Sketch (width, depth): (" << CMS_WIDTH << ", " << CMS_DEPTH << ")" << std::endl;
        metrics_file << "Total Inserciones (N): " << total_insertions << std::endl;
        metrics_file << "Tiempo (segundos): " << total_time_seconds << std::endl;
        metrics_file << "Peak Memoria (MB): " << peak_mem_mb << std::endl;
        metrics_file.close();
        std::cout << "Métricas guardadas." << std::endl;
    } else {
        std::cout << "Error al guardar las métricas." << std::endl;
    }

    std::string test_query = "AAAAATTTTTCCCCCGGGGGAT"; 
    std::cout << "Frecuencia estimada de un k-mer de prueba: " << cms.query(test_query) << std::endl;

    return 0;
}