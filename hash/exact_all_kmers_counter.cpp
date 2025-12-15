#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <filesystem>
#include <cstdint>
#include <limits>
#include <unordered_map> 
#include <chrono>

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


void process_sequence(const std::string& sequence, int k, 
                      std::unordered_map<std::string, uint64_t>& counts) { 
    if (sequence.length() < k) {
        return;
    }

    for (int i = 0; i <= (int)sequence.length() - k; ++i) {
        std::string kmer = sequence.substr(i, k);

        if (!is_valid_kmer(kmer)) {
            continue; 
        }

        std::string canonical = get_canonical_kmer(kmer);
        

        counts[canonical]++;
    }
}

void process_fasta_file(const std::filesystem::path& filepath, int k, 
                        std::unordered_map<std::string, uint64_t>& counts) { 
    std::ifstream file(filepath);
    if (!file.is_open()) {  }

    std::string line;
    std::stringstream seq_buffer;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            std::string sequence = seq_buffer.str();
            if (!sequence.empty()) {
                process_sequence(sequence, k, counts); 
            }
            seq_buffer.str(""); 
        } else {
            if (!line.empty() && line.back() == '\r') { line.pop_back(); }
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
        process_sequence(sequence, k, counts); 
    }
}

int main() {
    std::string INPUT_DIR = "../Genomas_Helicobacter";
    std::string OUTPUT_COUNTS_FILE = "./exact_all_kmers_countsT.txt";
    std::string METRICS_FILE = "./exact_all_kmers_metrics.txtT";
    int K_MER = 21;
    
    std::cout << "Iniciando Pipeline C++ (Contador Exacto de TODOS los K-mers)..." << std::endl;
    std::cout << "Parametros: k=" << K_MER << std::endl;
    std::cout << "Usando: std::unordered_map (¡ALTO USO DE RAM!)" << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    std::unordered_map<std::string, uint64_t> exact_counts;
    
    int processed_files = 0;
    int total_files = 0;
    
    if (!std::filesystem::exists(INPUT_DIR)) { }
    for (const auto& entry : std::filesystem::directory_iterator(INPUT_DIR)) {
        std::string ext = entry.path().extension().string();
        if (ext == ".fna" || ext == ".fasta") { total_files++; }
    }
    if (total_files == 0) {  }

    for (const auto& entry : std::filesystem::directory_iterator(INPUT_DIR)) {
        std::string ext = entry.path().extension().string();
        if (ext == ".fna" || ext == ".fasta") {
            processed_files++;
            std::cout << "[" << processed_files << "/" << total_files << "] Procesando: " 
                      << entry.path().filename() << "...\r" << std::flush;
            
            process_fasta_file(entry.path(), K_MER, exact_counts); 
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    
    long peak_mem_kb = get_peak_memory_rss();
    double peak_mem_mb = (peak_mem_kb > 0) ? (peak_mem_kb / 1024.0) : 0.0;
    uint64_t total_unique_kmers = exact_counts.size();

    std::cout << "\n\n--- Procesamiento de TODOS los K-mers Completado ---" << std::endl;
    std::cout << "Total de k-mers UNICOS encontrados: " << total_unique_kmers << std::endl;
    std::cout << "Tiempo total de procesamiento: " << duration.count() << " segundos." << std::endl;
    std::cout << "Peak de memoria (RSS HWM): " << peak_mem_mb << " MB." << std::endl;

    std::cout << "Guardando conteos exactos en " << OUTPUT_COUNTS_FILE << "..." << std::endl;
    std::ofstream counts_file(OUTPUT_COUNTS_FILE);
    if (counts_file.is_open()) {
        for (const auto& pair : exact_counts) {
            counts_file << pair.first << "\t" << pair.second << "\n";
        }
        counts_file.close();
    } else {
        std::cerr << "Error al guardar los conteos." << std::endl;
    }

    std::cout << "Guardando métricas en " << METRICS_FILE << "..." << std::endl;
    std::ofstream metrics_file(METRICS_FILE);
    if (metrics_file.is_open()) {
        metrics_file << "--- Métricas de Procesamiento (Exacto, Todos los K-mers) ---" << std::endl;
        metrics_file << "Total Archivos FNA: " << total_files << std::endl;
        metrics_file << "Parametros (k): (" << K_MER << ")" << std::endl;
        metrics_file << "Total K-mers Unicos: " << total_unique_kmers << std::endl;
        metrics_file << "Tiempo (segundos): " << duration.count() << std::endl;
        metrics_file << "Peak Memoria (MB): " << peak_mem_mb << std::endl;
        metrics_file.close();
    } else {
        std::cerr << "Error al guardar las métricas." << std::endl;
    }

    return 0;
}