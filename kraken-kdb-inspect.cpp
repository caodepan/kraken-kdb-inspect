#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "quickfile.hpp"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
using namespace kraken;

class KrakenDBExtractor {
public:
    static string decode_kmer(uint64_t kmer, int k) {
        string seq(k, ' ');
        for (int i = k - 1; i >= 0; i--) {
            int base = kmer & 3;
            seq[i] = "ACGT"[base];
            kmer >>= 2;
        }
        return seq;
    }

    static void extract_kmers_with_taxid(const string& db_filename, uint32_t target_taxid) {
        QuickFile db_file;
        db_file.open_file(db_filename);
        KrakenDB database(db_file.ptr());

        char *pair_ptr = database.get_pair_ptr();
        uint64_t key_ct = database.get_key_ct();
        uint64_t key_len = database.get_key_len();
        uint64_t val_len = database.get_val_len();
        uint64_t pair_sz = database.pair_size();
        uint8_t k = database.get_k();

        uint32_t counter = 1; // Counter for sequence numbering
        
        for (uint64_t i = 0; i < key_ct; i++) {
            char *current_pair = pair_ptr + i * pair_sz;
            
            // Extract k-mer (key)
            uint64_t kmer = 0;
            memcpy(&kmer, current_pair, key_len);
            kmer &= (1ULL << (k * 2)) - 1;  // trim to actual k-mer size
            
            // Extract taxid (value)
            uint32_t taxid;
            memcpy(&taxid, current_pair + key_len, val_len);
            
            if (taxid == target_taxid) {
                string kmer_seq = decode_kmer(kmer, k);
                // Output in FASTA format with header: >taxidX_Y
                cout << ">taxid" << target_taxid << "_" << counter << "\n" << kmer_seq << "\n";
                counter++;
            }
        }
    }
};

int main(int argc, char *argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <database.kdb> <taxid>" << endl;
        return 1;
    }
    
    string db_filename = argv[1];
    uint32_t target_taxid = atoi(argv[2]);
    
    KrakenDBExtractor::extract_kmers_with_taxid(db_filename, target_taxid);
    
    return 0;
}
