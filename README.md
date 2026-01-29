# kraken-kdb-inspect
Inpect k-mers sequences of specified taxid from Kraken/KrakenUniq's *.kdb file

Put kraken-kdb-inspect.cpp to kraken/src/, then compile it by "g++ -O3 -std=c++11 -o kraken-kdb-inspect kraken-kdb-inspect.cpp krakendb.cpp krakenutil.cpp quickfile.cpp seqreader.cpp -lz"

How to use:
./extract_kmers /path/to/database.kdb 1234 > taxid_1234_kmers.fasta
