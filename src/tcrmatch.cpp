#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <string>
#include <tuple>
#include <unistd.h>
#include <vector>

std::array<std::array<float, 20>, 20> k1;
int p_kmin = 1;
int p_kmax = 30;
float p_beta = 0.11387;
// Hardcoded because parsing + computing matrix is annoying
float blm_qij[20][20] = {
    {0.0215, 0.0023, 0.0019, 0.0022, 0.0016, 0.0019, 0.003,
     0.0058, 0.0011, 0.0032, 0.0044, 0.0033, 0.0013, 0.0016,
     0.0022, 0.0063, 0.0037, 0.0004, 0.0013, 0.0051},
    {0.0023, 0.0178, 0.002,  0.0016, 0.0004, 0.0025, 0.0027,
     0.0017, 0.0012, 0.0012, 0.0024, 0.0062, 0.0008, 0.0009,
     0.001,  0.0023, 0.0018, 0.0003, 0.0009, 0.0016},
    {0.0019, 0.002,  0.0141, 0.0037, 0.0004, 0.0015, 0.0022,
     0.0029, 0.0014, 0.001,  0.0014, 0.0024, 0.0005, 0.0008,
     0.0009, 0.0031, 0.0022, 0.0002, 0.0007, 0.0012},
    {0.0022, 0.0016, 0.0037, 0.0213, 0.0004, 0.0016, 0.0049,
     0.0025, 0.001,  0.0012, 0.0015, 0.0024, 0.0005, 0.0008,
     0.0012, 0.0028, 0.0019, 0.0002, 0.0006, 0.0013},
    {0.0016, 0.0004, 0.0004, 0.0004, 0.0119, 0.0003, 0.0004,
     0.0008, 0.0002, 0.0011, 0.0016, 0.0005, 0.0004, 0.0005,
     0.0004, 0.001,  0.0009, 0.0001, 0.0003, 0.0014},
    {0.0019, 0.0025, 0.0015, 0.0016, 0.0003, 0.0073, 0.0035,
     0.0014, 0.001,  0.0009, 0.0016, 0.0031, 0.0007, 0.0005,
     0.0008, 0.0019, 0.0014, 0.0002, 0.0007, 0.0012},
    {0.003,  0.0027, 0.0022, 0.0049, 0.0004, 0.0035, 0.0161,
     0.0019, 0.0014, 0.0012, 0.002,  0.0041, 0.0007, 0.0009,
     0.0014, 0.003,  0.002,  0.0003, 0.0009, 0.0017},
    {0.0058, 0.0017, 0.0029, 0.0025, 0.0008, 0.0014, 0.0019,
     0.0378, 0.001,  0.0014, 0.0021, 0.0025, 0.0007, 0.0012,
     0.0014, 0.0038, 0.0022, 0.0004, 0.0008, 0.0018},
    {0.0011, 0.0012, 0.0014, 0.001,  0.0002, 0.001,  0.0014,
     0.001,  0.0093, 0.0006, 0.001,  0.0012, 0.0004, 0.0008,
     0.0005, 0.0011, 0.0007, 0.0002, 0.0015, 0.0006},
    {0.0032, 0.0012, 0.001,  0.0012, 0.0011, 0.0009, 0.0012,
     0.0014, 0.0006, 0.0184, 0.0114, 0.0016, 0.0025, 0.003,
     0.001,  0.0017, 0.0027, 0.0004, 0.0014, 0.012},
    {0.0044, 0.0024, 0.0014, 0.0015, 0.0016, 0.0016, 0.002,
     0.0021, 0.001,  0.0114, 0.0371, 0.0025, 0.0049, 0.0054,
     0.0014, 0.0024, 0.0033, 0.0007, 0.0022, 0.0095},
    {0.0033, 0.0062, 0.0024, 0.0024, 0.0005, 0.0031, 0.0041,
     0.0025, 0.0012, 0.0016, 0.0025, 0.0161, 0.0009, 0.0009,
     0.0016, 0.0031, 0.0023, 0.0003, 0.001,  0.0019},
    {0.0013, 0.0008, 0.0005, 0.0005, 0.0004, 0.0007, 0.0007,
     0.0007, 0.0004, 0.0025, 0.0049, 0.0009, 0.004,  0.0012,
     0.0004, 0.0009, 0.001,  0.0002, 0.0006, 0.0023},
    {0.0016, 0.0009, 0.0008, 0.0008, 0.0005, 0.0005, 0.0009,
     0.0012, 0.0008, 0.003,  0.0054, 0.0009, 0.0012, 0.0183,
     0.0005, 0.0012, 0.0012, 0.0008, 0.0042, 0.0026},
    {0.0022, 0.001,  0.0009, 0.0012, 0.0004, 0.0008, 0.0014,
     0.0014, 0.0005, 0.001,  0.0014, 0.0016, 0.0004, 0.0005,
     0.0191, 0.0017, 0.0014, 0.0001, 0.0005, 0.0012},
    {0.0063, 0.0023, 0.0031, 0.0028, 0.001,  0.0019, 0.003,
     0.0038, 0.0011, 0.0017, 0.0024, 0.0031, 0.0009, 0.0012,
     0.0017, 0.0126, 0.0047, 0.0003, 0.001,  0.0024},
    {0.0037, 0.0018, 0.0022, 0.0019, 0.0009, 0.0014, 0.002,
     0.0022, 0.0007, 0.0027, 0.0033, 0.0023, 0.001,  0.0012,
     0.0014, 0.0047, 0.0125, 0.0003, 0.0009, 0.0036},
    {0.0004, 0.0003, 0.0002, 0.0002, 0.0001, 0.0002, 0.0003,
     0.0004, 0.0002, 0.0004, 0.0007, 0.0003, 0.0002, 0.0008,
     0.0001, 0.0003, 0.0003, 0.0065, 0.0009, 0.0004},
    {0.0013, 0.0009, 0.0007, 0.0006, 0.0003, 0.0007, 0.0009,
     0.0008, 0.0015, 0.0014, 0.0022, 0.001,  0.0006, 0.0042,
     0.0005, 0.001,  0.0009, 0.0009, 0.0102, 0.0015},
    {0.0051, 0.0016, 0.0012, 0.0013, 0.0014, 0.0012, 0.0017,
     0.0018, 0.0006, 0.012,  0.0095, 0.0019, 0.0023, 0.0026,
     0.0012, 0.0024, 0.0036, 0.0004, 0.0015, 0.0196}};

struct peptide {
  std::string seq;
  int len;
  float aff;
  std::vector<int> i;
};

std::vector<std::string> read_IEDB_data(std::string IEDB_data_file) {
  std::vector<std::string> iedb_data;
  std::ifstream iedb_file(IEDB_data_file);
  std::string line;
  while (getline(iedb_file, line)) {
    std::stringstream ss(line);
    std::string sequence;
    ss >> sequence;
    if (sequence != "trimmed_seq") {
      iedb_data.push_back(sequence);
    }
  }
  return iedb_data;
}

std::array<std::array<float, 20>, 20> fmatrix_k1() {
  // Calculates the modified (normalized) blosum62 matrix

  int k, j;
  float marg[20];
  float sum;

  // initialize margin array
  for (int i = 0; i < 20; i++) {
    marg[i] = 0.0;
  }
  // initialize k1
  for (int i = 0; i < 20; i++) {
    for (int j = 0; j < 20; j++) {
      k1[i][j] = 0.0;
    }
  }
  // normalize matrix by marginal frequencies
  for (j = 0; j < 20; j++) {
    sum = 0;
    for (k = 0; k < 20; k++)
      sum += blm_qij[j][k];
    marg[j] = sum;
  }
  // calculate K1
  for (j = 0; j < 20; j++) {
    for (k = 0; k < 20; k++) {
      k1[j][k] = blm_qij[j][k] / (marg[j] * marg[k]);
      k1[j][k] = pow(k1[j][k], p_beta);
    }
  }

  return (k1);
}

float k3_sum(peptide pep1, peptide pep2) {
  // Recursively calculate Kernel 3 using Kernel 1 lookups
  float k2, term, k3 = 0.0;
  int start1, start2;
  int k, j1, j2;

  float k2_prod_save[31][31][31];

  for (k = p_kmin; k <= p_kmax; k++) {
    for (start1 = 0; start1 <= pep1.len - k; start1++) {
      for (start2 = 0; start2 <= pep2.len - k; start2++) {

        j1 = pep1.i[start1 + k - 1];
        j2 = pep2.i[start2 + k - 1];
        term = k1[j1][j2];

        if (k == 1) {
          k2 = term;
        } else {
          k2 = k2_prod_save[k - 1][start1][start2] * term;
        }

        k2_prod_save[k][start1][start2] = k2;
        k3 += k2;
      }
    }
  }
  return (k3);
}

void multi_calc_k3(std::vector<peptide> peplist1, std::vector<peptide> peplist2,
                   float threshold) {
  // Simple method to calculate pairwise TCRMatch scores using two peptide
  // vectors
  std::vector<std::tuple<std::string, std::string, float>>
      results[omp_get_max_threads()];
#pragma omp parallel for
  for (int i = 0; i < peplist1.size(); i++) {
    for (int j = 0; j < peplist2.size(); j++) {
      peptide pep1 = peplist1[i];
      peptide pep2 = peplist2[j];
      float score = 0.0;
      score = k3_sum(pep1, pep2) / sqrt(pep1.aff * pep2.aff);
      if (score > threshold) {
        int tid = omp_get_thread_num();
        results[tid].push_back(make_tuple(pep1.seq, pep2.seq, score));
      }
    }
  }
  for (int i = 0; i < omp_get_max_threads(); i++) {
    for (auto &tuple : results[i]) {
      std::cout << std::fixed << std::setprecision(2) << std::get<0>(tuple)
                << " " << std::get<1>(tuple) << " " << std::get<2>(tuple)
                << std::endl;
    }
  }
}

// Move this to outside -> import everything you need
int main(int argc, char *argv[]) {
  int opt;
  int n_threads;
  float threshold;
  std::string in_file;
  std::string iedb_file = "data/IEDB_data.tsv";
  int i_flag = -1;
  int t_flag = -1;
  int thresh_flag = -1;

  // Command line argument parsing
  while ((opt = getopt(argc, argv, "t:i:s:d:")) != -1) {
    switch (opt) {
    case 't':
      n_threads = std::stoi(optarg);
      t_flag = 1;
      break;
    case 'i':
      in_file = optarg;
      i_flag = 1;
      break;
    case 's':
      threshold = std::stof(optarg);
      thresh_flag = 1;
      break;
    case 'd':
      iedb_file = optarg;
      break;
    default:
      std::cerr << "Usage: ./tcrmatch -i infile_name.txt -t num_threads -s "
                   "score_threshold -d /path/to/database"
                << std::endl;
      return EXIT_FAILURE;
    }
  }
  
  // Check that required parameters are there + error correcting
  if (i_flag == -1 || t_flag == -1) {
    std::cerr << "Missing mandatory parameters" << std::endl
              << "Usage: ./tcrmatch -i infile_name.txt -t num_threads"
              << std::endl;
    return EXIT_FAILURE;
  }
  if (n_threads < 1) {
    std::cerr << "Number of threads cannot be less than one." << std::endl;
    return EXIT_FAILURE;
  }
  if (thresh_flag == -1) {
    threshold = .97;
  }
  if (t_flag == -1){
    n_threads = 1;
  }
  if (threshold < 0 || threshold > 1) {
    std::cerr << "Threshold must be between 0 and 1" << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<std::string> iedb_data = read_IEDB_data(iedb_file);
  std::ifstream file1(in_file);
  std::string line;
  std::string alphabet;
  std::vector<peptide> peplist1;
  std::vector<peptide> peplist2;

  omp_set_num_threads(n_threads);

  alphabet = "ARNDCQEGHILKMFPSTWYV";
  k1 = fmatrix_k1();
  while (getline(file1, line)) {
    std::vector<int> int_vec;
    for (int i = 0; i < line.length(); i++) {
      if (alphabet.find(line[i]) == -1) {
        std::cerr << "Invalid amino acid found in " << line << " at position "
                  << i + 1 << std::endl;
        return EXIT_FAILURE;
      }
    }
    peplist1.push_back({line, int(line.length()), -99.9, int_vec});
  }
  file1.close();

// Calculate the normalization score (aff) (kernel 3 self vs self) list 1
#pragma omp parallel for
  for (int i = 0; i < peplist1.size(); i++) {
    peptide *pep_ptr = &peplist1[i];
    for (int x = 0; x < pep_ptr->len; x++) {
      pep_ptr->i.push_back(alphabet.find(pep_ptr->seq[x]));
    }
    pep_ptr->aff = k3_sum(*pep_ptr, *pep_ptr);
  }

  // change to IEDB data
  for (std::vector<std::string>::iterator it = iedb_data.begin();
       it != iedb_data.end(); it++) {
    std::vector<int> int_vec;
    for (int i = 0; i < (*it).length(); i++) {
      if (alphabet.find((*it)[i]) == -1) {
        std::cerr << "Invalid amino acid found in " << *it << " at position "
                  << i + 1 << std::endl;
        return EXIT_FAILURE;
      }
    }
    peplist2.push_back({*it, int((*it).length()), -99.9, int_vec});
  }

// Calculate the normalization score (aff) (kernel 3 self vs self) for list 2
#pragma omp parallel for
  for (int i = 0; i < peplist2.size(); i++) {
    peptide *pep_ptr = &peplist2[i];
    for (int x = 0; x < pep_ptr->len; x++) {
      pep_ptr->i.push_back(alphabet.find(pep_ptr->seq[x]));
    }
    pep_ptr->aff = k3_sum(*pep_ptr, *pep_ptr);
  }
  multi_calc_k3(peplist1, peplist2, threshold);

  return 0;
}
