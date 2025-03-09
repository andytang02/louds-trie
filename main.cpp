#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <vector>
#include <sys/resource.h>

#include "louds-trie.hpp"

using namespace std;
using namespace std::chrono;

long getMemoryUsage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss; // Memory usage in kilobytes
}

void test1() {
    louds::Trie trie1;
    vector<string> keys1 = {"b", "ba", "bs", "bat", "bar", "baj", "bsj", "bst", "c", "ca", "cs", "cat", "car", "caj", "csj", "cst"};
    sort(keys1.begin(), keys1.end());
  
    for(auto & key: keys1) {
      trie1.add(key);
    }
  
    trie1.build();
  
    trie1.print();
  
    louds::Trie trie2;
    vector<string> keys2 = {"c", "ca", "cat", "ba", "car", "d", "da", "dar", "dart", "cart", "cats"};
    sort(keys2.begin(), keys2.end());
  
    for(auto & key: keys2) {
      trie2.add(key);
    }
  
    trie2.build();
    
    /*
    louds::Trie* trie3 = louds::merge_naive2(trie1, trie2);
  
    trie3->build();
  
    trie3->print();
  
    vector<string> extracted_keys = trie3->extract_keys();
  
    cout << "MERGED TRIE" << endl;
    for (string &s: extracted_keys) {
      cout << s << endl;
    }
          */
  
    louds::Trie* trie4 = louds::Trie::merge_naive(trie1, trie2);
  
    trie4->build();
  
    trie4->print();
  
    vector<string> extracted_keys = trie4->extract_keys();
  
    cout << "MERGED TRIE" << endl;
    for (string &s: extracted_keys) {
      cout << s << endl;
    }


    louds::Trie* trie5 = louds::Trie::merge_efficient(trie1, trie2);
  
    trie5->build();
  
    trie5->print();
    
    extracted_keys = trie5->extract_keys();
  
    cout << "MERGED TRIE" << endl;
    for (string &s: extracted_keys) {
      cout << s << endl;
    }
}

void test2() {
    ifstream file1("data/enwiki-latest-all-titles-in-ns0.sorted");

    string line;
    vector<string> keys1;
    while (getline(file1, line)) {
        keys1.push_back(line);
    }
    file1.close();

    ifstream file2("data/jawiki-latest-all-titles-in-ns0.sorted");
    vector<string> keys2;
    while (getline(file2, line)) {
        keys2.push_back(line);
    }
    file2.close();

    louds::Trie trie1;
    for (uint64_t i = 0; i < keys1.size(); ++i) {
      trie1.add(keys1[i]);
    }
    trie1.build();

    louds::Trie trie2;
    for (uint64_t i = 0; i < keys2.size(); ++i) {
      trie2.add(keys2[i]);
    }
    trie2.build();

    high_resolution_clock::time_point begin = high_resolution_clock::now();
    long before = getMemoryUsage();

    louds::Trie* trie3 = louds::Trie::merge_naive(trie1, trie2);

    long after = getMemoryUsage();
    std::cout << "Memory used: " << (after - before) * 1024 << " bytes" << std::endl;

    high_resolution_clock::time_point end = high_resolution_clock::now();
    double elapsed = (double)duration_cast<milliseconds>(end - begin).count();
    cout << "Naive Implementation elapsed time: " << elapsed << endl;

    trie3->build();

    for (uint64_t i = 0; i < keys1.size(); ++i) {
        assert(trie3->lookup(keys1[i]) != -1);
    }

    for (uint64_t i = 0; i < keys2.size(); ++i) {
        assert(trie3->lookup(keys2[i]) != -1);
    }
}

void test3() {
    ifstream file1("data/enwiki-latest-all-titles-in-ns0.sorted");

    string line;
    vector<string> keys1;
    while (getline(file1, line)) {
        keys1.push_back(line);
    }
    file1.close();

    ifstream file2("data/jawiki-latest-all-titles-in-ns0.sorted");
    vector<string> keys2;
    while (getline(file2, line)) {
        keys2.push_back(line);
    }
    file2.close();

    louds::Trie trie1;
    for (uint64_t i = 0; i < keys1.size(); ++i) {
      trie1.add(keys1[i]);
    }
    trie1.build();

    louds::Trie trie2;
    for (uint64_t i = 0; i < keys2.size(); ++i) {
      trie2.add(keys2[i]);
    }
    trie2.build();

    high_resolution_clock::time_point begin = high_resolution_clock::now();
    long before = getMemoryUsage();

    louds::Trie* trie3 = louds::Trie::merge_efficient(trie1, trie2);

    long after = getMemoryUsage();
    std::cout << "Memory used: " << (after - before) * 1024 << " bytes" << std::endl;

    high_resolution_clock::time_point end = high_resolution_clock::now();
    double elapsed = (double)duration_cast<milliseconds>(end - begin).count();
    cout << "Efficient Implementation elapsed time: " << elapsed << endl;

    trie3->build();

    for (uint64_t i = 0; i < keys1.size(); ++i) {
        assert(trie3->lookup(keys1[i]) != -1);
    }

    for (uint64_t i = 0; i < keys2.size(); ++i) {
        assert(trie3->lookup(keys2[i]) != -1);
    }
}

int main() {
  test3();
}
