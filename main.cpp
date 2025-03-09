#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>

#include "louds-trie.hpp"

using namespace std;
using namespace std::chrono;

int main() {
  /*
  louds::Trie trie;

  vector<string> keys = {"car", "cat", "cap", "cot", "cop", "bap", "bat", "ba"};

  sort(keys.begin(), keys.end());

  for(auto & key: keys) {
    trie.add(key);
  }

  trie.build();

  trie.print();

  vector<string> extracted_keys = trie.extract_keys();

  for (string &s: extracted_keys) {
    cout << s << endl;
  }
  */

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

  louds::Trie* trie3 = louds::merge_optimal(trie1, trie2);

  trie3->build();

  trie3->print();

  vector<string> extracted_keys = trie3->extract_keys();

  cout << "MERGED TRIE" << endl;
  for (string &s: extracted_keys) {
    cout << s << endl;
  }

  louds::Trie* trie4 = louds::merge_naive(trie1, trie2);

  trie4->build();

  trie4->print();

  extracted_keys = trie4->extract_keys();

  cout << "MERGED TRIE" << endl;
  for (string &s: extracted_keys) {
    cout << s << endl;
  }
}
