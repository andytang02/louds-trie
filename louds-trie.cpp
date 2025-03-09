#include "louds-trie.hpp"

#ifdef _MSC_VER
 #include <intrin.h>
 #include <immintrin.h>
#else  // _MSC_VER
 #include <x86intrin.h>
#endif  // _MSC_VER

#include <cassert>
#include <vector>
#include <iostream>
#include <bitset>
#include <cstdint>
#include <string>

namespace louds {
namespace {

uint64_t Popcnt(uint64_t x) {
#ifdef _MSC_VER
  return __popcnt64(x);
#else  // _MSC_VER
  return __builtin_popcountll(x);
#endif  // _MSC_VER
}

uint64_t Ctz(uint64_t x) {
#ifdef _MSC_VER
  return _tzcnt_u64(x);
#else  // _MSC_VER
  return __builtin_ctzll(x);
#endif  // _MSC_VER
}

struct BitVector {
  struct Rank {
    uint32_t abs_hi;
    uint8_t abs_lo;
    uint8_t rels[3];

    uint64_t abs() const {
      return ((uint64_t)abs_hi << 8) | abs_lo;
    }
    void set_abs(uint64_t abs) {
      abs_hi = (uint32_t)(abs >> 8);
      abs_lo = (uint8_t)abs;
    }
  };

  vector<uint64_t> words;
  vector<Rank> ranks;
  vector<uint32_t> selects;
  uint64_t n_bits;

  BitVector() : words(), ranks(), selects(), n_bits(0) {}

  uint64_t get(uint64_t i) const {
    return (words[i / 64] >> (i % 64)) & 1UL;
  }
  void set(uint64_t i, uint64_t bit) {
    if (bit) {
      words[i / 64] |= (1UL << (i % 64));
    } else {
      words[i / 64] &= ~(1UL << (i % 64));
    }
  }

  void add(uint64_t bit) {
    if (n_bits % 256 == 0) {
      words.resize((n_bits + 256) / 64);
    }
    set(n_bits, bit);
    ++n_bits;
  }
  // build builds indexes for rank and select.
  void build() {
    uint64_t n_blocks = words.size() / 4;
    uint64_t n_ones = 0;
    ranks.resize(n_blocks + 1);
    for (uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
      ranks[block_id].set_abs(n_ones);
      for (uint64_t j = 0; j < 4; ++j) {
        if (j != 0) {
          uint64_t rel = n_ones - ranks[block_id].abs();
          ranks[block_id].rels[j - 1] = rel;
        }

        uint64_t word_id = (block_id * 4) + j;
        uint64_t word = words[word_id];
        uint64_t n_pops = Popcnt(word);
        uint64_t new_n_ones = n_ones + n_pops;
        if (((n_ones + 255) / 256) != ((new_n_ones + 255) / 256)) {
          uint64_t count = n_ones;
          while (word != 0) {
            uint64_t pos = Ctz(word);
            if (count % 256 == 0) {
              selects.push_back(((word_id * 64) + pos) / 256);
              break;
            }
            word ^= 1UL << pos;
            ++count;
          }
        }
        n_ones = new_n_ones;
      }
    }
    ranks.back().set_abs(n_ones);
    selects.push_back(words.size() * 64 / 256);
  }

  // rank returns the number of 1-bits in the range [0, i).
  uint64_t rank(uint64_t i) const {
    uint64_t word_id = i / 64;
    uint64_t bit_id = i % 64;
    uint64_t rank_id = word_id / 4;
    uint64_t rel_id = word_id % 4;
    uint64_t n = ranks[rank_id].abs();
    if (rel_id != 0) {
      n += ranks[rank_id].rels[rel_id - 1];
    }
    n += Popcnt(words[word_id] & ((1UL << bit_id) - 1));
    return n;
  }
  // select returns the position of the (i+1)-th 1-bit.
  uint64_t select(uint64_t i) const {
    const uint64_t block_id = i / 256;
    uint64_t begin = selects[block_id];
    uint64_t end = selects[block_id + 1] + 1UL;
    if (begin + 10 >= end) {
      while (i >= ranks[begin + 1].abs()) {
        ++begin;
      }
    } else {
      while (begin + 1 < end) {
        const uint64_t middle = (begin + end) / 2;
        if (i < ranks[middle].abs()) {
          end = middle;
        } else {
          begin = middle;
        }
      }
    }
    const uint64_t rank_id = begin;
    i -= ranks[rank_id].abs();

    uint64_t word_id = rank_id * 4;
    if (i < ranks[rank_id].rels[1]) {
      if (i >= ranks[rank_id].rels[0]) {
        word_id += 1;
        i -= ranks[rank_id].rels[0];
      }
    } else if (i < ranks[rank_id].rels[2]) {
      word_id += 2;
      i -= ranks[rank_id].rels[1];
    } else {
      word_id += 3;
      i -= ranks[rank_id].rels[2];
    }
    return (word_id * 64) + Ctz(_pdep_u64(1UL << i, words[word_id]));
  }

  uint64_t size() const {
    return sizeof(uint64_t) * words.size()
      + sizeof(Rank) * ranks.size()
      + sizeof(uint32_t) * selects.size();
  }
};

struct Level {
  BitVector louds;
  BitVector outs;
  vector<uint8_t> labels;
  uint64_t offset;

  Level() : louds(), outs(), labels(), offset(0) {}

  uint64_t size() const;

  void print();
};

void Level::print() {
  cout << louds.n_bits << endl;
  for (uint64_t x: louds.words) {
    bitset<64> b(x);
    cout << b << " " << endl;
  }

  cout << outs.n_bits << endl;
  for (uint64_t x: outs.words) {
    bitset<64> b(x);
    cout << b << " " << endl;
  }

  for (uint8_t l: labels) {
    cout << (char) l << " ";
  }
  cout << endl;
}

uint64_t Level::size() const {
  return louds.size() + outs.size() + labels.size();
}

} // namespace

class TrieImpl {
 public:
  TrieImpl();
  ~TrieImpl() {}

  void add(const string &key);
  void build();
  int64_t lookup(const string &query) const;

  uint64_t n_keys() const {
    return n_keys_;
  }
  uint64_t n_nodes() const {
    return n_nodes_;
  }
  uint64_t size() const {
    return size_;
  }

  void print();

  void extract_keys(uint64_t level_i, uint64_t node_id, string prefix, vector<string> &keys);

 private:
  vector<Level> levels_;
  uint64_t n_keys_;
  uint64_t n_nodes_;
  uint64_t size_;
  string last_key_;
};

void TrieImpl::print() {
  for (auto &l: levels_) {
    l.print();
  }
}

TrieImpl::TrieImpl()
  : levels_(2), n_keys_(0), n_nodes_(1), size_(0), last_key_() {
  levels_[0].louds.add(0);
  levels_[0].louds.add(1);
  levels_[1].louds.add(1);
  levels_[0].outs.add(0);
  levels_[0].labels.push_back(' ');
}

void TrieImpl::add(const string &key) {
  assert(key > last_key_);
  if (key.empty()) {
    levels_[0].outs.set(0, 1);
    ++levels_[1].offset;
    ++n_keys_;
    return;
  }
  if (key.length() + 1 >= levels_.size()) {
    levels_.resize(key.length() + 2);
  }
  uint64_t i = 0;
  for ( ; i < key.length(); ++i) {
    Level &level = levels_[i + 1];
    uint8_t byte = key[i];
    if ((i == last_key_.length()) || (byte != level.labels.back())) {
      level.louds.set(levels_[i + 1].louds.n_bits - 1, 0);
      level.louds.add(1);
      level.outs.add(0);
      level.labels.push_back(key[i]);
      ++n_nodes_;
      break;
    }
  }
  for (++i; i < key.length(); ++i) {
    Level &level = levels_[i + 1];
    level.louds.add(0);
    level.louds.add(1);
    level.outs.add(0);
    level.labels.push_back(key[i]);
    ++n_nodes_;
  }
  levels_[key.length() + 1].louds.add(1);
  ++levels_[key.length() + 1].offset;
  levels_[key.length()].outs.set(levels_[key.length()].outs.n_bits - 1, 1);
  ++n_keys_;
  last_key_ = key;
}

void TrieImpl::build() {
  uint64_t offset = 0;
  for (uint64_t i = 0; i < levels_.size(); ++i) {
    Level &level = levels_[i];
    level.louds.build();
    level.outs.build();
    offset += levels_[i].offset;
    level.offset = offset;
    size_ += level.size();
  }
}

int64_t TrieImpl::lookup(const string &query) const {
  if (query.length() >= levels_.size()) {
    return false;
  }
  uint64_t node_id = 0;
  for (uint64_t i = 0; i < query.length(); ++i) {
    const Level &level = levels_[i + 1];
    uint64_t node_pos;
    if (node_id != 0) {
      node_pos = level.louds.select(node_id - 1) + 1;
      node_id = node_pos - node_id;
    } else {
      node_pos = 0;
    }

    // Linear search implementation
    // for (uint8_t byte = query[i]; ; ++node_pos, ++node_id) {
    //   if (level.louds.get(node_pos) || level.labels[node_id] > byte) {
    //     return -1;
    //   }
    //   if (level.labels[node_id] == byte) {
    //     break;
    //   }
    // }

    // Binary search implementation
    uint64_t end = node_pos;
    uint64_t word = level.louds.words[end / 64] >> (end % 64);
    if (word == 0) {
      end += 64 - (end % 64);
      word = level.louds.words[end / 64];
      while (word == 0) {
        end += 64;
        word = level.louds.words[end / 64];
      }
    }
    end += Ctz(word);
    uint64_t begin = node_id;
    end = begin + end - node_pos;

    uint8_t byte = query[i];
    while (begin < end) {
      node_id = (begin + end) / 2;
      if (byte < level.labels[node_id]) {
        end = node_id;
      } else if (byte > level.labels[node_id]) {
        begin = node_id + 1;
      } else {
        break;
      }
    }
    if (begin >= end) {
      return -1;
    }
  }
  const Level &level = levels_[query.length()];
  if (!level.outs.get(node_id)) {
    return false;
  }
  return level.offset + level.outs.rank(node_id);
}

void TrieImpl::extract_keys(uint64_t level_i, uint64_t node_id, string prefix, vector<string> &keys) {
  if (prefix.length() >= levels_.size() - 2) {
    return;
  }

  const Level &level = levels_[level_i];

  uint64_t node_pos;
  if (node_id != 0) {
    node_pos = level.louds.select(node_id - 1) + 1;
    node_id = node_pos - node_id;
  } else {
    node_pos = 0;
  }

  // Caculate [begin, end) of children
  uint64_t end = node_pos;
  uint64_t word = level.louds.words[end / 64] >> (end % 64);
  if (word == 0) {
    end += 64 - (end % 64);
    word = level.louds.words[end / 64];
    while (word == 0) {
      end += 64;
      word = level.louds.words[end / 64];
    }
  }
  end += Ctz(word);
  uint64_t begin = node_id;
  end = begin + end - node_pos;

  for (uint64_t i = begin; i < end; i ++) {
    string next = prefix + (char) level.labels[i];
    if (level.outs.get(i) == 1) {
      keys.push_back(next);
    }
    extract_keys(level_i + 1, i, next, keys);
  }
}

Trie::Trie() : impl_(new TrieImpl) {}

Trie::~Trie() {
  delete impl_;
}

void Trie::add(const string &key) {
  return impl_->add(key);
}

void Trie::build() {
  impl_->build();
}

int64_t Trie::lookup(const string &query) const {
  return impl_->lookup(query);
}

uint64_t Trie::n_keys() const {
  return impl_->n_keys();
}

uint64_t Trie::n_nodes() const {
  return impl_->n_nodes();
}

uint64_t Trie::size() const {
  return impl_->size();
}

void Trie::print() {
  impl_->print();
}

vector<string> Trie::extract_keys() const {
  vector<string> keys;

  impl_->extract_keys(1, 0, "", keys);

  return keys;
}

Trie* merge(const Trie& trie1, const Trie& trie2) {
  Trie* merged_trie = new Trie();

  vector<string> keys1 = trie1.extract_keys();
  vector<string> keys2 = trie2.extract_keys();

  cout << "TRIE 1" << endl;
  for (string & s1: keys1) {
    cout << s1 << endl;
  }

  cout << "TRIE 2" << endl;
  for (string & s2: keys2) {
    cout << s2 << endl;
  }

  uint64_t i = 0;
  uint64_t j = 0;
  string last_string = "";

  while (i < keys1.size() && j < keys2.size()) {
    if (keys1[i] < keys2[j]) {
      if (keys1[i] == last_string) {
        i ++;
      } else {
        last_string = keys1[i ++];
        merged_trie->add(last_string);
      }
    } else {
      if (keys2[j] == last_string) {
        j ++;
      } else {
        last_string = keys2[j ++];
        merged_trie->add(last_string);
      }
    }
  }

  while (i < keys1.size()) {
    if (keys1[i] == last_string) {
      i ++;
    } else {
      last_string = keys1[i ++];
      merged_trie->add(last_string);
    }
  }

  while (j < keys2.size()) {
    if (keys2[j] == last_string) {
      j ++;
    } else {
      last_string = keys2[j ++];
      merged_trie->add(last_string);
    }
  }

  merged_trie->build();

  return merged_trie;
} 

}  // namespace louds

#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>

// #include "louds-trie.hpp"

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
  vector<string> keys1 = {"b", "ba", "bat", "bar"};
  sort(keys1.begin(), keys1.end());

  for(auto & key: keys1) {
    trie1.add(key);
  }

  trie1.build();

  louds::Trie trie2;
  vector<string> keys2 = {"c", "ca", "cat", "ba", "car"};
  sort(keys2.begin(), keys2.end());

  for(auto & key: keys2) {
    trie2.add(key);
  }

  trie2.build();

  louds::Trie* trie3 = louds::merge(trie1, trie2);

  vector<string> extracted_keys = trie3->extract_keys();

  cout << "MERGED TRIE" << endl;
  for (string &s: extracted_keys) {
    cout << s << endl;
  }
}
