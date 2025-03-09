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
  cout << "Louds Bits: " << louds.n_bits << endl;
  for (uint64_t x: louds.words) {
    bitset<64> b(x);
    cout << b << " " << endl;
  }

  cout << "Outs Bits: " << outs.n_bits << endl;
  for (uint64_t x: outs.words) {
    bitset<64> b(x);
    cout << b << " " << endl;
  }

  for (uint8_t l: labels) {
    cout << (char) l << " ";
  }
  cout << endl;

  cout << "Offset: " << offset << endl;
}

uint64_t Level::size() const {
  return louds.size() + outs.size() + labels.size();
}

} // namespace

class TrieImpl {
 public:
  TrieImpl();
  TrieImpl(bool initialize_first);
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

  static TrieImpl* merge_optimal(TrieImpl& trie1, TrieImpl &trie2);

 private:
  vector<Level> levels_;
  uint64_t n_keys_;
  uint64_t n_nodes_;
  uint64_t size_;
  string last_key_;
};

void TrieImpl::print() {
  cout << "NUM KEYS: " << n_keys_ << endl;
  cout << "NUM NODES: " << n_nodes_ << endl;
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

TrieImpl::TrieImpl(bool initialize_first)
: levels_(2), n_keys_(0), n_nodes_(1), size_(0), last_key_() {
  levels_[0].louds.add(0);
  levels_[0].louds.add(1);
  if (initialize_first) {
    levels_[1].louds.add(1);
  }
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

TrieImpl* TrieImpl::merge_optimal(TrieImpl &trie1, TrieImpl &trie2) {
  TrieImpl* merged_trie = new TrieImpl(false);

  if (trie1.levels_.size() > trie2.levels_.size()) {
    swap(trie1, trie2);
  }

  merged_trie->levels_.resize(trie2.levels_.size());

  vector<string> prefix1 = {""};
  vector<string> prefix2 = {""};

  vector<string> next_prefix1;
  vector<string> next_prefix2;

  for (uint64_t i = 0; i < trie1.levels_.size() - 1; i ++) {
    // Level &par_level1 = trie1.levels_[i];
    // Level &par_level2 = trie2.levels_[i];

    Level &level1 = trie1.levels_[i+1];
    Level &level2 = trie2.levels_[i+1];

    uint64_t j1 = 0;
    uint64_t j2 = 0;
    
    for (string &s1: prefix1) {
      cout << s1 << " ";
    }
    cout << endl;

    for (string &s2: prefix2) {
      cout << s2 << " ";
    }
    cout << endl;

    while (j1 < prefix1.size() && j2 < prefix2.size()) {
      if (prefix1[j1] == prefix2[j2]) {
        uint64_t start1 = j1 == 0 ? 0 : level1.louds.select(j1 - 1) + 1 - j1;
        uint64_t end1 = level1.louds.select(j1) - j1;

        uint64_t start2 = j2 == 0 ? 0 : level2.louds.select(j2 - 1) + 1 - j2;
        uint64_t end2 = level2.louds.select(j2) - j2;

        uint64_t k1 = start1;
        uint64_t k2 = start2;

        // cout << start1 << " " << end1 << endl;

        while (k1 < end1 && k2 < end2) {
          if (level1.labels[k1] < level2.labels[k2]) {
            merged_trie->levels_[i+1].louds.add(0);
            merged_trie->levels_[i+1].outs.add(level1.outs.get(k1));
            merged_trie->levels_[i+1].labels.push_back(level1.labels[k1]);
  
            next_prefix1.push_back(prefix1[j1] + (char) level1.labels[k1]);

            merged_trie->n_nodes_ ++;
            if (level1.outs.get(k1) == 1) {
              merged_trie->n_keys_ ++;
              merged_trie->levels_[i+2].offset ++;
            }

            k1 ++;
          }
          else if(level1.labels[k1] > level2.labels[k2]) {
            merged_trie->levels_[i+1].louds.add(0);
            merged_trie->levels_[i+1].outs.add(level2.outs.get(k2));
            merged_trie->levels_[i+1].labels.push_back(level2.labels[k2]);
  
            next_prefix2.push_back(prefix2[j2] + (char) level2.labels[k2]);

            merged_trie->n_nodes_ ++;
            if (level2.outs.get(k2) == 1) {
              merged_trie->n_keys_ ++;
              merged_trie->levels_[i+2].offset ++;
            }

            k2 ++;
          }
          else {
            merged_trie->levels_[i+1].louds.add(0);
            merged_trie->levels_[i+1].outs.add(level1.outs.get(k1) | level2.outs.get(k2));
            merged_trie->levels_[i+1].labels.push_back(level1.labels[k1]);

            next_prefix1.push_back(prefix1[j1] + (char) level1.labels[k1]);
            next_prefix2.push_back(prefix2[j2] + (char) level2.labels[k2]);

            merged_trie->n_nodes_ ++;
            if ((level1.outs.get(k1) | level2.outs.get(k2)) == 1) {
              merged_trie->n_keys_ ++;
              merged_trie->levels_[i+2].offset ++;
            }

            k1 ++;
            k2 ++;
          }
        }

        while (k1 < end1) {
          merged_trie->levels_[i+1].louds.add(0);
          merged_trie->levels_[i+1].outs.add(level1.outs.get(k1));
          merged_trie->levels_[i+1].labels.push_back(level1.labels[k1]);

          next_prefix1.push_back(prefix1[j1] + (char) level1.labels[k1]);

          merged_trie->n_nodes_ ++;
          if (level1.outs.get(k1) == 1) {
            merged_trie->n_keys_ ++;
            merged_trie->levels_[i+2].offset ++;
          }

          k1 ++;
        }

        while (k2 < end2) {
          merged_trie->levels_[i+1].louds.add(0);
          merged_trie->levels_[i+1].outs.add(level2.outs.get(k2));
          merged_trie->levels_[i+1].labels.push_back(level2.labels[k2]);

          next_prefix2.push_back(prefix2[j2] + (char) level2.labels[k2]);
          
          merged_trie->n_nodes_ ++;
          if (level2.outs.get(k2) == 1) {
            merged_trie->n_keys_ ++;
            merged_trie->levels_[i+2].offset ++;
          }

          k2 ++;
        }

        merged_trie->levels_[i+1].louds.add(1);

        j1 ++;
        j2 ++;
      }
      else if (prefix1[j1] < prefix2[j2]) {
        uint64_t start = j1 == 0 ? 0 : level1.louds.select(j1 - 1) + 1 - j1;
        uint64_t end = level1.louds.select(j1) - j1;

        for (uint64_t k = start; k < end; k ++) {
          merged_trie->levels_[i+1].louds.add(0);
          merged_trie->levels_[i+1].outs.add(level1.outs.get(k));
          merged_trie->levels_[i+1].labels.push_back(level1.labels[k]);

          merged_trie->n_nodes_ ++;
          if (level1.outs.get(k) == 1) {
            merged_trie->n_keys_ ++;
            merged_trie->levels_[i+2].offset ++;
          }

          next_prefix1.push_back(prefix1[j1] + (char) level1.labels[k]);
        }

        merged_trie->levels_[i+1].louds.add(1);

        j1 ++;
      }
      else {
        uint64_t start = j2 == 0 ? 0 : level2.louds.select(j2 - 1) + 1 - j2;
        uint64_t end = level2.louds.select(j2) - j2;

        for (uint64_t k = start; k < end; k ++) {
          merged_trie->levels_[i+1].louds.add(0);
          merged_trie->levels_[i+1].outs.add(level2.outs.get(k));
          merged_trie->levels_[i+1].labels.push_back(level2.labels[k]);

          merged_trie->n_nodes_ ++;
          if (level2.outs.get(k) == 1) {
            merged_trie->n_keys_ ++;
            merged_trie->levels_[i+2].offset ++;
          }

          next_prefix2.push_back(prefix2[j2] + (char) level2.labels[k]);
        }

        merged_trie->levels_[i+1].louds.add(1);

        j2 ++;
      }
    }

    while (j1 < prefix1.size()) {
      uint64_t start = j1 == 0 ? 0 : level1.louds.select(j1 - 1) + 1 - j1;
      uint64_t end = level1.louds.select(j1) - j1;

      for (uint64_t k = start; k < end; k ++) {
        merged_trie->levels_[i+1].louds.add(0);
        merged_trie->levels_[i+1].outs.add(level1.outs.get(k));
        merged_trie->levels_[i+1].labels.push_back(level1.labels[k]);

        merged_trie->n_nodes_ ++;
        if (level1.outs.get(k) == 1) {
          merged_trie->n_keys_ ++;
          merged_trie->levels_[i+2].offset ++;
        }

        next_prefix1.push_back(prefix1[j1] + (char) level1.labels[k]);
      }

      merged_trie->levels_[i+1].louds.add(1);

      j1 ++;
    }

    while (j2 < prefix2.size()) {
      uint64_t start = j2 == 0 ? 0 : level2.louds.select(j2 - 1) + 1 - j2;
      uint64_t end = level2.louds.select(j2) - j2;

      for (uint64_t k = start; k < end; k ++) {
        merged_trie->levels_[i+1].louds.add(0);
        merged_trie->levels_[i+1].outs.add(level2.outs.get(k));
        merged_trie->levels_[i+1].labels.push_back(level2.labels[k]);

        merged_trie->n_nodes_ ++;
        if (level2.outs.get(k) == 1) {
          merged_trie->n_keys_ ++;
          merged_trie->levels_[i+2].offset ++;
        }

        next_prefix2.push_back(prefix2[j2] + (char) level2.labels[k]);
      }

      merged_trie->levels_[i+1].louds.add(1);

      j2 ++;
    }

    prefix1 = next_prefix1;
    prefix2 = next_prefix2;

    next_prefix1.clear();
    next_prefix2.clear();
  }


  for (uint64_t i = trie1.levels_.size() - 1; i < trie2.levels_.size() - 1; i ++) {
    Level &level2 = trie2.levels_[i+1];

    uint64_t j2 = 0;

    for (string &s1: prefix1) {
      cout << s1 << " ";
    }
    cout << endl;
  
    for (string &s2: prefix2) {
      cout << s2 << " ";
    }
    cout << endl;

    while (j2 < prefix2.size()) {
      uint64_t start = j2 == 0 ? 0 : level2.louds.select(j2 - 1) + 1 - j2;
      uint64_t end = level2.louds.select(j2) - j2;

      for (uint64_t k = start; k < end; k ++) {
        merged_trie->levels_[i+1].louds.add(0);
        merged_trie->levels_[i+1].outs.add(level2.outs.get(k));
        merged_trie->levels_[i+1].labels.push_back(level2.labels[k]);

        merged_trie->n_nodes_ ++;
        if (level2.outs.get(k) == 1) {
          merged_trie->n_keys_ ++;
          merged_trie->levels_[i+2].offset ++;
        }

        next_prefix2.push_back(prefix2[j2] + (char) level2.labels[k]);
      }

      merged_trie->levels_[i+1].louds.add(1);

      j2 ++;
    }

    prefix2 = next_prefix2;

    next_prefix2.clear();
  }

  for (string &s1: prefix1) {
    cout << s1 << " ";
  }
  cout << endl;

  for (string &s2: prefix2) {
    cout << s2 << " ";
  }
  cout << endl;

  /*
  for (uint64_t i = 0; i < merged_trie->levels_[merged_trie->levels_.size() - 2].outs.n_bits; i ++) {
    merged_trie->levels_[merged_trie->levels_.size() - 1].louds.add(1);
  }
  */

  return merged_trie;
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

  // uint64_t original_id = node_id;

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

  // cout << end << " " << level.louds.select(original_id) - original_id << endl;

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

Trie* merge_naive(const Trie& trie1, const Trie& trie2) {
  Trie* merged_trie = new Trie();

  vector<string> keys1 = trie1.extract_keys();
  vector<string> keys2 = trie2.extract_keys();

  /*
  cout << "TRIE 1" << endl;
  for (string & s1: keys1) {
    cout << s1 << endl;
  }

  cout << "TRIE 2" << endl;
  for (string & s2: keys2) {
    cout << s2 << endl;
  }
  */

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

  // merged_trie->build();

  return merged_trie;
} 

Trie* merge_optimal(const Trie& trie1, const Trie& trie2) {
  Trie* merged_trie = new Trie();

  merged_trie->impl_ = TrieImpl::merge_optimal(*trie1.impl_,*trie2.impl_);

  return merged_trie;
}

}  // namespace louds