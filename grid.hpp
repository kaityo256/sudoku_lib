#pragma once
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#ifdef _WIN32
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

typedef unsigned __int128 mbit;

const static mbit mask81 = (mbit(1) << 81) - 1;
inline std::ostream &operator<<(std::ostream &dest, mbit v) {
  char buffer[81];
  for (int i = 0; i < 81; i++) {
    buffer[i] = (v & 1) + '0';
    v = v >> 1;
  }
  dest.rdbuf()->sputn(buffer, 81);
  return dest;
}

inline mbit str2mbit(const std::string &str) {
  mbit t = 0;
  for (unsigned int i = 0; i < 81 && i < str.size(); i++) {
    t ^= mbit((str[i] != '0')) << i;
  }
  return t;
}

inline int popcnt_u128(const mbit &s) {
  const uint64_t *m = (uint64_t *)&s;
  return _mm_popcnt_u64(m[0]) + _mm_popcnt_u64(m[1]);
}

inline int bitpos(mbit s) {
  return popcnt_u128(s - 1);
}

class Grid {
private:
public:
  // 場所とボックスインデックス
  constexpr static int box_index[81] = {
      0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 0, 0,
      1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 3, 3, 3, 4, 4, 4,
      5, 5, 5, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8,
      6, 6, 6, 7, 7, 7, 8, 8, 8, 6, 6, 6, 7, 7, 7, 8, 8, 8};

  // 行、列、ボックスの位置を表すマスク
  static mbit unit_mask[27];

  //そのマスに置かれたら、どのcell_maskを消すべきかANDするマスク
  static mbit kill_cell_mask[81];

  // マスクの自動初期化用クラス
  class GridInitializer {
  public:
    GridInitializer() {
      Grid::init_masks();
    };
  };
  static GridInitializer si;

  int _rest;         // 残りマス
  bool _valid;       // 正常な状態かどうか
  int data[81];      //現在決定している数字
  mbit cell_mask[9]; // 各数字ごとにおける可能性のあるセル

  /*
   * ユニットごとにどの場所に置けるかを表現するマスク
   * 引数は場所
   * 1がおける行が9ビット、2がおける行が9ビット...と81ビット使う
   */
  mbit row_mask[9];
  mbit column_mask[9];
  mbit box_mask[9];

  bool solved_squares(void); // Naked Singles
  bool hidden_singles(void); // Hidden Singles
  // 現在の状態が正常かどうか
  bool is_valid(void) {
    if (!_valid)
      return false;
    mbit data_mask = mbit(0);
    for (int i = 0; i < 81; i++) {
      if (data[i] == 0) {
        data_mask |= mbit(1) << i;
      }
    }
    mbit r = mbit(0);
    for (auto &m : cell_mask) {
      r |= m;
    }
    return r == data_mask;
  }

  static void init_masks(void);

  void init() {
    std::fill(&cell_mask[0], &cell_mask[9], mask81);
    std::fill(&column_mask[0], &column_mask[9], mask81);
    std::fill(&data[0], &data[81], 0);
    _rest = 81;
    _valid = true;
  }

  Grid() {
    init();
  }

  // その数字の列マスクを削除する
  void kill_column(int index, int num) {
    num = num - 1;
    int ri = index / 9;
    int ci = index % 9;
    mbit mask1 = mask81 ^ (((mbit(1) << 9) - 1) << (num * 9));
    column_mask[ri] &= mask1;
    mbit mask2 = mask81 ^ (mbit(1) << (num * 9 + ci));
    for (int i = 0; i < 9; i++) {
      column_mask[i] &= mask2;
    }
    for (int n = 0; n < 9; n++) {
      column_mask[ri] &= mask81 ^ (mbit(1) << (n * 9 + ci));
    }
    // ボックスを消す
    int bx = (ci / 3) * 3;
    int by = (ri / 3) * 3;
    mbit mask3 = mask81 ^ ((mbit(1) << 3) - 1) << (num * 9 + bx);
    column_mask[by + 0] &= mask3;
    column_mask[by + 1] &= mask3;
    column_mask[by + 2] &= mask3;
  }

  Grid(const std::string &str) {
    init();
    for (int i = 0; i < 81; i++) {
      const int n = str[i] - '0';
      if (n == 0)
        continue;
      if (!can_put(i, n)) {
        _valid = false;
      }
      put(i, n);
    }
  }

  mbit find_single2(void) {
    const mbit *g = cell_mask;
    const mbit g01a = g[0] ^ g[1];
    const mbit g01b = ~(g[0] | g[1]);
    const mbit g23a = g[2] ^ g[3];
    const mbit g23b = ~(g[2] | g[3]);
    const mbit g03a = (g01a & g23b) | (g01b & g23a);
    const mbit g03b = g01b & g23b;
    const mbit g45a = g[4] ^ g[5];
    const mbit g45b = ~(g[4] | g[5]);
    const mbit g67a = g[6] ^ g[7];
    const mbit g67b = ~(g[6] | g[7]);
    const mbit g47a = (g45a & g67b) | (g45b & g67a);
    const mbit g47b = g45b & g67b;
    const mbit g07a = (g03a & g47b) | (g03b & g47a);
    const mbit g07b = g03b & g47b;
    mbit b = (g07a & ~g[8]) | (g07b & g[8]);
    return b;
  }

  mbit find_single_kawai2(void) {
    const mbit *g = cell_mask;
    mbit ba = g[0] & g[1];
    mbit sa = g[0] | g[1];
    ba |= sa & g[2];
    sa |= g[2];
    ba |= sa & g[3];
    sa |= g[3];

    mbit bb = g[4] & g[5];
    mbit sb = g[4] | g[5];
    bb |= sb & g[6];
    sb |= g[6];
    bb |= sb & g[7];
    sb |= g[7];
    bb |= sb & g[8];
    sb |= g[8];
    ba |= bb | (sa & sb);
    sa |= sb;
    return sa ^ ba;
  }
  mbit find_single_kawai3(void) {
    const mbit *g = cell_mask;
    mbit b = g[0] & g[1];
    mbit s = g[0] | g[1];
    for (int i = 2; i < 9; i++) {
      b |= s & g[i];
      s |= g[i];
    }
    return s ^ b;
  }

  mbit find_single_kawai(void) {
    const mbit *g = cell_mask;
    mbit b = g[0] & g[1];
    mbit s = g[0] | g[1];
    b |= s & g[2];
    s |= g[2];
    b |= s & g[3];
    s |= g[3];
    b |= s & g[4];
    s |= g[4];
    b |= s & g[5];
    s |= g[5];
    b |= s & g[6];
    s |= g[6];
    b |= s & g[7];
    s |= g[7];
    b |= s & g[8];
    s |= g[8];
    return s ^ b;
  }

  mbit find_single_org(void) {
    const mbit *g = cell_mask;
    const mbit x1 = (g[0] ^ g[1] ^ g[2]);
    const mbit x2 = (g[3] ^ g[4] ^ g[5]);
    const mbit x3 = (g[6] ^ g[7] ^ g[8]);
    const mbit a1 = (g[0] & g[1] & g[2]);
    const mbit a2 = (g[3] & g[4] & g[5]);
    const mbit a3 = (g[6] & g[7] & g[8]);
    const mbit o1 = (g[0] | g[1] | g[2]);
    const mbit o2 = (g[3] | g[4] | g[5]);
    const mbit o3 = (g[6] | g[7] | g[8]);
    const mbit b1 = x1 ^ a1;
    const mbit b2 = x2 ^ a2;
    const mbit b3 = x3 ^ a3;
    const mbit c1 = b1 & (mask81 ^ (o2 | o3));
    const mbit c2 = b2 & (mask81 ^ (o3 | o1));
    const mbit c3 = b3 & (mask81 ^ (o1 | o2));
    mbit b = c1 | c2 | c3;
    return b;
  }

  // セル内二択を探す
  mbit find_two(void) {
    mbit *b = cell_mask;
    // 0..3が1ビット
    mbit b0_3a = (b[0] ^ b[1]) & ~(b[2] | b[3]);
    mbit b0_3b = (b[2] ^ b[3]) & ~(b[0] | b[1]);
    mbit b0_3_1 = b0_3a | b0_3b;

    // 0..3が2ビット
    mbit b0_3_0 = ~(b[0] | b[1] | b[2] | b[3]); // 0-3が0
    mbit b0_3c = ~(b[0] ^ b[1] ^ b[2] ^ b[3]);  // 0か2か4
    mbit b0_3d = (b[0] & b[1] & b[2] & b[3]);   // 4
    mbit b0_3_2 = b0_3c & (~b0_3d) & (~b0_3_0);

    // 4..7が1ビット
    mbit b4_7a = (b[4] ^ b[5]) & ~(b[6] | b[7]);
    mbit b4_7b = (b[6] ^ b[7]) & ~(b[4] | b[5]);
    mbit b4_7_1 = b4_7a | b4_7b;

    // 4..7が2ビット
    mbit b4_7_0 = ~(b[4] | b[5] | b[6] | b[7]); // 0-3が0
    mbit b4_7c = ~(b[4] ^ b[5] ^ b[6] ^ b[7]);  // 0か2か4
    mbit b4_7d = (b[4] & b[5] & b[6] & b[7]);   // 4
    mbit b4_7_2 = b4_7c & (~b4_7d) & (~b4_7_0);

    mbit b2 = (b0_3_2 & b4_7_0 & ~b[8])    // 2, 0, 0
              | (b0_3_1 & b4_7_1 & ~b[8])  // 1,1,0
              | (b0_3_1 & b4_7_0 & b[8])   // 1,0,1
              | (b0_3_0 & b4_7_1 & b[8])   // 0,1,1
              | (b0_3_0 & b4_7_2 & ~b[8]); // 0,2,0

    return b2;
  }

  //その場所で置ける数字のリストを返す
  std::vector<int> get_possibles(int index) {
    std::vector<int> v;
    mbit m = mbit(1) << index;
    for (int i = 0; i < 9; i++) {
      if (cell_mask[i] & m) {
        v.push_back(i + 1);
      }
    }
    return v;
  }

  // 数字をマスに置き、マスクの対応するビットを削除
  void put(int i, int n) {
    cell_mask[n - 1] &= kill_cell_mask[i];
    mbit mm = mask81 ^ (mbit(1) << i);
    for (auto &m : cell_mask) {
      m &= mm;
    }
    // kill_column(i, n);
    data[i] = n;
    _rest--;
  }

  // 現在のマスクの表示
  void show_mask() {
    for (int i = 0; i < 9; i++) {
      std::cout << cell_mask[i] << std::endl;
    }
    std::cout << std::endl;
  }

  //現在のユニットマスクの表示
  void show_unit_mask() {
    for (int i = 0; i < 9; i++) {
      std::cout << column_mask[i] << std::endl;
    }
    std::cout << std::endl;
  }

  // その場所にその数字がおけるか
  bool can_put(const int index, const int n) {
    return (cell_mask[n - 1] & (mbit(1) << index));
  }

  unsigned int solve_internal(std::string &answer);
  unsigned int solve_unit(std::string &answer);
  static void solve(std::string &str);

  std::string get_data() {
    std::stringstream ss;
    for (int i = 0; i < 81; i++) {
      ss << data[i];
    }
    return ss.str();
  }

  bool is_unique() {
    std::string ans;
    int n = solve_internal(ans);
    return (n == 1);
  }

  // 指定された数字以下についてユニット内二択の数を返す
  int find_alt_unit(int v) {
    int sum = 0;
    for (int n = 0; n < v; n++) {
      for (int i = 0; i < 27; i++) {
        mbit m = cell_mask[n] & unit_mask[i];
        if (popcnt_u128(m) == 2) {
          sum++;
        }
      }
    }
    return sum;
  }

  // 現状、Hidden Singlesがあるかどうか返す
  // Naked Singlesのチェックはしていない
  bool has_singles(void) {
    for (int i = 0; i < 9; i++) {
      if (cell_mask[i] == mbit(0))
        continue;
      for (const auto &m : unit_mask) {
        const mbit p = cell_mask[i] & m;
        if ((popcnt_u128(p) == 1)) {
          return true;
        }
      }
    }
    return false;
  }

  void show() {
    for (int i = 0; i < 81; i++) {
      std::cout << data[i];
    }
    std::cout << std::endl;
  }
};
