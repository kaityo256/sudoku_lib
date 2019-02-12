#pragma once
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#ifdef _WIN32
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

#include "../stopwatch/stopwatch.hpp"
#include "mbit.hpp"

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

  int _rest;    // 残りマス
  bool _valid;  // 正常な状態かどうか
  int data[81]; //現在決定している数字
  mbit data_mask;
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
  // 現在の状態が正常かどうか
  bool is_valid(void) {
    if (!_valid) return false;
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
    data_mask = mask81;
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
      if (n == 0) continue;
      if (!can_put(i, n)) {
        _valid = false;
      }
      put(i, n);
    }
  }
  // 9個のビット列のうち、同じ桁で1つだけ立っているビットを返す
  static mbit find_single(const mbit *g) {
    mbit b = g[0] & g[1];
    mbit s = g[0] | g[1];
    for (int i = 2; i < 9; i++) {
      b |= s & g[i];
      s |= g[i];
    }
    return s ^ b;
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
  void put(int pos, int n) {
    cell_mask[n - 1] &= kill_cell_mask[pos];
    mbit mm = mask81 ^ (mbit(1) << pos);
    for (auto &m : cell_mask) {
      m &= mm;
    }
    // kill_column(i, n);
    data[pos] = n;
    data_mask ^= mbit(1) << pos;
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

  const std::string str() const {
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

  static bool is_unique(std::string &str) {
    Grid g(str);
    return g.is_unique();
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

  bool hidden_singles_row(void) {
    bool hit = false;
    mbit mm[9] = {};
    for (int i = 0; i < 9; i++) {
      for (int j = 0; j < 9; j++) {
        mm[i] |= (mbit(1) << (j * 9 + i));
      }
    }
    mbit m_row[9] = {};
    for (int n = 0; n < 9; n++) {
      mbit m = cell_mask[n];
      for (int i = 0; i < 9; i++) {
        m_row[i] |= (((m & mm[i]) >> i) << n);
      }
    }
    mbit gs = Grid::find_single(m_row);
    while (gs) {
      mbit v = gs & -gs;
      int n = bitpos(v) % 9;
      int r = bitpos(v) / 9;
      for (int i = 0; i < 9; i++) {
        if (v & m_row[i]) {
          put(r * 9 + 1, n + 1);
          hit = true;
        }
      }
      gs ^= v;
    }
    return hit;
  }

  bool hidden_singles(void) {
    static stopwatch::timer<> timer("hidden_singles");
    // hidden_singles_row();
    bool hit = false;
    timer.start();
    static const mbit mzero = mbit(0);
    for (int i = 0; i < 9; i++) {
      if (cell_mask[i] == mzero) continue;
      for (auto m : unit_mask) {
        const mbit p = cell_mask[i] & m;
        if ((popcnt_u128(p) == 1)) {
          put(bitpos(p), i + 1);
          // printf("puts %d on %d\n", i + 1, bitpos(p));
          hit = true;
        }
      }
    }
    timer.stop();
    return hit;
  }
  // 現状、Hidden Singlesがあるかどうか返す
  // Naked Singlesのチェックはしていない
  bool has_singles(void) {
    for (int i = 0; i < 9; i++) {
      if (cell_mask[i] == mbit(0)) continue;
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
