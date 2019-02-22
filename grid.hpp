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

  //その数字がその場所に置かれたら、ユニットマスクのどこを消すべきかANDするマスク
  // 数字、場所、ユニット内の場所
  static mbit kill_row_mask[9][81][9];
  static mbit kill_column_mask[9][81][9];
  static mbit kill_box_mask[9][81][9];

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
  mbit cell_mask[9];            // 各数字ごとにおける可能性のあるセル(naked singles)
  mbit remained_row_mask[9];    // 各行ごとのユニットマスク(hidden_singles)
  mbit remained_column_mask[9]; // 各行ごとのユニットマスク(hidden_singles)
  mbit remained_box_mask[9];    // 各行ごとのユニットマスク(hidden_singles)

  /*
   * ユニットごとにどの場所に置けるかを表現するマスク
   * 引数はユニット内の位置
   * ヒント数字ごとに9ビットずつ使う
   * row_mask[3]は、三行目のマスク。最初の9ビットがヒント数字1の、次の9ビットが2のマスク
   */
  mbit row_mask[9];
  mbit column_mask[9];
  mbit box_mask[9];

#if 1
  bool solved_squares(void); // Naked Singles
#else
  bool solved_squares(void) {
    mbit b = find_single(cell_mask);
    bool flag = false;
    while (b) {
      const mbit p = (b & -b);
      for (int i = 0; i < 9; i++) {
        if (p & cell_mask[i]) {
          put(bitpos(p), i + 1);
          flag = true;
        }
      }
      b ^= p;
    }
    return flag;
  }
#endif

  // 現在の状態が正常かどうか
  bool is_valid(void) {
    if (!_valid)
      return false;
    mbit r = mbit(0);
    for (auto &m : cell_mask) {
      r |= m;
    }
    return r == data_mask;
  }

  static void init_masks(void);

  void init() {
    std::fill(&cell_mask[0], &cell_mask[9], mask81);
    std::fill(&remained_row_mask[0], &remained_row_mask[9], mask81);
    std::fill(&remained_column_mask[0], &remained_column_mask[9], mask81);
    std::fill(&remained_box_mask[0], &remained_box_mask[9], mask81);
    std::fill(&data[0], &data[81], 0);
    data_mask = mask81;
    _rest = 81;
    _valid = true;
  }

  Grid() {
    init();
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

  // 二択になっている場所を返す
  mbit find_two(mbit *b) {
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

  void kill_unit_mask(int pos, int n) {
    for (int i = 0; i < 9; i++) {
      remained_row_mask[i] &= kill_row_mask[n][pos][i];
      remained_column_mask[i] &= kill_column_mask[n][pos][i];
      remained_box_mask[i] &= kill_box_mask[n][pos][i];
    }
  }

  // 数字をマスに置き、マスクの対応するビットを削除
  void put(int pos, int n) {
    cell_mask[n - 1] &= kill_cell_mask[pos];
    mbit mm = mask81 ^ (mbit(1) << pos);
    for (auto &m : cell_mask) {
      m &= mm;
    }
    kill_unit_mask(pos, n - 1);

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

  // その場所にその数字がおけるか
  bool can_put(unsigned int pos, const int n) {
    return (cell_mask[n - 1] & (mbit(1) << pos));
  }

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

  // TODO: 現在のcell_maskからkill maskを作る
  // m_rowなどを受け取りにする
  // 後でリファクタリング
  // m_row などは0初期化されていることが想定されている
  void get_kill_mask(mbit m_row[9], mbit m_column[9], mbit m_box[9]) {
    for (int n = 0; n < 9; n++) {
      for (int i = 0; i < 81; i++) {
        if (cell_mask[n] & (mbit(1) << i)) {
          int r = i % 9;
          int c = i / 9;
          int br = r % 3;
          int bc = c % 3;
          int b = br + bc * 3;            //ボックス内インデックス
          int bi = (r / 3) + (c / 3) * 3; //ボックスのインデックス
          m_row[r] |= mbit(1) << (c + n * 9);
          m_column[c] |= mbit(1) << (r + n * 9);
          m_box[b] |= mbit(1) << (bi + n * 9);
        }
      }
    }
  }

  // マスクによるhidden singlesの探索
  bool hidden_singles_mask(void) {
    // Hidden singles in rows
    mbit gs;
    bool hit = false;
    gs = find_single(remained_row_mask);
    while (gs) {
      hit = true;
      mbit v = gs & -gs;
      int n = bitpos(v) / 9 + 1;
      int r = bitpos(v) % 9;
      mbit vv = cell_mask[n - 1] & unit_mask[r];
      if (!vv) break;
      int pos2 = bitpos(vv);
      if (can_put(pos2, n)) {
        put(pos2, n);
      }
      gs ^= v;
    }
    // Hidden singles in columns
    gs = find_single(remained_column_mask);
    while (gs) {
      mbit v = gs & -gs;
      int n = bitpos(v) / 9 + 1;
      int c = bitpos(v) % 9;
      mbit vv = cell_mask[n - 1] & unit_mask[c + 9];
      if (!vv) break;
      int pos2 = bitpos(vv);
      if (can_put(pos2, n)) {
        put(pos2, n);
      }
      gs ^= v;
    }
    // Hidden singles in boxes
    gs = find_single(remained_box_mask);
    while (gs) {
      mbit v = gs & -gs;
      int n = bitpos(v) / 9 + 1;  //どの数字か
      int bindex = bitpos(v) % 9; //どのボックスか
      mbit vv = cell_mask[n - 1] & unit_mask[bindex + 18];
      if (!vv) break;
      int pos2 = bitpos(vv);
      if (can_put(pos2, n)) {
        put(pos2, n);
      }
      gs ^= v;
    }
    return hit;
  }

  bool hidden_singles(void) {
    bool hit = false;
    static const mbit mzero = mbit(0);
    for (int i = 0; i < 9; i++) {
      if (cell_mask[i] == mzero) continue;
      for (auto m : unit_mask) {
        const mbit p = cell_mask[i] & m;
        if ((popcnt_u128(p) == 1)) {
          put(bitpos(p), i + 1);
          //printf("puts %d on %d\n", i + 1, bitpos(p));
          hit = true;
        }
      }
    }
    return hit;
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

  //行内二択
  unsigned int row_alt(mbit &mtwo, std::string &answer) {
    mbit v = mtwo & (-mtwo);
    int n = bitpos(v) / 9;
    int r = bitpos(v) % 9;
    mbit mpos = cell_mask[n] & unit_mask[r];
    mbit mpos1 = mpos & -mpos;
    mbit mpos2 = mpos ^ mpos1;
    int pos1 = bitpos(mpos1);
    int pos2 = bitpos(mpos2);
    int sum = 0;
    Grid g1 = (*this);
    g1.put(pos1, n + 1);
    Grid g2 = (*this);
    g2.put(pos2, n + 1);
    sum += g1.solve_internal(answer);
    if (sum > 1) return sum;
    sum += g2.solve_internal(answer);
    return sum;
  }

  //ボックス内二択
  unsigned int box_alt(mbit &mtwo, std::string &answer) {
    mbit v = mtwo & (-mtwo);
    int n = bitpos(v) / 9;
    int b = bitpos(v) % 9;
    mbit mpos = cell_mask[n] & unit_mask[b + 18];
    mbit mpos1 = mpos & -mpos;
    mbit mpos2 = mpos ^ mpos1;
    int pos1 = bitpos(mpos1);
    int pos2 = bitpos(mpos2);
    int sum = 0;
    Grid g1 = (*this);
    g1.put(pos1, n + 1);
    Grid g2 = (*this);
    g2.put(pos2, n + 1);
    sum += g1.solve_internal(answer);
    if (sum > 1) return sum;
    sum += g2.solve_internal(answer);
    return sum;
  }

  // セル内二択による再帰
  unsigned int cell_alt(mbit &mtwo, std::string &answer) {
    mtwo = mtwo & (-mtwo);
    int pos = bitpos(mtwo);
    int n1 = -1;
    int n2 = -1;
    for (int i = 0; i < 9; i++) {
      if (mtwo & cell_mask[i]) {
        if (n1 == -1) {
          n1 = i + 1;
        } else {
          n2 = i + 1;
        }
      }
    }

    int sum = 0;
    Grid g1 = (*this);
    g1.put(pos, n1);
    Grid g2 = (*this);
    g2.put(pos, n2);
    sum += g1.solve_internal(answer);
    if (sum > 1) return sum;
    sum += g2.solve_internal(answer);
    return sum;
  }

  template <class T>
  static void show(T &tt) {
    for (auto m : tt) {
      std::cout << m << std::endl;
    }
    std::cout << std::endl;
  }

  void show() {
    for (int i = 0; i < 81; i++) {
      std::cout << data[i];
    }
    std::cout << std::endl;
  }

//なぜかヘッダで宣言すると遅くなる。要検証。
#if 1
  unsigned int solve_internal(std::string &answer);
#else
  // 解の数を返す
  // 0: 解なし
  // 1: 唯一解あり
  // 2以上: 複数解あり
  unsigned int solve_internal(std::string &answer) {
    bool hit = true;
    // Naked/Hidden singlesで解けるだけ解く
    while (hit) {
      hit = false;
      if (solved_squares()) hit = true;
      //if (hidden_singles()) hit = true;
      if (hidden_singles_mask()) hit = true;
    }
    if (!is_valid()) {
      return 0;
    }

    if (_rest == 0) {
      // 解けたので解答をセット
      answer.resize(81);
      for (int i = 0; i < 81; i++) {
        answer[i] = '0' + data[i];
      }
      return 1;
    }

    //セル内二択
    mbit mtwo = find_two(cell_mask);
    if (mtwo) {
      return cell_alt(mtwo, answer);
    } else {
      return solve_unit(answer);
    }
  }
#endif
};
