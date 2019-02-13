#include "grid.hpp"

Grid::GridInitializer Grid::si; //マスクの初期化
mbit Grid::unit_mask[27];
mbit Grid::kill_cell_mask[81];
const int Grid::box_index[81];
mbit Grid::kill_row_mask[9][81][9];
mbit Grid::kill_column_mask[9][81][9];
mbit Grid::kill_box_mask[9][81][9];

// マスクの初期化
void Grid::init_masks(void) {
  mbit s = (1 << 9) - 1;
  mbit c = 1;
  mbit b = 7 | (mbit(7) << 9) | (mbit(7) << 18);
  for (int i = 0; i < 8; i++) {
    c = c << 9;
    c = c | 1;
  }
  for (int i = 0; i < 9; i++) {
    unit_mask[i] = s;
    unit_mask[i + 9] = c;
    s = s << 9;
    c = c << 1;
  }
  for (int i = 0; i < 9; i++) {
    unit_mask[i + 18] = b << (((i % 3) * 3) + (i / 3) * 27);
  }
  // kill_cell_maskの作成
  for (int i = 0; i < 81; i++) {
    const int ri = i / 9;
    const int ci = i % 9;
    const int bi = (ci / 3) + (ri / 3) * 3;
    const mbit m = unit_mask[ri] | unit_mask[ci + 9] | unit_mask[bi + 18];
    kill_cell_mask[i] = ((mbit(1) << 81) - 1) & (~m);
  }
  // ユニットごとのkillマスク作成
  for (int n = 0; n < 9; n++) {
    for (int pos = 0; pos < 81; pos++) {
      Grid g;
      g.put(pos, n + 1);
      mbit m_row[9] = {};
      mbit m_column[9] = {};
      mbit m_box[9] = {};
      g.get_kill_mask(m_row, m_column, m_box);
      for (int j = 0; j < 9; j++) {
        kill_row_mask[n][pos][j] = m_row[j];
        kill_column_mask[n][pos][j] = m_column[j];
        kill_box_mask[n][pos][j] = m_box[j];
      }
    }
  }
  int n = 0;
  for (int i = 0; i < 9; i++) {
    std::cout << kill_box_mask[n][0][i] << std::endl;
  }
}

bool Grid::solved_squares(void) {
  static stopwatch::timer<> timer("solved_squares");
  timer.start();
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
  timer.stop();
  return flag;
}

void Grid::solve(std::string &str) {
  Grid g(str);
  std::string answer;
  int n = g.solve_internal(answer);
  if (n == 0) {
    std::cout << "no solution" << std::endl;
  } else if (n > 1) {
    std::cout << "multiple solutions" << std::endl;
  } else {
    std::cout << answer << std::endl;
  }
}

unsigned int Grid::solve_unit(std::string &answer) {
  int min = 9;
  int min_index = -1;
  mbit um = 0;
  for (int i = 0; i < 9; i++) {
    const mbit nm = cell_mask[i];
    if (nm == mbit(0))
      continue;
    for (const auto &m : unit_mask) {
      const int n = popcnt_u128(nm & m);
      // assert(n!=1);
      if (n != 0 && n < min) {
        min_index = i;
        um = m;
        if (n <= 2) {
          goto break_loop;
        }
      }
    }
  }
break_loop:
  mbit v = (cell_mask[min_index] & um);
  int sum = 0;
  while (v) {
    const mbit p = (v & -v);
    const int n = bitpos(p);
    Grid g2 = (*this);
    g2.put(n, min_index + 1);
    sum = sum + g2.solve_internal(answer);
    if (sum > 1)
      return sum;
    v ^= p;
  }
  return sum;
}
// 解の数を返す
// 0: 解なし
// 1: 唯一解あり
// 2以上: 複数解あり
unsigned int Grid::solve_internal(std::string &answer) {
  bool flag = true;
  // Naked/Hidden singlesで解けるだけ解く
  while (flag) {
    flag = false;
    if (solved_squares())
      flag = true;
    if (hidden_singles())
      flag = true;
  }
  if (!is_valid())
    return 0;

  if (_rest == 0) {
    // 解けたので解答をセット
    answer.resize(81);
    for (int i = 0; i < 81; i++) {
      answer[i] = '0' + data[i];
    }
    return 1;
  }
  mbit mtwo = find_two();
  if (mtwo == mbit(0)) {
    return solve_unit(answer);
  }
  mtwo = mtwo & (-mtwo);
  int pos = bitpos(mtwo);
  int sum = 0;
  for (int i = 0; i < 9; i++) {
    if (!(mtwo & cell_mask[i]))
      continue;
    Grid g2 = (*this);
    g2.put(pos, i + 1);
    sum = sum + g2.solve_internal(answer);
    if (sum > 1)
      return sum;
  }

  return sum;
}
