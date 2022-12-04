#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "s21_matrix.h"

#define ERROR 1
#define ERROR_RESULT 2

void init_matrix_by_number(double number, matrix_t *A) {
  for (int i = 0; i < A->rows; ++i) {
    for (int j = 0; j < A->columns; ++j) {
      A->matrix[i][j] = number;
      number += 1.0;
    }
  }
}

START_TEST(s21_create_matrix_01) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(1, 1, &A);
  ck_assert_int_eq(res, OK);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_create_matrix_02) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(2, 2, &A);
  ck_assert_int_eq(res, OK);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_create_matrix_03) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(3, 2, &A);
  ck_assert_int_eq(res, OK);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_create_matrix_04) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(2, 4, &A);
  ck_assert_int_eq(res, OK);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_create_matrix_05) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(10, 15, &A);
  ck_assert_int_eq(res, OK);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_create_matrix_06) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(0, 0, &A);
  ck_assert_int_eq(res, ERROR);
}
END_TEST

START_TEST(s21_create_matrix_07) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(1, 0, &A);
  ck_assert_int_eq(res, ERROR);
}
END_TEST

START_TEST(s21_create_matrix_08) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(0, 1, &A);
  ck_assert_int_eq(res, ERROR);
}
END_TEST

START_TEST(s21_remove_matrix_01) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(2, 2, &A);
  s21_remove_matrix(&A);
  if (A.matrix != NULL) {
    res = ERROR;
  }
  ck_assert_int_eq(res, OK);
}
END_TEST

START_TEST(s21_remove_matrix_02) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(1, 1, &A);
  s21_remove_matrix(&A);
  if (A.matrix != NULL) {
    res = ERROR;
  }
  ck_assert_int_eq(res, OK);
}
END_TEST

START_TEST(s21_remove_matrix_03) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(3, 2, &A);
  s21_remove_matrix(&A);
  if (A.matrix != NULL) {
    res = ERROR;
  }
  ck_assert_int_eq(res, OK);
}
END_TEST

START_TEST(s21_remove_matrix_04) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(3, 4, &A);
  s21_remove_matrix(&A);
  if (A.matrix != NULL) {
    res = ERROR;
  }
  ck_assert_int_eq(res, OK);
}
END_TEST

START_TEST(s21_remove_matrix_05) {
  int res = 0;
  matrix_t A = {0};

  res = s21_create_matrix(3, 2, &A);
  s21_remove_matrix(&A);
  if (A.matrix != NULL) {
    res = ERROR;
  }
  ck_assert_int_eq(res, OK);
}
END_TEST

START_TEST(s21_eq_matrix_01) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(3, 2, &A);
  init_matrix_by_number(1.0, &A);

  s21_create_matrix(3, 2, &B);
  init_matrix_by_number(1.0, &B);

  res = s21_eq_matrix(&A, &B);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_eq_matrix_02) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(3, 2, &A);
  init_matrix_by_number(1.0, &A);

  s21_create_matrix(3, 2, &B);
  init_matrix_by_number(2.0, &B);

  res = s21_eq_matrix(&B, &A);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  ck_assert_int_eq(res, FAILURE);
}
END_TEST

START_TEST(s21_eq_matrix_03) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(3, 4, &A);
  init_matrix_by_number(-5.0, &A);

  s21_create_matrix(3, 4, &B);
  init_matrix_by_number(1.0, &B);

  res = s21_eq_matrix(&A, &B);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  ck_assert_int_eq(res, FAILURE);
}
END_TEST

START_TEST(s21_eq_matrix_04) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(3, 2, &A);
  init_matrix_by_number(10.0, &A);

  s21_create_matrix(3, 2, &B);
  init_matrix_by_number(-10.0, &B);

  res = s21_eq_matrix(&A, &B);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  ck_assert_int_eq(res, FAILURE);
}
END_TEST

START_TEST(s21_eq_matrix_05) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(3, 2, &A);
  init_matrix_by_number(1.0, &A);
  s21_create_matrix(3, 3, &B);
  init_matrix_by_number(1.0, &B);
  res = s21_eq_matrix(&A, &B);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  ck_assert_int_eq(res, FAILURE);
}
END_TEST

START_TEST(s21_eq_matrix_06) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(10, 10, &A);
  init_matrix_by_number(1.0, &A);

  s21_create_matrix(10, 10, &B);
  init_matrix_by_number(1.0, &B);

  res = s21_eq_matrix(&A, &B);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_eq_matrix_07) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(0, 0, &A);

  s21_create_matrix(0, 0, &B);

  res = s21_eq_matrix(&A, &B);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  ck_assert_int_eq(res, FAILURE);
}
END_TEST

START_TEST(s21_eq_matrix_08) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 1;

  s21_create_matrix(1, 1, &B);
  B.matrix[0][0] = 1;

  res = s21_eq_matrix(&A, &B);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sum_matrix_01) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  init_matrix_by_number(1.0, &A);

  s21_create_matrix(4, 4, &B);
  init_matrix_by_number(1.0, &B);

  s21_sum_matrix(&A, &B, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 2.0;
  X.matrix[0][1] = 4.0;
  X.matrix[0][2] = 6.0;
  X.matrix[0][3] = 8.0;
  X.matrix[1][0] = 10.0;
  X.matrix[1][1] = 12.0;
  X.matrix[1][2] = 14.0;
  X.matrix[1][3] = 16.0;
  X.matrix[2][0] = 18.0;
  X.matrix[2][1] = 20.0;
  X.matrix[2][2] = 22.0;
  X.matrix[2][3] = 24.0;
  X.matrix[3][0] = 26.0;
  X.matrix[3][1] = 28.0;
  X.matrix[3][2] = 30.0;
  X.matrix[3][3] = 32.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sum_matrix_02) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  init_matrix_by_number(0.0, &A);

  s21_create_matrix(4, 4, &B);
  init_matrix_by_number(2.0, &B);

  s21_sum_matrix(&A, &B, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 2.0;
  X.matrix[0][1] = 4.0;
  X.matrix[0][2] = 6.0;
  X.matrix[0][3] = 8.0;
  X.matrix[1][0] = 10.0;
  X.matrix[1][1] = 12.0;
  X.matrix[1][2] = 14.0;
  X.matrix[1][3] = 16.0;
  X.matrix[2][0] = 18.0;
  X.matrix[2][1] = 20.0;
  X.matrix[2][2] = 22.0;
  X.matrix[2][3] = 24.0;
  X.matrix[3][0] = 26.0;
  X.matrix[3][1] = 28.0;
  X.matrix[3][2] = 30.0;
  X.matrix[3][3] = 32.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sum_matrix_03) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  init_matrix_by_number(-3.0, &A);

  s21_create_matrix(4, 4, &B);
  init_matrix_by_number(5.0, &B);

  s21_sum_matrix(&A, &B, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 2.0;
  X.matrix[0][1] = 4.0;
  X.matrix[0][2] = 6.0;
  X.matrix[0][3] = 8.0;
  X.matrix[1][0] = 10.0;
  X.matrix[1][1] = 12.0;
  X.matrix[1][2] = 14.0;
  X.matrix[1][3] = 16.0;
  X.matrix[2][0] = 18.0;
  X.matrix[2][1] = 20.0;
  X.matrix[2][2] = 22.0;
  X.matrix[2][3] = 24.0;
  X.matrix[3][0] = 26.0;
  X.matrix[3][1] = 28.0;
  X.matrix[3][2] = 30.0;
  X.matrix[3][3] = 32.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sum_matrix_04) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  init_matrix_by_number(5.0, &A);

  s21_create_matrix(4, 4, &B);
  init_matrix_by_number(-3.0, &B);

  s21_sum_matrix(&A, &B, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 2.0;
  X.matrix[0][1] = 4.0;
  X.matrix[0][2] = 6.0;
  X.matrix[0][3] = 8.0;
  X.matrix[1][0] = 10.0;
  X.matrix[1][1] = 12.0;
  X.matrix[1][2] = 14.0;
  X.matrix[1][3] = 16.0;
  X.matrix[2][0] = 18.0;
  X.matrix[2][1] = 20.0;
  X.matrix[2][2] = 22.0;
  X.matrix[2][3] = 24.0;
  X.matrix[3][0] = 26.0;
  X.matrix[3][1] = 28.0;
  X.matrix[3][2] = 30.0;
  X.matrix[3][3] = 32.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sum_matrix_05) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 2, &A);
  init_matrix_by_number(1.0, &A);

  s21_create_matrix(4, 2, &B);
  init_matrix_by_number(1.0, &B);

  s21_sum_matrix(&A, &B, &Z);

  s21_create_matrix(4, 2, &X);

  X.matrix[0][0] = 2.0;
  X.matrix[0][1] = 4.0;
  X.matrix[1][0] = 6.0;
  X.matrix[1][1] = 8.0;
  X.matrix[2][0] = 10.0;
  X.matrix[2][1] = 12.0;
  X.matrix[3][0] = 14.0;
  X.matrix[3][1] = 16.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sum_matrix_06) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 2, &A);
  init_matrix_by_number(-1.0, &A);

  s21_create_matrix(4, 2, &B);
  init_matrix_by_number(3.0, &B);

  s21_sum_matrix(&A, &B, &Z);

  s21_create_matrix(4, 2, &X);

  X.matrix[0][0] = 2.0;
  X.matrix[0][1] = 4.0;
  X.matrix[1][0] = 6.0;
  X.matrix[1][1] = 8.0;
  X.matrix[2][0] = 10.0;
  X.matrix[2][1] = 12.0;
  X.matrix[3][0] = 14.0;
  X.matrix[3][1] = 16.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sum_matrix_07) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 2, &A);
  init_matrix_by_number(4.0, &A);

  s21_create_matrix(4, 2, &B);
  init_matrix_by_number(-2.0, &B);

  s21_sum_matrix(&A, &B, &Z);

  s21_create_matrix(4, 2, &X);

  X.matrix[0][0] = 2.0;
  X.matrix[0][1] = 4.0;
  X.matrix[1][0] = 6.0;
  X.matrix[1][1] = 8.0;
  X.matrix[2][0] = 10.0;
  X.matrix[2][1] = 12.0;
  X.matrix[3][0] = 14.0;
  X.matrix[3][1] = 16.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sum_matrix_08) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};

  s21_create_matrix(4, 2, &A);

  s21_create_matrix(4, 3, &B);

  res = s21_sum_matrix(&A, &B, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, ERROR_RESULT);
}
END_TEST

START_TEST(s21_sum_matrix_09) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};

  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 1;

  s21_create_matrix(1, 1, &B);
  B.matrix[0][0] = 1;

  res = s21_sum_matrix(&A, &B, &Z);
  ck_assert_double_eq(Z.matrix[0][0], 2.0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, OK);
}
END_TEST

START_TEST(s21_sub_matrix_01) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  init_matrix_by_number(1.0, &A);

  s21_create_matrix(4, 4, &B);
  init_matrix_by_number(1.0, &B);

  s21_sub_matrix(&A, &B, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 0.0;
  X.matrix[0][1] = 0.0;
  X.matrix[0][2] = 0.0;
  X.matrix[0][3] = 0.0;
  X.matrix[1][0] = 0.0;
  X.matrix[1][1] = 0.0;
  X.matrix[1][2] = 0.0;
  X.matrix[1][3] = 0.0;
  X.matrix[2][0] = 0.0;
  X.matrix[2][1] = 0.0;
  X.matrix[2][2] = 0.0;
  X.matrix[2][3] = 0.0;
  X.matrix[3][0] = 0.0;
  X.matrix[3][1] = 0.0;
  X.matrix[3][2] = 0.0;
  X.matrix[3][3] = 0.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sub_matrix_02) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  init_matrix_by_number(0.0, &A);

  s21_create_matrix(4, 4, &B);
  init_matrix_by_number(2.0, &B);

  s21_sub_matrix(&A, &B, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = -2.0;
  X.matrix[0][1] = -2.0;
  X.matrix[0][2] = -2.0;
  X.matrix[0][3] = -2.0;
  X.matrix[1][0] = -2.0;
  X.matrix[1][1] = -2.0;
  X.matrix[1][2] = -2.0;
  X.matrix[1][3] = -2.0;
  X.matrix[2][0] = -2.0;
  X.matrix[2][1] = -2.0;
  X.matrix[2][2] = -2.0;
  X.matrix[2][3] = -2.0;
  X.matrix[3][0] = -2.0;
  X.matrix[3][1] = -2.0;
  X.matrix[3][2] = -2.0;
  X.matrix[3][3] = -2.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sub_matrix_03) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  init_matrix_by_number(-3.0, &A);

  s21_create_matrix(4, 4, &B);
  init_matrix_by_number(5.0, &B);

  s21_sub_matrix(&A, &B, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = -8.0;
  X.matrix[0][1] = -8.0;
  X.matrix[0][2] = -8.0;
  X.matrix[0][3] = -8.0;
  X.matrix[1][0] = -8.0;
  X.matrix[1][1] = -8.0;
  X.matrix[1][2] = -8.0;
  X.matrix[1][3] = -8.0;
  X.matrix[2][0] = -8.0;
  X.matrix[2][1] = -8.0;
  X.matrix[2][2] = -8.0;
  X.matrix[2][3] = -8.0;
  X.matrix[3][0] = -8.0;
  X.matrix[3][1] = -8.0;
  X.matrix[3][2] = -8.0;
  X.matrix[3][3] = -8.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sub_matrix_04) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  init_matrix_by_number(5.0, &A);

  s21_create_matrix(4, 4, &B);
  init_matrix_by_number(-3.0, &B);

  s21_sub_matrix(&A, &B, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 8.0;
  X.matrix[0][1] = 8.0;
  X.matrix[0][2] = 8.0;
  X.matrix[0][3] = 8.0;
  X.matrix[1][0] = 8.0;
  X.matrix[1][1] = 8.0;
  X.matrix[1][2] = 8.0;
  X.matrix[1][3] = 8.0;
  X.matrix[2][0] = 8.0;
  X.matrix[2][1] = 8.0;
  X.matrix[2][2] = 8.0;
  X.matrix[2][3] = 8.0;
  X.matrix[3][0] = 8.0;
  X.matrix[3][1] = 8.0;
  X.matrix[3][2] = 8.0;
  X.matrix[3][3] = 8.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sub_matrix_05) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 2, &A);
  init_matrix_by_number(1.0, &A);

  s21_create_matrix(4, 2, &B);
  init_matrix_by_number(1.0, &B);

  s21_sub_matrix(&A, &B, &Z);

  s21_create_matrix(4, 2, &X);

  X.matrix[0][0] = 0.0;
  X.matrix[0][1] = 0.0;
  X.matrix[1][0] = 0.0;
  X.matrix[1][1] = 0.0;
  X.matrix[2][0] = 0.0;
  X.matrix[2][1] = 0.0;
  X.matrix[3][0] = 0.0;
  X.matrix[3][1] = 0.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sub_matrix_06) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 2, &A);
  init_matrix_by_number(-1.0, &A);

  s21_create_matrix(4, 2, &B);
  init_matrix_by_number(3.0, &B);

  s21_sub_matrix(&A, &B, &Z);

  s21_create_matrix(4, 2, &X);

  X.matrix[0][0] = -4.0;
  X.matrix[0][1] = -4.0;
  X.matrix[1][0] = -4.0;
  X.matrix[1][1] = -4.0;
  X.matrix[2][0] = -4.0;
  X.matrix[2][1] = -4.0;
  X.matrix[3][0] = -4.0;
  X.matrix[3][1] = -4.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sub_matrix_07) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 2, &A);
  init_matrix_by_number(4.0, &A);

  s21_create_matrix(4, 2, &B);
  init_matrix_by_number(-2.0, &B);

  s21_sub_matrix(&A, &B, &Z);

  s21_create_matrix(4, 2, &X);

  X.matrix[0][0] = 6.0;
  X.matrix[0][1] = 6.0;
  X.matrix[1][0] = 6.0;
  X.matrix[1][1] = 6.0;
  X.matrix[2][0] = 6.0;
  X.matrix[2][1] = 6.0;
  X.matrix[3][0] = 6.0;
  X.matrix[3][1] = 6.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_sub_matrix_08) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};

  s21_create_matrix(4, 2, &A);

  s21_create_matrix(4, 3, &B);

  res = s21_sub_matrix(&A, &B, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, ERROR_RESULT);
}
END_TEST

START_TEST(s21_sub_matrix_09) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};

  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 1;

  s21_create_matrix(1, 1, &B);
  B.matrix[0][0] = 1;

  res = s21_sub_matrix(&A, &B, &Z);
  ck_assert_double_eq(Z.matrix[0][0], 0.0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, OK);
}
END_TEST

START_TEST(s21_mult_number_01) {
  int res = 0;
  double number = 3.0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  init_matrix_by_number(1.0, &A);

  s21_mult_number(&A, number, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 3.0;
  X.matrix[0][1] = 6.0;
  X.matrix[0][2] = 9.0;
  X.matrix[0][3] = 12.0;
  X.matrix[1][0] = 15.0;
  X.matrix[1][1] = 18.0;
  X.matrix[1][2] = 21.0;
  X.matrix[1][3] = 24.0;
  X.matrix[2][0] = 27.0;
  X.matrix[2][1] = 30.0;
  X.matrix[2][2] = 33.0;
  X.matrix[2][3] = 36.0;
  X.matrix[3][0] = 39.0;
  X.matrix[3][1] = 42.0;
  X.matrix[3][2] = 45.0;
  X.matrix[3][3] = 48.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_number_02) {
  int res = 0;
  double number = 5.0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  init_matrix_by_number(1.0, &A);

  s21_mult_number(&A, number, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 5.0;
  X.matrix[0][1] = 10.0;
  X.matrix[0][2] = 15.0;
  X.matrix[0][3] = 20.0;
  X.matrix[1][0] = 25.0;
  X.matrix[1][1] = 30.0;
  X.matrix[1][2] = 35.0;
  X.matrix[1][3] = 40.0;
  X.matrix[2][0] = 45.0;
  X.matrix[2][1] = 50.0;
  X.matrix[2][2] = 55.0;
  X.matrix[2][3] = 60.0;
  X.matrix[3][0] = 65.0;
  X.matrix[3][1] = 70.0;
  X.matrix[3][2] = 75.0;
  X.matrix[3][3] = 80.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_number_03) {
  int res = 0;
  double number = 2.0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 2, &A);
  init_matrix_by_number(2.0, &A);

  s21_mult_number(&A, number, &Z);

  s21_create_matrix(4, 2, &X);

  X.matrix[0][0] = 4.0;
  X.matrix[0][1] = 6.0;
  X.matrix[1][0] = 8.0;
  X.matrix[1][1] = 10.0;
  X.matrix[2][0] = 12.0;
  X.matrix[2][1] = 14.0;
  X.matrix[3][0] = 16.0;
  X.matrix[3][1] = 18.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_number_04) {
  int res = 0;
  double number = -2.0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 2, &A);
  init_matrix_by_number(2.0, &A);

  s21_mult_number(&A, number, &Z);

  s21_create_matrix(4, 2, &X);

  X.matrix[0][0] = -4.0;
  X.matrix[0][1] = -6.0;
  X.matrix[1][0] = -8.0;
  X.matrix[1][1] = -10.0;
  X.matrix[2][0] = -12.0;
  X.matrix[2][1] = -14.0;
  X.matrix[3][0] = -16.0;
  X.matrix[3][1] = -18.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_number_05) {
  int res = 0;
  double number = 2.0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(1, 1, &A);
  init_matrix_by_number(2.0, &A);

  s21_mult_number(&A, number, &Z);

  s21_create_matrix(1, 1, &X);

  X.matrix[0][0] = 4.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_number_06) {
  int res = 0;
  double number = 2.0;
  matrix_t A = {0};
  matrix_t Z = {0};

  s21_create_matrix(0, 0, &A);

  res = s21_mult_number(&A, number, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, ERROR);
}
END_TEST

START_TEST(s21_mult_matrix_01) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  init_matrix_by_number(1.0, &A);

  s21_create_matrix(4, 4, &B);
  init_matrix_by_number(1.0, &B);

  s21_mult_matrix(&A, &B, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 90.0;
  X.matrix[0][1] = 100.0;
  X.matrix[0][2] = 110.0;
  X.matrix[0][3] = 120.0;
  X.matrix[1][0] = 202.0;
  X.matrix[1][1] = 228.0;
  X.matrix[1][2] = 254.0;
  X.matrix[1][3] = 280.0;
  X.matrix[2][0] = 314.0;
  X.matrix[2][1] = 356.0;
  X.matrix[2][2] = 398.0;
  X.matrix[2][3] = 440.0;
  X.matrix[3][0] = 426.0;
  X.matrix[3][1] = 484.0;
  X.matrix[3][2] = 542.0;
  X.matrix[3][3] = 600.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_matrix_02) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  init_matrix_by_number(1.0, &A);

  s21_create_matrix(4, 4, &B);
  init_matrix_by_number(-5.0, &B);

  s21_mult_matrix(&A, &B, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 30.0;
  X.matrix[0][1] = 40.0;
  X.matrix[0][2] = 50.0;
  X.matrix[0][3] = 60.0;
  X.matrix[1][0] = 46.0;
  X.matrix[1][1] = 72.0;
  X.matrix[1][2] = 98.0;
  X.matrix[1][3] = 124.0;
  X.matrix[2][0] = 62.0;
  X.matrix[2][1] = 104.0;
  X.matrix[2][2] = 146.0;
  X.matrix[2][3] = 188.0;
  X.matrix[3][0] = 78.0;
  X.matrix[3][1] = 136.0;
  X.matrix[3][2] = 194.0;
  X.matrix[3][3] = 252.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_matrix_03) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  init_matrix_by_number(-5.0, &A);

  s21_create_matrix(4, 4, &B);
  init_matrix_by_number(1.0, &B);

  s21_mult_matrix(&A, &B, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = -78.0;
  X.matrix[0][1] = -92.0;
  X.matrix[0][2] = -106.0;
  X.matrix[0][3] = -120.0;
  X.matrix[1][0] = 34.0;
  X.matrix[1][1] = 36.0;
  X.matrix[1][2] = 38.0;
  X.matrix[1][3] = 40.0;
  X.matrix[2][0] = 146.0;
  X.matrix[2][1] = 164.0;
  X.matrix[2][2] = 182.0;
  X.matrix[2][3] = 200.0;
  X.matrix[3][0] = 258.0;
  X.matrix[3][1] = 292.0;
  X.matrix[3][2] = 326.0;
  X.matrix[3][3] = 360.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_matrix_04) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 2, &A);
  init_matrix_by_number(-5.0, &A);

  s21_create_matrix(2, 4, &B);
  init_matrix_by_number(1.0, &B);

  s21_mult_matrix(&A, &B, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = -25.0;
  X.matrix[0][1] = -34.0;
  X.matrix[0][2] = -43.0;
  X.matrix[0][3] = -52.0;
  X.matrix[1][0] = -13.0;
  X.matrix[1][1] = -18.0;
  X.matrix[1][2] = -23.0;
  X.matrix[1][3] = -28.0;
  X.matrix[2][0] = -1.0;
  X.matrix[2][1] = -2.0;
  X.matrix[2][2] = -3.0;
  X.matrix[2][3] = -4.0;
  X.matrix[3][0] = 11.0;
  X.matrix[3][1] = 14.0;
  X.matrix[3][2] = 17.0;
  X.matrix[3][3] = 20.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_matrix_05) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(1, 1, &A);
  init_matrix_by_number(-5.0, &A);

  s21_create_matrix(1, 1, &B);
  init_matrix_by_number(1.0, &B);

  s21_mult_matrix(&A, &B, &Z);

  s21_create_matrix(1, 1, &X);

  X.matrix[0][0] = -5.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_mult_matrix_06) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};

  res = s21_mult_matrix(&A, &B, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, ERROR);
}
END_TEST

START_TEST(s21_mult_matrix_07) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t Z = {0};

  s21_create_matrix(2, 3, &A);
  init_matrix_by_number(-5.0, &A);

  s21_create_matrix(2, 2, &B);
  init_matrix_by_number(1.0, &B);

  res = s21_mult_matrix(&A, &B, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, ERROR_RESULT);
}
END_TEST

START_TEST(s21_transpose_01) {
  int res = 0;
  matrix_t A = {0};

  res = s21_transpose(&A, &A);
  ck_assert_int_eq(res, ERROR);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_transpose_02) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(3, 4, &A);
  init_matrix_by_number(1.0, &A);

  res = s21_transpose(&A, &B);
  res = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(res, FAILURE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_transpose_03) {
  int res = 0;
  matrix_t A = {0};
  matrix_t T = {0};
  matrix_t B = {0};

  s21_create_matrix(2, 3, &A);
  init_matrix_by_number(1.0, &A);
  s21_create_matrix(3, 2, &B);
  B.matrix[0][0] = 1.0;
  B.matrix[0][1] = 4.0;
  B.matrix[1][0] = 2.0;
  B.matrix[1][1] = 5.0;
  B.matrix[2][0] = 3.0;
  B.matrix[2][1] = 6.0;

  res = s21_transpose(&A, &T);
  res = s21_eq_matrix(&T, &B);
  ck_assert_int_eq(res, SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&T);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_transpose_04) {
  int res = 0;
  matrix_t A = {0};
  matrix_t B = {0};

  s21_create_matrix(1, 1, &A);
  init_matrix_by_number(1.0, &A);

  res = s21_transpose(&A, &B);
  res = s21_eq_matrix(&A, &B);

  ck_assert_int_eq(res, SUCCESS);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_determinant_01) {
  int res = 0;
  double determ = 0.0;
  matrix_t A = {0};

  res = s21_determinant(&A, &determ);
  ck_assert_int_eq(res, ERROR);
  s21_remove_matrix(&A);
}
END_TEST

// START_TEST(s21_determinant_02) {
//     int res = 0;
//     double determ = 0.0;
//     matrix_t A = {0};
//     s21_create_matrix(1, 2, &A);

//     res = s21_determinant(&A, &determ);
//     ck_assert_int_eq(res, ERROR_RESULT);
//     s21_remove_matrix(&A);
// }
// END_TEST

START_TEST(s21_determinant_03) {
  double determ = 0.0;
  matrix_t A = {0};

  s21_create_matrix(1, 1, &A);
  init_matrix_by_number(1.0, &A);

  s21_determinant(&A, &determ);
  ck_assert_double_eq(determ, 1.0);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_04) {
  double determ = 0.0;
  matrix_t A = {0};

  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 4000.87;
  A.matrix[0][1] = 454.0;
  A.matrix[0][2] = 4.0;
  A.matrix[1][0] = 2.0;
  A.matrix[1][1] = 1.897;
  A.matrix[1][2] = 3.0;
  A.matrix[2][0] = 0.000087;
  A.matrix[2][1] = 2.5668;
  A.matrix[2][2] = 1.0;

  s21_determinant(&A, &determ);
  ck_assert_double_eq(determ, -24105.996724156);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_05) {
  double determ = 0.0;
  matrix_t A = {0};

  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = -677700.0;
  A.matrix[0][1] = 654.0;
  A.matrix[1][0] = 2.0;
  A.matrix[1][1] = -0.00001;

  s21_determinant(&A, &determ);
  ck_assert_double_eq(determ, -1301.223);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_06) {
  double determ = 0.0;
  matrix_t A = {0};

  s21_create_matrix(8, 8, &A);
  init_matrix_by_number(1.0, &A);

  s21_determinant(&A, &determ);
  ck_assert_double_eq(determ, 0.0);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_07) {
  double determ = 0.0;
  matrix_t A = {0};

  s21_create_matrix(1, 1, &A);
  init_matrix_by_number(1.0, &A);

  s21_determinant(&A, &determ);
  ck_assert_double_eq(determ, 1.0);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_calc_complements_01) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  A.matrix[0][0] = 112.0;
  A.matrix[0][1] = 243.0;
  A.matrix[0][2] = 3.0;
  A.matrix[0][3] = -8.0;
  A.matrix[1][0] = 47.0;
  A.matrix[1][1] = 51.0;
  A.matrix[1][2] = -66.0;
  A.matrix[1][3] = 99.0;
  A.matrix[2][0] = -74.0;
  A.matrix[2][1] = 85.0;
  A.matrix[2][2] = 97.0;
  A.matrix[2][3] = 63.0;
  A.matrix[3][0] = -13.0;
  A.matrix[3][1] = 79.0;
  A.matrix[3][2] = -99.0;
  A.matrix[3][3] = -121.0;

  s21_calc_complements(&A, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = -2879514.0;
  X.matrix[0][1] = -1236631.0;
  X.matrix[0][2] = -1685096.0;
  X.matrix[0][3] = 880697.0;
  X.matrix[1][0] = 1162090.0;
  X.matrix[1][1] = -714015.0;
  X.matrix[1][2] = 4046255.0;
  X.matrix[1][3] = -3901600.0;
  X.matrix[2][0] = 4362897.0;
  X.matrix[2][1] = -2049432.0;
  X.matrix[2][2] = -532912.0;
  X.matrix[2][3] = -1370781.0;
  X.matrix[3][0] = 3412773.0;
  X.matrix[3][1] = -1569493.0;
  X.matrix[3][2] = 3144517.0;
  X.matrix[3][3] = 1284666.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_calc_complements_02) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  A.matrix[0][0] = -1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 7.0;
  A.matrix[0][3] = 9.0;
  A.matrix[1][0] = 1.0;
  A.matrix[1][1] = 0.0;
  A.matrix[1][2] = 0.0;
  A.matrix[1][3] = 0.0;
  A.matrix[2][0] = 47.0;
  A.matrix[2][1] = 13.0;
  A.matrix[2][2] = 17.0;
  A.matrix[2][3] = 21.0;
  A.matrix[3][0] = 22.0;
  A.matrix[3][1] = 7.0;
  A.matrix[3][2] = 1.0;
  A.matrix[3][3] = 3.0;

  s21_calc_complements(&A, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 0.0;
  X.matrix[0][1] = -30.0;
  X.matrix[0][2] = -108.0;
  X.matrix[0][3] = 106.0;
  X.matrix[1][0] = 138.0;
  X.matrix[1][1] = -726.0;
  X.matrix[1][2] = -1137.0;
  X.matrix[1][3] = 1061.0;
  X.matrix[2][0] = 0.0;
  X.matrix[2][1] = 12.0;
  X.matrix[2][2] = 57.0;
  X.matrix[2][3] = -47.0;
  X.matrix[3][0] = 0.0;
  X.matrix[3][1] = 6.0;
  X.matrix[3][2] = -75.0;
  X.matrix[3][3] = 57.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

// START_TEST(s21_calc_complements_03) {
//     int res = 0;
//     matrix_t A = {0};
//     matrix_t Z = {0};

//     s21_create_matrix(4, 3, &A);
//     A.matrix[0][0] = -1.0;
//     A.matrix[0][1] = 2.0;
//     A.matrix[0][2] = 7.0;
//     A.matrix[1][0] = 1.0;
//     A.matrix[1][1] = 0.0;
//     A.matrix[1][2] = 0.0;
//     A.matrix[2][0] = 47.0;
//     A.matrix[2][1] = 13.0;
//     A.matrix[2][2] = 17.0;
//     A.matrix[3][0] = 22.0;
//     A.matrix[3][1] = 7.0;
//     A.matrix[3][2] = 1.0;

//     res = s21_calc_complements(&A, &Z);

//     s21_remove_matrix(&A);

//     ck_assert_int_eq(res, ERROR_RESULT);
// }
// END_TEST

START_TEST(s21_calc_complements_04) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};

  s21_create_matrix(2, 2, &A);
  s21_remove_matrix(&A);

  res = s21_calc_complements(&A, &Z);

  ck_assert_int_eq(res, ERROR);
}
END_TEST

START_TEST(s21_inverse_matrix_01) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};
  matrix_t X = {0};

  s21_create_matrix(4, 4, &A);
  A.matrix[0][0] = -1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 7.0;
  A.matrix[0][3] = 9.0;
  A.matrix[1][0] = 1.0;
  A.matrix[1][1] = 0.0;
  A.matrix[1][2] = 0.0;
  A.matrix[1][3] = 0.0;
  A.matrix[2][0] = 47.0;
  A.matrix[2][1] = 13.0;
  A.matrix[2][2] = 17.0;
  A.matrix[2][3] = 21.0;
  A.matrix[3][0] = 22.0;
  A.matrix[3][1] = 7.0;
  A.matrix[3][2] = 1.0;
  A.matrix[3][3] = 3.0;

  s21_inverse_matrix(&A, &Z);

  s21_create_matrix(4, 4, &X);

  X.matrix[0][0] = 0.0;
  X.matrix[0][1] = 1.0;
  X.matrix[0][2] = 0.0;
  X.matrix[0][3] = 0.0;
  X.matrix[1][0] = -5.0 / 23.0;
  X.matrix[1][1] = -121.0 / 23.0;
  X.matrix[1][2] = 2.0 / 23.0;
  X.matrix[1][3] = 1.0 / 23.0;
  X.matrix[2][0] = -18.0 / 23.0;
  X.matrix[2][1] = -379.0 / 46.0;
  X.matrix[2][2] = 19.0 / 46.0;
  X.matrix[2][3] = -25.0 / 46.0;
  X.matrix[3][0] = 53.0 / 69.0;
  X.matrix[3][1] = 1061.0 / 138.0;
  X.matrix[3][2] = -47.0 / 138.0;
  X.matrix[3][3] = 19.0 / 46.0;

  res = s21_eq_matrix(&X, &Z);

  s21_remove_matrix(&A);
  s21_remove_matrix(&X);
  s21_remove_matrix(&Z);

  ck_assert_int_eq(res, SUCCESS);
}
END_TEST

START_TEST(s21_inverse_matrix_02) {
  int res;
  matrix_t A = {0};
  matrix_t Z = {0};

  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 2;
  res = s21_inverse_matrix(&A, &Z);
  ck_assert_int_eq((Z.matrix[0][0] == (1.0 / A.matrix[0][0])), 1);
  ck_assert_int_eq(res, OK);
  s21_remove_matrix(&A);
  s21_remove_matrix(&Z);
}
END_TEST

START_TEST(s21_inverse_matrix_03) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};

  s21_create_matrix(3, 3, &A);
  init_matrix_by_number(1.0, &A);

  res = s21_inverse_matrix(&A, &Z);
  ck_assert_int_eq(res, ERROR_RESULT);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_inverse_matrix_04) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};

  s21_create_matrix(3, 2, &A);
  init_matrix_by_number(1.0, &A);

  res = s21_inverse_matrix(&A, &Z);
  ck_assert_int_eq(res, ERROR_RESULT);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_inverse_matrix_05) {
  int res = 0;
  matrix_t A = {0};
  matrix_t Z = {0};

  s21_create_matrix(3, 3, &A);
  init_matrix_by_number(1.0, &A);
  s21_remove_matrix(&A);

  res = s21_inverse_matrix(&A, &Z);
  ck_assert_int_eq(res, ERROR);
}
END_TEST

Suite *s21_create_matrix_test(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("s21_create_matrix_test");
  tc_core = tcase_create("s21_create_matrix_test");

  tcase_add_test(tc_core, s21_create_matrix_01);
  tcase_add_test(tc_core, s21_create_matrix_02);
  tcase_add_test(tc_core, s21_create_matrix_03);
  tcase_add_test(tc_core, s21_create_matrix_04);
  tcase_add_test(tc_core, s21_create_matrix_05);
  tcase_add_test(tc_core, s21_create_matrix_06);
  tcase_add_test(tc_core, s21_create_matrix_07);
  tcase_add_test(tc_core, s21_create_matrix_08);

  suite_add_tcase(s, tc_core);

  return s;
}

Suite *s21_remove_matrix_test(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("s21_remove_matrix_test");
  tc_core = tcase_create("s21_remove_matrix_test");

  tcase_add_test(tc_core, s21_remove_matrix_01);
  tcase_add_test(tc_core, s21_remove_matrix_02);
  tcase_add_test(tc_core, s21_remove_matrix_03);
  tcase_add_test(tc_core, s21_remove_matrix_04);
  tcase_add_test(tc_core, s21_remove_matrix_05);

  suite_add_tcase(s, tc_core);

  return s;
}

Suite *s21_eq_matrix_test(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("s21_eq_matrix_test");
  tc_core = tcase_create("s21_eq_matrix_test");

  tcase_add_test(tc_core, s21_eq_matrix_01);
  tcase_add_test(tc_core, s21_eq_matrix_02);
  tcase_add_test(tc_core, s21_eq_matrix_03);
  tcase_add_test(tc_core, s21_eq_matrix_04);
  tcase_add_test(tc_core, s21_eq_matrix_05);
  tcase_add_test(tc_core, s21_eq_matrix_06);
  tcase_add_test(tc_core, s21_eq_matrix_07);
  tcase_add_test(tc_core, s21_eq_matrix_08);

  suite_add_tcase(s, tc_core);

  return s;
}

Suite *s21_sum_matrix_test(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("s21_sum_matrix_test");
  tc_core = tcase_create("s21_sum_matrix_test");

  tcase_add_test(tc_core, s21_sum_matrix_01);
  tcase_add_test(tc_core, s21_sum_matrix_02);
  tcase_add_test(tc_core, s21_sum_matrix_03);
  tcase_add_test(tc_core, s21_sum_matrix_04);
  tcase_add_test(tc_core, s21_sum_matrix_05);
  tcase_add_test(tc_core, s21_sum_matrix_06);
  tcase_add_test(tc_core, s21_sum_matrix_07);
  tcase_add_test(tc_core, s21_sum_matrix_08);
  tcase_add_test(tc_core, s21_sum_matrix_09);

  suite_add_tcase(s, tc_core);

  return s;
}

Suite *s21_sub_matrix_test(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("s21_sub_matrix_test");
  tc_core = tcase_create("s21_sub_matrix_test");

  tcase_add_test(tc_core, s21_sub_matrix_01);
  tcase_add_test(tc_core, s21_sub_matrix_02);
  tcase_add_test(tc_core, s21_sub_matrix_03);
  tcase_add_test(tc_core, s21_sub_matrix_04);
  tcase_add_test(tc_core, s21_sub_matrix_05);
  tcase_add_test(tc_core, s21_sub_matrix_06);
  tcase_add_test(tc_core, s21_sub_matrix_07);
  tcase_add_test(tc_core, s21_sub_matrix_08);
  tcase_add_test(tc_core, s21_sub_matrix_09);

  suite_add_tcase(s, tc_core);

  return s;
}

Suite *s21_mult_number_test(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("s21_mult_number_test");
  tc_core = tcase_create("s21_mult_number_test");

  tcase_add_test(tc_core, s21_mult_number_01);
  tcase_add_test(tc_core, s21_mult_number_02);
  tcase_add_test(tc_core, s21_mult_number_03);
  tcase_add_test(tc_core, s21_mult_number_04);
  tcase_add_test(tc_core, s21_mult_number_05);
  tcase_add_test(tc_core, s21_mult_number_06);

  suite_add_tcase(s, tc_core);

  return s;
}

Suite *s21_mult_matrix_test(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("s21_mult_matrix_test");
  tc_core = tcase_create("s21_mult_matrix_test");

  tcase_add_test(tc_core, s21_mult_matrix_01);
  tcase_add_test(tc_core, s21_mult_matrix_02);
  tcase_add_test(tc_core, s21_mult_matrix_03);
  tcase_add_test(tc_core, s21_mult_matrix_04);
  tcase_add_test(tc_core, s21_mult_matrix_05);
  tcase_add_test(tc_core, s21_mult_matrix_06);
  tcase_add_test(tc_core, s21_mult_matrix_07);

  suite_add_tcase(s, tc_core);

  return s;
}

Suite *s21_transpose_test(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("s21_transpose_test");
  tc_core = tcase_create("s21_transpose_test");

  tcase_add_test(tc_core, s21_transpose_01);
  tcase_add_test(tc_core, s21_transpose_02);
  tcase_add_test(tc_core, s21_transpose_03);
  tcase_add_test(tc_core, s21_transpose_04);

  suite_add_tcase(s, tc_core);

  return s;
}

Suite *s21_determinant_test(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("s21_determinant_test");
  tc_core = tcase_create("s21_determinant_test");

  tcase_add_test(tc_core, s21_determinant_01);
  // tcase_add_test(tc_core, s21_determinant_02);
  tcase_add_test(tc_core, s21_determinant_03);
  tcase_add_test(tc_core, s21_determinant_04);
  tcase_add_test(tc_core, s21_determinant_05);
  tcase_add_test(tc_core, s21_determinant_06);
  tcase_add_test(tc_core, s21_determinant_07);

  suite_add_tcase(s, tc_core);

  return s;
}

Suite *s21_calc_complements_test(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("s21_calc_complements_test");
  tc_core = tcase_create("s21_calc_complements_test");

  tcase_add_test(tc_core, s21_calc_complements_01);
  tcase_add_test(tc_core, s21_calc_complements_02);
  // tcase_add_test(tc_core, s21_calc_complements_03);
  tcase_add_test(tc_core, s21_calc_complements_04);

  suite_add_tcase(s, tc_core);

  return s;
}

Suite *s21_inverse_matrix_test(void) {
  Suite *s;
  TCase *tc_core;

  s = suite_create("s21_inverse_matrix_test");
  tc_core = tcase_create("s21_inverse_matrix_test");

  tcase_add_test(tc_core, s21_inverse_matrix_01);
  tcase_add_test(tc_core, s21_inverse_matrix_02);
  tcase_add_test(tc_core, s21_inverse_matrix_03);
  tcase_add_test(tc_core, s21_inverse_matrix_04);
  tcase_add_test(tc_core, s21_inverse_matrix_05);

  suite_add_tcase(s, tc_core);

  return s;
}

int main(void) {
  int no_failed = 0;
  Suite *s;
  SRunner *runner;

  s = s21_create_matrix_test();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_NORMAL);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  s = s21_remove_matrix_test();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_NORMAL);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  s = s21_eq_matrix_test();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_NORMAL);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  s = s21_sum_matrix_test();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_NORMAL);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  s = s21_sub_matrix_test();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_NORMAL);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  s = s21_mult_number_test();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_NORMAL);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  s = s21_mult_matrix_test();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_NORMAL);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  s = s21_transpose_test();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_NORMAL);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  s = s21_determinant_test();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_NORMAL);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  s = s21_calc_complements_test();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_NORMAL);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  s = s21_inverse_matrix_test();
  runner = srunner_create(s);
  srunner_run_all(runner, CK_NORMAL);
  no_failed += srunner_ntests_failed(runner);
  srunner_free(runner);

  return (no_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
