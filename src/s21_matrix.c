#include "s21_matrix.h"

int s21_create_matrix(const int rows, const int columns, matrix_t *result) {
  int ret_code = INCORRECT_MATRIX;

  if (rows > 0 && columns > 0) {
    result->columns = columns;
    result->rows = rows;
    result->matrix = calloc(rows, sizeof(double *));
    ret_code = MALLOC_FAILED;

    if (result->matrix) {
      for (int i = 0; i < rows; i++) {
        result->matrix[i] = calloc(columns, sizeof(double));
        if (!result->matrix[i]) break;
      }
      ret_code = OK;
    }
  }

  return (ret_code);
}

void s21_remove_matrix(matrix_t *const A) {
  if (A) {
    for (int i = 0; i < A->rows; i++) free(A->matrix[i]);

    free(A->matrix);
    A->matrix = NULL;
    A->rows = 0;
    A->columns = 0;
  }
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A->rows != B->rows || A->columns != B->columns) return CALC_ERROR;

  int ret_code = s21_create_matrix(A->rows, A->columns, result);

  if (ret_code == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
    }
  }
  return ret_code;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A->rows != B->rows || A->columns != B->columns) return CALC_ERROR;

  int ret_code = s21_create_matrix(A->rows, A->columns, result);

  if (ret_code == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  }
  return ret_code;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int ret_code = s21_create_matrix(A->rows, A->columns, result);
  if (ret_code == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }
  return ret_code;
}

int s21_check_matrix(matrix_t *check) {
  int ret_code = OK;
  if (!check || check->rows < 1 || check->columns < 1) {
    ret_code = INCORRECT_MATRIX;
  }

  return ret_code;
}

// int s21_check_matrix(matrix_t *A) {
//     int res;

//     if (A->matrix == NULL) {
//         res = ERROR;
//     } else if (A->rows == A->columns && A->rows >= 1) {
//         res = OK;
//     } else if ((A->rows != A->columns) && (A->rows && A->columns)) {
//         res = ERROR_RESULT;
//     } else {
//         res = ERROR;
//     }
//     return res;
// }

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (s21_check_matrix(A) || s21_check_matrix(B)) return INCORRECT_MATRIX;

  if (A->columns != B->rows) return CALC_ERROR;

  int ret_code = s21_create_matrix(A->rows, B->columns, result);

  if (ret_code == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns;
           j++) {  // |3 5|       |8 2 3|     |3*8 + 5*1  3*2 + 5*7  3*3 + 5*2|
        result->matrix[i][j] = 0;  // RESULT = A * B =   |2 1|   *   |1 7 2|  =
                                   // |2*8 + 1*1  2*2 + 1*7  2*3 + 1*2|
        for (int m = 0; m < B->rows; m++) {
          result->matrix[i][j] += A->matrix[i][m] * B->matrix[m][j];
        }
      }
    }
  }

  return ret_code;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  if (s21_check_matrix(A)) return INCORRECT_MATRIX;

  int ret_code = s21_create_matrix(A->columns, A->rows, result);

  if (ret_code == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }
  return ret_code;
}

void init_minor(int i, int j, matrix_t *A, matrix_t *minor) {
  int minor_i = 0, minor_j = 0;

  for (int this_i = 0; this_i < A->rows; ++this_i) {
    for (int this_j = 0; this_j < A->columns; ++this_j) {
      if (this_i != i && this_j != j) {
        minor->matrix[minor_i][minor_j] = A->matrix[this_i][this_j];
        ++minor_j;
      }
    }
    minor_j = 0;
    if (this_i != i) ++minor_i;
  }
}

double calculate_determinant(matrix_t *A) {
  double det = 0.0;

  if (A->rows == 2) {
    det = A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
  } else {
    for (int j = 0; j < A->columns; ++j) {
      matrix_t minor = {0};
      s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
      init_minor(0, j, A, &minor);
      if (j % 2) {
        det -= A->matrix[0][j] * calculate_determinant(&minor);
      } else {
        det += A->matrix[0][j] * calculate_determinant(&minor);
      }
      s21_remove_matrix(&minor);
    }
  }
  return det;
}

double get_minor(int i, int j, matrix_t *A) {
  double det_minor = 0.0;
  matrix_t minor = {0};

  s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
  if (A->rows == 2) {
    init_minor(i, j, A, &minor);
    det_minor = minor.matrix[0][0];
  } else {
    init_minor(i, j, A, &minor);
    s21_determinant(&minor, &det_minor);
  }
  s21_remove_matrix(&minor);
  return det_minor;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int ret_code;

  if ((ret_code = s21_check_matrix(A)) == OK) {
    ret_code = s21_create_matrix(A->columns, A->rows, result);
    if (A->rows == 1 && A->columns == 1) {
      result->matrix[0][0] = A->matrix[0][0];
    } else {
      for (int i = 0; i < A->rows; ++i) {
        for (int j = 0; j < A->columns; ++j) {
          if ((i + j) % 2) {
            result->matrix[i][j] = -(get_minor(i, j, A));
          } else {
            result->matrix[i][j] = get_minor(i, j, A);
          }
        }
      }
    }
  }
  return ret_code;
}

int s21_determinant(matrix_t *A, double *result) {
  int ret_code;

  if (((ret_code = s21_check_matrix(A)) == OK)) {
    if (A->rows >= 2) {
      *result = calculate_determinant(A);
    } else {
      *result = A->matrix[0][0];
    }
  }
  return ret_code;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int ret_code;
  double det = 0.0;
  matrix_t calc_comp = {0}, trans = {0};

  if (((ret_code = s21_check_matrix(A)) == OK)) {
    s21_determinant(A, &det);
    if (det != 0.0) {
      if (A->rows == 1 && A->columns == 1) {
        s21_create_matrix(1, 1, result);
        result->matrix[0][0] = 1.0 / A->matrix[0][0];
      } else {
        s21_calc_complements(A, &calc_comp);
        s21_transpose(&calc_comp, &trans);
        s21_mult_number(&trans, (1.0 / det), result);
        s21_remove_matrix(&calc_comp);
        s21_remove_matrix(&trans);
      }
    } else {
      ret_code = CALC_ERROR;
    }
  }
  return ret_code;
}

// int main(void) {
//     // double **A[3][5] = {
//     //     {5.0, 5.0, 23.0, 56.0, 4.0},
//     //     {4.0, 3.0, 7.0, 13.0, 65.0},
//     //     {34.0, 1.0, 23.0, 0.0, 1.0}
//     // };

//     matrix_t X = {0};
//     matrix_t B = {0};
//     s21_create_matrix(4, 4, &X);
//     s21_create_matrix(4, 4, &B);
//     X.matrix[0][0] = 2.0;
//     X.matrix[0][1] = 4.0;
//     X.matrix[0][2] = 6.0;
//     X.matrix[0][3] = 8.0;
//     X.matrix[1][0] = 10.0;
//     X.matrix[1][1] = 12.0;
//     X.matrix[1][2] = 14.0;
//     X.matrix[1][3] = 16.0;
//     X.matrix[2][0] = 18.0;
//     X.matrix[2][1] = 20.0;
//     X.matrix[2][2] = 22.0;
//     X.matrix[2][3] = 24.0;
//     X.matrix[3][0] = 26.0;
//     X.matrix[3][1] = 28.0;
//     X.matrix[3][2] = 30.0;
//     X.matrix[3][3] = 32.0;

//     s21_transpose(&X, &B);
//     for (int i = 0; i < B.rows; i++) {
//         for (int j = 0; j < B.columns; j++) {
//             printf("%1.f ", B.matrix[i][j]);
//             if (j == B.rows - 1) {
//                 printf("\n");
//             }
//         }
//     }
// }

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  if (s21_check_matrix(A) != OK || s21_check_matrix(B) != OK ||
      A->rows != B->rows || A->columns != B->columns)
    return FAILURE;

  for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->columns; j++)
      if (fabs(A->matrix[i][j] - B->matrix[i][j]) > EPS) return FAILURE;

  return SUCCESS;
}