#include "s21_matrix.h"

int s21_is_equal(matrix_t *A, matrix_t *B) {
  int ret = FAILURE;
  if (A->rows == B->rows && A->columns == B->columns) ret = SUCCESS;
  return ret;
}

int s21_is_correct(matrix_t *m) {
  int ret = FAILURE;
  if (m && m->matrix && m->rows >= 1 && m->columns >= 1) ret = SUCCESS;
  return ret;
}

double s21_determinant_recursive(matrix_t *A) {
  double result = 0;
  if (A->rows == 2) {
    result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
  } else {
    for (int i = 0; i < A->rows; i++) {
      matrix_t minor;
      Minor(1, i + 1, A, &minor);
      result +=
          pow((-1), i) * A->matrix[0][i] * s21_determinant_recursive(&minor);
      s21_remove_matrix(&minor);
    }
  }
  return result;
}

int Minor(int row, int column, matrix_t *A, matrix_t *result) {
  int err = 1;
  if (A->matrix != NULL) {
    err = s21_create_matrix(A->rows - 1, A->columns - 1, result);
    if (err == 0) {
      int m, n;
      for (int i = 0; i < A->rows; i++) {
        m = i;
        if (i > row - 1) {
          m--;
        }
        for (int j = 0; j < A->columns; j++) {
          n = j;
          if (j > column - 1) {
            n--;
          }
          if (i != row - 1 && j != column - 1) {
            result->matrix[m][n] = A->matrix[i][j];
          }
        }
      }
    }
  }
  return err;
}

int calc_helper(matrix_t *A, matrix_t *result) {
  int err = 0;
  result->matrix[0][0] = 1;
  if (A->rows != 1) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        double deted;
        matrix_t minored;
        err = Minor(i + 1, j + 1, A, &minored);
        if (err == 0) {
          err = s21_determinant(&minored, &deted);
          if (err == 0) {
            result->matrix[i][j] = pow((-1), i + j) * deted;
          }
        }
        s21_remove_matrix(&minored);
      }
    }
  }
  return err;
}