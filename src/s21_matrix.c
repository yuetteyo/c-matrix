#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int err = OK;
  if (!result || rows < 1 || columns < 1) {
    err = INCORRECT_MATRIX;
  } else {
    result->rows = rows;  // (*result).rows
    result->columns = columns;
    result->matrix = calloc(rows, sizeof(double *));
    for (int i = 0; i < rows; i++) {
      result->matrix[i] = calloc(columns, sizeof(double));
    }
  }
  return err;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix) {
    for (int i = 0; i < A->rows; i++) {
      if (A->matrix[i]) {
        free(A->matrix[i]);
        A->matrix[i] = NULL;
      }
    }
    free(A->matrix);
    A->matrix = NULL;
  }
  if (A->rows) A->rows = 0;
  if (A->columns) A->columns = 0;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int err = SUCCESS;
  if (A->matrix == NULL || B->matrix == NULL || !s21_is_equal(A, B)) {
    err = FAILURE;
  } else {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) >= EPS) err = FAILURE;
      }
    }
  }
  return err;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int err = OK;
  if (!s21_is_correct(A) || !s21_is_correct(B)) {
    err = INCORRECT_MATRIX;
  } else if (!s21_is_equal(A, B)) {
    err = CALCULATION_ERROR;
  } else if ((err = s21_create_matrix(A->rows, A->columns, result)) == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
    }
  }
  return err;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int err = OK;
  if (!s21_is_correct(A) || !s21_is_correct(B)) {
    err = INCORRECT_MATRIX;
  } else if (!s21_is_equal(A, B)) {
    err = CALCULATION_ERROR;
  } else if ((err = s21_create_matrix(A->rows, A->columns, result)) == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  }
  return err;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int err = OK;
  if (!s21_is_correct(A)) {
    err = INCORRECT_MATRIX;
  } else if ((err = s21_create_matrix(A->rows, A->columns, result)) == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }
  return err;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int err = OK, check = OK;
  if (!s21_is_correct(A) || !s21_is_correct(B)) {
    err = INCORRECT_MATRIX;
  } else if (A->columns != B->rows) {
    err = CALCULATION_ERROR;
  } else {
    check = s21_create_matrix(A->rows, B->columns, result);
    for (int i = 0; check == OK && i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        for (int k = 0; k < A->columns; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  }
  return err;
}

int s21_determinant(matrix_t *A, double *result) {
  int err = OK;
  if (!s21_is_correct(A)) {
    err = INCORRECT_MATRIX;
  } else if (A->rows != A->columns) {
    err = CALCULATION_ERROR;
  } else if (A->rows == 1) {
    *result = A->matrix[0][0];
  } else if (A->rows != 1) {
    *result = s21_determinant_recursive(A);
  }
  return err;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int err = 1;
  if (!s21_is_correct(A))
    err = INCORRECT_MATRIX;
  else if (A->rows != A->columns)
    err = CALCULATION_ERROR;
  else {
    err = s21_create_matrix(A->rows, A->columns, result);
    if (err == 0) err = calc_helper(A, result);
  }
  return err;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int err = OK;
  if (!s21_is_correct(A) || result == NULL) err = INCORRECT_MATRIX;
  if ((err = s21_create_matrix(A->columns, A->rows, result)) == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }
  return err;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int err = 1;
  if (s21_is_correct(A)) {
    err = 2;
    double det;
    s21_determinant(A, &det);
    if (fabs(det - 0) > 1e-6) {
      matrix_t tmp_calc;
      err = s21_calc_complements(A, &tmp_calc);
      if (err == 0) {
        matrix_t tmp_trans;
        err = s21_transpose(&tmp_calc, &tmp_trans);
        if (err == 0) {
          s21_mult_number(&tmp_trans, 1 / det, result);
        }
        s21_remove_matrix(&tmp_trans);
      }
      s21_remove_matrix(&tmp_calc);
    }
  }
  return err;
}
