/*
    Matrix - A simple matrix class.

    This class provides functions for a few basic matrix operations, such as
    sum, product and transpose. There are also functions to compute the Cholesky
    decomposition, its inverse and the square root.
    
    Most functions have two versions: one to create and return a new matrix
    ("new" prefix) and; one to perform the operation in-place on an existing 
    matrix ("this" prefix") in order to improve efficiency.

    Author: Roger Stuckey, 2010.
*/

static class Matrix {
  int ROWS, COLS; // storage size
  int rows, cols; // active size
  float[][] data;

  Matrix(int rows_, int cols_) {
    // Create a new Matrix of zeros
    ROWS = rows_;
    COLS = cols_;
    rows = rows_;
    cols = cols_;
    data = new float[ROWS][COLS];
  }

  Matrix(int rows_, int cols_, String s) {
    // Create a new Matrix of random floats
    ROWS = rows_;
    COLS = cols_;
    rows = rows_;
    cols = cols_;
    data = new float[rows][cols];
    for (int i = 0; i < rows; i++) {
      if (s.equals("identity")) {
        data[i][i] = 1.0;
      } else {
        for (int j = 0; j < cols; j++) {
          if (s.equals("zeros")) {
            data[i][j] = 0.0;
          } else if (s.equals("ones")) {
            data[i][j] = 1.0;
          } else if (s.equals("rand")) {
            data[i][j] = random(1.0);
          }
        }
      }
    }
  }

  Matrix(float[][] data_) {
    // Create a new Matrix from the array supplied
    ROWS = data_.length;
    COLS = data_[0].length;
    rows = data_.length;
    cols = data_[0].length;
    data = new float[rows][cols];
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        data[i][j] = data_[i][j];
      }
    }
  }

  Matrix newcopy() {
    // Return a copy of this Matrix
    mat = new Matrix(ROWS, COLS);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        mat.data[i][j] = data[i][j];
      }
    }
    return mat;
  }

  void copy(Matrix mat_) {
    // Copy the Matrix in place
    if ((mat_.rows > ROWS) || (mat_.cols > COLS)) { // expand the storage size
      this.resize(mat_.rows, mat_.cols);
    }
    for (int i = 0; i < mat_.rows; i++) {
      for (int j = 0; j < mat_.cols; j++) {
        data[i][j] = mat_.data[i][j];
      }
    }
    rows = mat_.rows;
    cols = mat_.cols;
  }

  Matrix newpcopy(int fromrow_, int torow_, int fromcol_, int tocol_) {
    // Return a partial copy of this Matrix
    mat = new Matrix(torow_ - fromrow_ + 1, tocol_ - fromcol_ + 1);
    for (int i = 0; i < torow_ - fromrow_ + 1; i++) {
      for (int j = 0; j < tocol_ - fromcol_ + 1; j++) {
        mat.data[i][j] = data[fromrow_ + i][fromcol_ + j];
      }
    }
    return mat;
  }

  void thispcopy(Matrix mat_, int fromrow_, int torow_, int fromcol_, int tocol_) {
    // Perform a partial copy of the Matrix in place
    if ((torow_ - fromrow_ + 1 > ROWS) || (tocol_ - fromcol_ + 1 > COLS)) { // expand the storage size
      this.resize(torow_ - fromrow_ + 1, tocol_ - fromcol_ + 1);
    }
    for (int i = 0; i < torow_ - fromrow_ + 1; i++) {
      for (int j = 0; j < tocol_ - fromcol_ + 1; j++) {
        data[i][j] = mat_.data[fromrow_ + i][fromcol_ + j];
      }
    }
    rows = torow_ - fromrow_ + 1;
    cols = tocol_ - fromcol_ + 1;
  }

  void replace(Matrix mat_, int row_, int col_) {
    // Replace a block of the Matrix with this
    for (int i = 0; i < mat_.rows; i++) {
      for (int j = 0; j < mat_.cols; j++) {
        data[row_ + i][col_ + j] = mat_.data[i][j];
      }
    }
  }

  void resize(int rows_, int cols_) {
    // Resize the active Matrix
    if ((rows_ == rows) && (cols_ == cols)) {
      return;
    }
    if ((rows_ > ROWS) || (cols_ > COLS)) { // expand the storage size
      ROWS = max(ROWS, rows_);
      COLS = max(COLS, cols_);
      data_ = new float[ROWS][COLS];
      // copy all old storage elements
      for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
          data_[i][j] = data[i][j];
        }
      }
      data = data_;
      rows = rows_;
      cols = cols_;
    } else { // change the active size
      rows = rows_;
      cols = cols_;
    }
  }

  float magnitude() {
    float s = 0.0;
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        s += data[i][j]*data[i][j];
      }
    }
    return sqrt(s);
  }

  void addi(int m) {
    // Add an int to each element of this Matrix
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        data[i][j] += m;
      }
    }
  }

  void addf(float m) {
    // Add a float to each element of this Matrix
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        data[i][j] += m;
      }
    }
  }

  void add(Matrix mat_) {
    // Add a Matrix to this
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        data[i][j] += mat_.data[i][j];
      }
    }
  }

  void subf(float m) {
    // Subtract a float from each element of this Matrix
    addf(-m);
  }

  void multf(float m) {
    // Multiply each element of this Matrix by a float
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        data[i][j] *= m;
      }
    }
  }

  void divf(float m) {
    // Divide each element of this Matrix by a float
    multf(1.0/m);
  }

  void normalize() {
    // Normalize each element of this Matrix
    float m = magnitude();
    if (m > 0.0) {
      divf(m);
    }
  }

  void limit(float max) {
    mag = magnitude();
    if (mag > max) {
      multf(max/mag);
    }
  }

  Matrix newtrans() {
    // Return a transpose of this Matrix
    mat = new Matrix(cols, rows);
    for (int i = 0; i < cols; i++) {
      for (int j = 0; j < rows; j++) {
        mat.data[i][j] = data[j][i];
      }
    }
    return mat;
  }

  void trans() {
    // Transpose this Matrix in place
    float data_tmp;
    int rows_tmp = rows;
    int cols_tmp = cols;
    if ((cols_tmp > ROWS) || (rows_tmp > COLS)) { // expand the storage size
      this.resize(cols_tmp, rows_tmp);
    }
    for (int i = 0; i < rows; i++) {
      for (int j = i + 1; j < cols; j++) {
        if ((i < cols) && (j < rows)) {
          data_tmp = data[i][j];
        }
        data[i][j] = data[j][i];
        if ((i < cols) && (j < rows)) {
          data[j][i] = data_tmp;
        }
      }
    }
  }

  void thistrans(Matrix mat_) {
    // Transpose the Matrix
    this.resize(mat_.cols, mat_.rows);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        data[i][j] = mat_.data[j][i];
      }
    }
  }

  Matrix newsum(Matrix mat_) {
    // Return the addition of a Matrix to this
    mat = new Matrix(rows, cols);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        mat.data[i][j] = data[i][j] + mat_.data[i][j];
      }
    }
    return mat;
  }

  void thissum(Matrix mat1_, Matrix mat2_) {
    // Add two Matrices
    this.resize(mat1_.rows, mat1_.cols);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        data[i][j] = mat1_.data[i][j] + mat2_.data[i][j];
      }
    }
  }

  void sum(Matrix mat_) {
    // Add a Matrix to this in place
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        data[i][j] += mat_.data[i][j];
      }
    }
  }

  Matrix newprod(Matrix mat_) {
    // Return the (pre-) multiplication of a Matrix to this
    mat = new Matrix(rows, mat_.cols);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < mat_.cols; j++) {
        for (int k = 0; k < cols; k++) {
          mat.data[i][j] += data[i][k]*mat_.data[k][j];
        }
      }
    }
    return mat;
  }

  void thisprod(Matrix mat1_, Matrix mat2_) {
    // Multiply two Matrices
    this.resize(mat1_.rows, mat2_.cols);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        data[i][j] = 0.0;
        for (int k = 0; k < mat1_.cols; k++) {
          data[i][j] += mat1_.data[i][k]*mat2_.data[k][j];
        }
      }
    }
  }

  Matrix newsqrtm22() {
    // Return the square root of this symmetric 2x2 matrix
    float m;
    float p;
    mat = new Matrix(2, 2);
    float tr_m = data[0][0] + data[1][1]; // trace
    float det_m = data[0][0]*data[1][1] - data[0][1]*data[1][0]; // determinant
    float q = sqrt(tr_m*tr_m - 4.0*det_m);
    // calculate eigenvalues
    float e1 = (tr_m + q)/2.0;
    float e2 = (tr_m - q)/2.0;
    if (e1 != e2) {
      m = (sqrt(e2) - sqrt(e1))/(e2 - e1);
      p = (e2*sqrt(e1) - e1*sqrt(e2))/(e2 - e1);
    } else {
      m = 1.0/(4.0*e1);
      p = sqrt(e1)/2.0;
    }
    mat.data[0][0] = m*data[0][0] + p;
    mat.data[1][0] = m*data[1][0];
    mat.data[0][1] = m*data[0][1];
    mat.data[1][1] = m*data[1][1] + p;

    return mat;
  }

  void thissqrtm22(Matrix mat_) {
    // Calculate the square root of the symmetric 2x2 matrix
    float m;
    float p;
    this.resize(2, 2);
    float tr_m = mat_.data[0][0] + mat_.data[1][1]; // trace
    float det_m = mat_.data[0][0]*mat_.data[1][1] - mat_.data[0][1]*mat_.data[1][0]; // determinant
    float q = sqrt(tr_m*tr_m - 4.0*det_m);
    // calculate eigenvalues
    float e1 = (tr_m + q)/2.0;
    float e2 = (tr_m - q)/2.0;
    if (e1 != e2) {
      m = (sqrt(e2) - sqrt(e1))/(e2 - e1);
      p = (e2*sqrt(e1) - e1*sqrt(e2))/(e2 - e1);
    } else {
      m = 1.0/(4.0*e1);
      p = sqrt(e1)/2.0;
    }
    data[0][0] = m*mat_.data[0][0] + p;
    data[1][0] = m*mat_.data[1][0];
    data[0][1] = m*mat_.data[0][1];
    data[1][1] = m*mat_.data[1][1] + p;
  }

  Matrix newchol(String t) {
    // Return the Cholesky decomposition
    mat = new Matrix(rows, rows);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j <= i; j++) {
        float s = 0.0;
        for (int k = 0; k < j; k++) {
          if (t.equals("lower")) {
            s += mat.data[i][k]*mat.data[j][k];
          } else {
            s += mat.data[k][i]*mat.data[k][j];
          }
        }
        if (i == j) {
          mat.data[i][i] = sqrt(data[i][i] - s);
        } else {
          if (t.equals("lower")) {
            mat.data[i][j] = (data[i][j] - s)/mat.data[j][j];
          } else {
            mat.data[j][i] = (data[j][i] - s)/mat.data[j][j];
          }
        }
      }
    }
    return mat;
  }

  void thischol(Matrix mat_, String t) {
    // Calculate the Cholesky decomposition: lower-diagonal form
    this.resize(mat_.rows, mat_.rows);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j <= i; j++) {
        float s = 0.0;
        for (int k = 0; k < j; k++) {
          if (t.equals("lower")) {
            s += data[i][k]*data[j][k];
          } else {
            s += data[k][i]*data[k][j];
          }
        }
        if (i == j) {
          data[i][i] = sqrt(mat_.data[i][i] - s);
        } else {
          if (t.equals("lower")) {
            data[i][j] = (mat_.data[i][j] - s)/data[j][j];
          } else {
            data[j][i] = (mat_.data[j][i] - s)/data[j][j];
          }
        }
      }
    }
  }

  Matrix newcholinv(String t) {
    // Return the inverse of the Cholesky decomposition: lower-diagonal form
    mat = this.newcopy(); // the cholesky decomposition
    for (int i = 0; i < rows; i++) {
      mat.data[i][i] = 1.0/this.data[i][i]; //diag
      for (int j = i + 1; j < rows; j++) {
        float s = 0.0;
        for (int k = i; k < j; k++) {
          if (t.equals("lower")) {
            s -= mat.data[j][k]*mat.data[k][i];
          } else {
            s -= mat.data[k][j]*mat.data[i][k];
          }
        }
        if (t.equals("lower")) {
          mat.data[j][i] = s/this.data[j][j];
        } else {
          mat.data[i][j] = s/this.data[j][j];
        }
      }
    }
    return mat;
  }

  void thischolinv(Matrix mat_, String t) {
    // Calculate the inverse of the Cholesky decomposition: lower-diagonal form
    this.copy(mat_);
    for (int i = 0; i < rows; i++) {
      data[i][i] = 1.0/mat_.data[i][i]; //diag
      for (int j = i + 1; j < rows; j++) {
        float s = 0.0;
        for (int k = i; k < j; k++) {
          if (t.equals("lower")) {
            s -= data[j][k]*data[k][i];
          } else {
            s -= data[k][j]*data[i][k];
          }
        }
        if (t.equals("lower")) {
          data[j][i] = s/mat_.data[j][j];
        } else {
          data[i][j] = s/mat_.data[j][j];
        }
      }
    }
  }

  void pprintm(int fromrow_, int torow_, int fromcol_, int tocol_) {
    // Print this active sub-matrix
    for (int i = fromrow_; i < torow_; i++) {
      for (int j = fromcol_; j < tocol_; j++) {
        print(" " + data[i][j]);
      }
      print("\n");
    }
    println("");
  }

  void printm() {
    // Print this active matrix
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        print(" " + data[i][j]);
      }
      print("\n");
    }
    println("");
  }

  void printmall() {
    // Print this storage Matrix
    for (int i = 0; i < ROWS; i++) {
      for (int j = 0; j < COLS; j++) {
        print(" " + data[i][j]);
      }
      print("\n");
    }
    println("");
  }
}
