#include <iostream>
#include <sstream>
#include <iomanip> // IO Manipulator
#include <cstdlib>
#include <vector>
#include <ctime>
#include <cmath>

#define DISPLAY_CONST 2

using namespace std;

bool findnegative();
int getspacer(double, bool);
void swap(double&, double&);
string ftos(double);

class Matrix
{
      // I won't try acessing this member outside the scope of the class relative methods, so it's private.
      vector<vector<double> > matrix;
      public:
        int row, column;
        // empty constructor
        Matrix() : row(0), column(0) {}
        
        // constructor
        Matrix(int x, int y)
        {
            if(x < 0 || y < 0)
            {
                cout << "Error: Size of matrix can not be less than 0.";
                row = column = 0;
            }
            else
            {
                row = x; column = y;
                matrix.resize(row); // resize makes a matrix of this size
                for (size_t i = 0; i < row; i++)
                  matrix[i].resize(column); // this one resizes the tiny matrixes inside, to create a vector of vectors acting like a 2d array.             
            }
        }
        
        // Takes values from user to fill
        void fill();
        void fill_random();
        void show_content();

        // To avoid const compile errors, when assigning a return value in a function that uses m(x,y), you get non-const assignment error.
        // The const after function name prevents it, as itself is const! Same was done on minor. As long as those values aren't modified.
        double access (const int x, const int y) const 
        {
            if(x > row || y > column)
            {
                cout << "Size index out of bounds!";
                return 0;
            }
            else
                return matrix[x][y]; 
        }

        // Easier access without constantly typing m.matrix[x][y]
        double & operator() (const int x, const int y) { return matrix[x][y]; }  
        
        // Multiplication with a scalar operators      
        Matrix & operator*(const double scalar);
        Matrix & operator*=(const double scalar);
        // Returns a new matrix without reference, as it's a local matrix in the function
        Matrix operator*(const Matrix& m);
        // This actually modifies your matrix, so it needs a reference.
        Matrix & operator*=(const Matrix& m);
        
        // Addition of matrices, 1st one creates a new matrix, 2nd one doesn't
        Matrix operator+(const Matrix& m);
        Matrix & operator+=(const Matrix& m);
        
        // Subtraction of matrices
        Matrix operator-(const Matrix& m);
        Matrix & operator-=(const Matrix& m);
        
        // Finds if there's a negative element in matrix -- needed for spacer
        bool findnegative();
        // Finds max element in matrix
        double findmax();
        double absfindmax();
        double findmax_length();
        // Transposes the matrix, saves result in the calling matrix
        void transpose();
        // Creates a diagonal matrix with ones on the diagonal of size (size).
        Matrix create_diagonal(const int size);
        // Finds minor of a matrix at user given slots, slots are from 1 to n+1 instead of 0 to n. (Easier to understand for user)
        Matrix minor(const int del_row, const int del_col) const;
        // Gets the Upper/Triangular matrixes together, return depends on t_type
        Matrix gettriangular(const Matrix& m, char t_type);
        // Finds upper/lower triangular matrixes of given matrix
        Matrix ltrianguar(const Matrix& m);
        Matrix utriangular(const Matrix& m);
        int checktriangular(const Matrix& m);
        // Exclusive to triangular matrices
        double det_exclusive(const Matrix& m);
        // Finds the actual determinant by LU decomposition
        double determinant(const Matrix& m);
        // Finds the inverse with Gaussian Elimination
        Matrix inverse(const Matrix& m);

        // Deconstructor
        ~Matrix() {};
};

Matrix& Matrix::operator*(const double scalar)
{
    size_t i, j;
    for(i = 0; i < row; i++)
      for(j = 0; j < column; j++)
        matrix[i][j] *= scalar;
    return *this;
}

Matrix& Matrix::operator*=(const double scalar)
{
    size_t i, j;
    for(i = 0; i < row; i++)
      for(j = 0; j < column; j++)
        matrix[i][j] *= scalar;
    return *this;
}

Matrix Matrix::operator*(const Matrix& m)
{
    double sum;
    Matrix temp(row, m.column);
    
    if(column != m.row)
      cout << "Matrix multiplication can't be performed, column of first matrix doesn't equal row of second.\n";
    
    for(size_t i = 0; i < row; i++)
    {
        for(size_t j = 0; j < m.column; j++)
        {
            for(size_t k = 0; k < m.row; k++)
                  sum += access(i,k)*m.access(k,j);
                  
            temp(i,j) = sum;
            sum = 0;
        }
    }
    return temp;
}

Matrix& Matrix::operator*=(const Matrix& m)
{
    double sum;
    Matrix temp = *this;
    if(column != m.row)
      cout << "Matrix multiplication can't be performed, column of first matrix doesn't equal row of second.\n";
    
    for(size_t i = 0; i < row; i++)
    {
        for(size_t j = 0; j < m.column; j++)
        {
            for(size_t k = 0; k < m.row; k++)
                  sum += temp.access(i,k)*m.access(k,j);
                  
            matrix[i][j] = sum;
            sum = 0;
        }
    }
    // Update the matrix
    column = m.column;
    matrix.resize(column);
    for(size_t i = 0; i < column; i++)
      matrix[i].resize(row);
          
    return *this;
}

Matrix Matrix::operator+(const Matrix& m)
{
    Matrix temp = *this;
    for(size_t i = 0; i < row; i++)
      for(size_t j = 0; j < column; j++)
        temp(i,j) += m.access(i,j);
    return temp;
}

Matrix& Matrix::operator+=(const Matrix& m)
{
    for(size_t i = 0; i < row; i++)
      for(size_t j = 0; j < column; j++)
        matrix[i][j] += m.access(i, j);
    return *this;
}

Matrix Matrix::operator-(const Matrix& m)
{
    Matrix temp = *this;
    for(size_t i = 0; i < row; i++)
      for(size_t j = 0; j < column; j++)
        temp(i,j) -= m.access(i,j);
    return temp;
}

Matrix& Matrix::operator-=(const Matrix& m)
{
    for(size_t i = 0; i < row; i++)
      for(size_t j = 0; j < column; j++)
        matrix[i][j] -= m.access(i,j);
    return *this;
}

void Matrix::fill()
{
    for (size_t i = 0; i < row; i++)
      for (size_t j = 0; j < column; j++)
      {
          cout << "Enter m_" << i+1 << j+1 << ": ";
          cin >> matrix[i][j];
      }
}

void Matrix::fill_random()
{
    for (size_t i = 0; i < row; i++)
      for (size_t j = 0; j < column; j++)
        matrix[i][j] = (rand()%101)*(2*(rand()%2)-1);
//          matrix[i][j] = double((rand()%1023+1))/double((rand()%1023+1));
}    

void Matrix::show_content()
{
     const int spacer = getspacer(absfindmax(), findnegative());
     for(size_t i = 0; i < row; i++)
     {
           for(size_t j = 0; j < column; j++)
             cout << setw(spacer+DISPLAY_CONST)<< setprecision(5) << matrix[i][j];
             
           cout << endl;
     }
     cout << endl;
}

double Matrix::findmax()
{
    double ext = access(0,0);
        for(size_t i = 0; i < row; i++)
          for(size_t j = 0; j < column; j++)
            if(ext < access(i,j))
               ext = access(i,j);
    return ext;        
}

double Matrix::absfindmax()
{
    double ext = matrix[0][0];
        for(size_t i = 0; i < row; i++)
          for(size_t j = 0; j < column; j++)
            if(ext < abs(matrix[i][j]))
               ext = abs(matrix[i][j]);
    return ext;        
}

double Matrix::findmax_length()
{
    double ext = matrix[0][0];
        for(size_t i = 0; i < row; i++)
          for(size_t j = 0; j < column; j++)
            if(ftos(ext).length() < ftos(access(i,j)).length())
               ext = access(i,j);
    return ext;
}

Matrix create_diagonal(const int size)
{
    Matrix temp = Matrix(size,size);
    for(size_t i = 0; i < size; i++)
        temp(i,i) = 1;
    return temp;
}

void Matrix::transpose()
{
    size_t i, j;
    int temp;
    vector<vector<double> > vtemp = matrix;
    matrix.resize(column);
    for(i = 0; i < column; i++)
      matrix[i].resize(row);
    
    // update elements
    for(i = 0; i < column; i++)
      for(j = 0; j < row; j++)
       matrix[i][j] = vtemp[j][i];
       
    // update members
    temp = column;
    column = row;
    row = temp;
}

// minor takes values from 1 to row+1, as opposed to other things. I don't know why I did it this way... Easier to understand I guess.
Matrix Matrix::minor(const int del_row, const int del_column) const
{
    Matrix temp;
    if(del_row > 0 && del_row <= row && del_column > 0 && del_column <= column)
    {
        temp = Matrix(row-1, column-1);
        
        // How this statement is created:
        // Assume matrix is 4x4, we want minor from the deletion of 2,2. i <= 4 - (2 >= 4) is 4. Same for column. So this will iterate as if it was i = 0 to < size, unless del_row == row in which case it will entirely skip the last iteration.
        // The actual index will become this => i = 1 => 1 - 1 - (1 > 2) = 0. Then, when i is 3 => 3 - 1 - (3 > 2) = 1. The 2nd index will be of 3rd index of my matrix!
        for(size_t i = 1; i <= (row - (del_row >= row)); i++)
          for(size_t j = 1; j <= (column - (del_column >= column)); j++)
          {
              // temp[i-1-(i>del_row][j-1-(j>del_column] = matrix[i-1][j-1];
              temp(i-1-(i > del_row), j-1-(j > del_column)) = access(i-1, j-1);
          }
    }
    else
    {
        cout << "Deletion index out of bounds!\n";
        return *this; // just to make sure it actually returns something...
    }
    return temp;
}

    // The logic really follows from multiplying this matrix, which equals m.
    /*
        [1      0   0   0] [u_11 u_12 u_13 u_14]
        [l_21   1   0   0] [0    u_22 u_23 u_24]
        [l_31 l_32  1   0] [0     0   u_33 u_34]
        [l_41 l_42 l_43 1] [0     0    0   u_44]
    */

Matrix gettriangular(const Matrix& m, char t_type)
{
    int max = m.row;
    Matrix temp = m;
    Matrix u(max, max); // I need to modify i >= j
    Matrix l(max, max);
    l = create_diagonal(max); // I need to modify j > i
    
    for(size_t i = 0; i < max; i++)
    {
        u(i,i) = temp(i,i); // by l assumption having 1 on all diagonals, u will have the same row as m's.
        // start at i+1 to fill unfilled entries
        for(size_t j = i+1; j < max; j++)
        {
            // ex: u(0,0) = m(0,0) => u(0,1) = m(0,1)/l(0,0) but l(0,0) = 1, so no need. This will be the same for all starting u.
            u(i,j) = temp(i,j);
            l(j,i) = temp(j,i)/u(i,i);
        }
        // This reduces the known values from temp for the matrix multiplication to recognize the newer values of u and l.
        for(size_t p = i+1; p < max; p++)
          for(size_t q = i+1; q < max; q++)
            temp(p,q) -= l(p,i)*u(i,q);        
    }
    if(t_type == 'u')
      return u;
    else
      return l;
}

Matrix ltriangular(const Matrix& m)
{
    Matrix temp = gettriangular(m, 'l');
    return temp;
}

Matrix utriangular(const Matrix& m)
{
    Matrix temp = gettriangular(m, 'u');
    return temp;
}

int checktriangular(const Matrix& m)
{
    int uresult, lresult;
    uresult = lresult = 1;
    int max = m.row;
    // checking the opposite case is better... We need zero on all (1,0), (2,0), (2,1) ... OR (0,1), (0, 2), (1,2)
    // uresult | lresult implies, if this fails BOTH tests for upper/lower triangularism, then it will fail. If it passes at least one test, the result will still be 1.
    for(size_t i = 1; i < max && (uresult | lresult); i++)
      for(size_t j = 0; j < i && (uresult | lresult); j++)
      {
          // First, check lower
          if(m.access(i,j) != 0)
            lresult = 0;
          // Then, check upper            
          if(m.access(j,i) != 0)
            uresult = 0;
      }
    return uresult|lresult;
}

double det_exclusive(const Matrix& m)
{
    double sum = 1;      
    for(size_t i = 0; i < m.row; i++)
      sum *= m.access(i,i);
    return sum;
}    

double determinant(const Matrix& m)
{
    double suml, sumu;
    if(m.row != m.column)
    {
        cout << "Only square matrices can have determinants.\n";
        return -1;
    }
    else
    {
        if(m.row == 1)
          return suml = m.access(0,0);
        else
        if(m.row == 2)        
          return suml = m.access(0,0)*m.access(1,1)-m.access(0,1)*m.access(1,0);
        else
        {
            // Now we have a 3x3 or higher matrix   
            // Shortcut if it's triangular (diagonal is also included)
            if(checktriangular(m))
                return det_exclusive(m);
            // sums from lower and upper triangular matrices
            int max = m.row;
            suml = sumu = 1;
            // get LU decomposition
            Matrix u = utriangular(m);
            Matrix l = ltriangular(m);
            
            for(size_t i = 0; i < max; i++)
            {
                // Another shortcut, don't bother if you find a 0.
                if(!u(i,i) || !l(i,i))
                  return 0; 
                sumu *= u(i,i);
                suml *= l(i,i);
            }
            
            return sumu*suml;
        }
    }
    
}
// The old, ugly determinant method
/*
double determinant(const Matrix& m)
{
    double sum = 0;
    if(m.row != m.column)
    {
        cout << "Only square matrices can have determinants.\n";
        return -1;
    }
    else
    {
        if(m.row == 1)
          sum = m.access(0,0);
        else
        if(m.row == 2)
          sum = m.access(0,0)*m.access(1,1)-m.access(0,1)*m.access(1,0);
        else
        {
            // Add shortcuts to check for diagonal matrix case.
            
            for(int i = 0; i < m.column; i++)
            {
                // start iterating from the top left corner, then recursively find the determinant by reducing the matrix
                Matrix temp = m.minor(1, i+1);
                sum += (2*((i + 1)%2) - 1)*m.access(0,i)*determinant(temp);
            }
        }
    }
    return sum;
}
*/

/*
            // Old ugly inverse method
            // Using adjugate matrix method
            det = determinant(m);
            for(int i = 0; i < m.row; i++) 
              for(int j = 0; j < m.column; j++)
              {
                  // 0.0 at the end takes care of the signed zero problem, which is -0 appearing as a result. A nice trick.
                  temp(i,j) = determinant(m.minor(i+1,j+1))*(1-(2*((i+j)%2)))+0.0;
              }
            temp *= 1/det; 
*/
Matrix inverse(const Matrix& m)
{
    Matrix temp;
    if(m.row != m.column)
    {
        cout << "Only square matrices can have inverses!\n";
        return temp;
    }
    else
    if(!determinant(m))
    {
        cout << "This is a singular matrix! It has determinant of 0.\n";
        return temp;
    }
    else
    {
        Matrix temp(m.row, m.column);
        double det;
        if(m.row == 1)
            temp(0,0) = 1/m.access(0,0);
        else
        if(m.row == 2)
        {
            // shorter method for 2x2 matrix
            double p;
            det = determinant(m);
            p = m.access(0,0);
            temp(0,0) = m.access(1,1);
            temp(1,1) = p;
            temp(0,1) = (-1)*m.access(0,1); temp(1,0) = (-1)*m.access(1,0);
            temp *= 1/det; 
        }
        else
        {
            // Using Gauss Elimination Method
            temp = create_diagonal(m.row); // create the identity matrix
            Matrix mtemp = m; // copy m
            int max = m.row;
            
            for(int j = 0; j < max; j++)
            {
                // element (j,j) should be non-zero, if not swap lower rows.
                int i;
                for(i = j; i < max && mtemp(i,j) == 0; i++) {} // skip when this condition is satisfied
                if(i != j)
                {
                    // swap rows
                    for(int p = 0; p < max; p++)
                    {
                        swap(mtemp(j,p), mtemp(i,p));
                        swap(temp(j,p), temp(i,p));
                    }
                }
                // eliminate non-zero values at other rows in column j
                double elim;
                for(i = 0; i < max; i++)
                {
                    if(i != j)
                    {
                        // eliminate value at column j row i
                        if(mtemp(i,j) != 0)
                        {
                            elim = - mtemp(i,j)/mtemp(j,j);
                            
                            // add this to eliminate the value
                            for(int p = 0; p < max; p++)
                            {
                                mtemp(i,p) += elim*mtemp(j,p);
                                temp(i,p) += elim*temp(j,p);
                            }
                        }
                    }
                    else
                    {
                        // make (j,j) = 1, to do that divide each value on rows with mtemp(j,j)
                        elim = mtemp(j,j);
                        for(int p = 0; p < max; p++)
                        {
                            mtemp(i,p) /= elim;
                            temp(i,p) /= elim;
                        }
                    }
                }
            }
        }
        return temp;
    }
}

string ftos(double value) 
{
  // ostringstream is under std namespace
  ostringstream o;
  if (!(o << value))
    return "";
  return o.str();
}

// Uses string conversion to detect the minus sign, otherwise doesn't really work for -0.000001 etc
bool Matrix::findnegative()
{
     string temp;
     for(size_t i = 0; i < row; i++)
       for(size_t j = 0; j < column; j++)
       {
           temp = ftos(access(i,j));
           if(temp[0] == '-')
             return true;
           else
             return false;
       }
}

// finds digits after the dot by putting the number in a string and getting it's length, but by creating a substring after the . is found
int getspacer(double num, bool isnegative)
{
    int digit = 0;
    double temp = num;
    // converts to string to detect the -
    string snum = ftos(num);
    // add in the - sign for spacer
    if(isnegative == true)
      digit += DISPLAY_CONST;    
    // this finds the amount of digits before the .
    while(temp >= 1)
    {
        temp /= 10;
        digit++;
    }
    // this finds the amount of digits after the .
    // convert the number to string, use find_next_of to search for . character, create a substring after it and get it's length.
    snum = ftos(temp);
    snum.substr(snum.find_first_of("."));
    digit += snum.length();
    return digit;
}

void swap(double &a, double &b)
{
    double temp = a;
    a = b;
    b = temp;
}

// MATRIX CLASS RELATED DEFINITIONS END HERE

int main(int x, int y)
{
    int option;
    srand(time(NULL));
    do
    {
        cout << "Type the number of the operation you'd like to do with the matrices." << endl;
        cout << "1.  Scalar multiplication\n2.  Matrix Multiplication\n3.  Matrix Addition\n4.  Matrix Subtraction\n5.  Transpose Matrix\n";
        cout << "6.  Minor\n7.  Lower Triangular Matrix\n8.  Upper Triangular Matrix\n9.  Check if triangular\n10. Determinant\n11. Inverse\n";
        cout << "12. Make random matrix, find determinant, transpose and inverse\n0.  Quit\nEnter your option: ";
        cin >> option;
        switch (option)
        {
               case 1:
               {
                    int scalar;
                    cout << "Create the size of your matrix: ";
                    cin >> x >> y;
                    Matrix mat1(x, y);
                    mat1.fill();
                    cout << "This is your matrix\n\n";
                    mat1.show_content();
                    cout << "\nEnter the scalar to multiply: ";
                    cin >> scalar;
                    mat1 *= scalar;
                    cout << "This is your matrix now\n\n";
                    mat1.show_content();
               }
               break;
               case 2:
               {
                    cout << "Create the size of the first matrix: ";
                    cin >> x >> y;
                    Matrix mat1(x, y);
                    cout << "Create the size of the second matrix: ";
                    cin >> x >> y;
                    Matrix mat2(x, y);
                    
                    cout << "Enter the elements of the first matrix" << endl;
                    mat1.fill();
                    cout << "Enter the elements of the second matrix" << endl;
                    mat2.fill();
                    mat1 *= mat2;
//                    Matrix result(mat1.row, mat2.column);
//                    result = mat1*mat2;
                    cout << "\nThe resulting matrix is:\n\n";
//                    result.show_content();
                    mat1.show_content();
               }
               break;
               case 3:
               {
                    cout << "Create the size of the first matrix: ";
                    cin >> x >> y;
                    Matrix mat1(x, y);
                    cout << "Create the size of the second matrix: ";
                    cin >> x >> y;
                    Matrix mat2(x, y);
                    
                    if(mat1.row != mat2.row || mat1.column != mat2.column)
                      cout << "Matrix addition can't be performed, matrices aren't equal sized.\n";
                    else
                    {
                        cout << "Enter the elements of the first matrix" << endl;
                        mat1.fill();
                        cout << "Enter the elements of the second matrix" << endl;
                        mat2.fill();
    //                    Matrix result(x, y);
    //                    result = mat1 + mat2;
    //                    result.show_content();
                        mat1 += mat2;
                        mat1.show_content();
                    }
               }
               break;
               case 4:
               {
                    cout << "Create the size of the first matrix: ";
                    cin >> x >> y;
                    Matrix mat1(x, y);
                    cout << "Create the size of the second matrix: ";
                    cin >> x >> y;
                    Matrix mat2(x, y);
                    
                    if(mat1.row != mat2.row || mat1.column != mat2.column)
                      cout << "Matrix subtraction can't be performed, matrices aren't equal sized.\n";
                    else
                    {
                        cout << "Enter the elements of the first matrix" << endl;
                        mat1.fill();
                        cout << "Enter the elements of the second matrix" << endl;
                        mat2.fill();
    //                    Matrix result(x, y);
    //                    result = mat1 - mat2;
    //                    result.show_content();
                        mat1 -= mat2;
                        mat1.show_content();
                    }
               }
               break;   
               case 5:      
               {
                    cout << "Create the size of the matrix: ";
                    cin >> x >> y;
                    Matrix mat1(x, y);
                    mat1.fill();
                    cout << "This is your matrix\n\n";
                    mat1.show_content();
                    mat1.transpose();
                    cout << "This is your matrix now\n\n";
                    mat1.show_content();               
               }
               break;
               case 6:
               {
                    cout << "Create the size of the matrix: ";
                    cin >> x >> y;
                    Matrix mat1(x, y);
                    mat1.fill();
                    cout << "This is your matrix\n\n";
                    mat1.show_content();
                    cout << "Enter the row and column to delete: ";
                    cin >> x >> y;
                    cout << "Your minor is\n\n";
                    Matrix mat2 = mat1.minor(x, y);
                    mat2.show_content();            
               }
               break;
               case 7:
               {
                    cout << "Create the size of the matrix: ";
                    cin >> x >> y;
                    Matrix mat1(x, y);
                    mat1.fill();
                    cout << "This is your matrix\n\n";
                    mat1.show_content();
                    Matrix result = ltriangular(mat1);
                    cout << "This is your lower triangular matrix\n\n";
                    result.show_content();                                
               }
               break;
               case 8:
               {
                    cout << "Create the size of the matrix: ";
                    cin >> x >> y;
                    Matrix mat1(x, y);
                    mat1.fill();
                    cout << "This is your matrix\n\n";
                    mat1.show_content();
                    Matrix result = utriangular(mat1);
                    cout << "This is your upper triangular matrix\n\n";
                    result.show_content();                                
               }
               break;
               case 9:
               {
                    cout << "Create the size of the matrix: ";
                    cin >> x >> y;
                    Matrix mat1(x, y);
                    mat1.fill();
                    cout << "This is your matrix\n\n";
                    mat1.show_content();
                    (checktriangular(mat1)) ? cout << "\nYour matrix is triangular!\n\n" : cout << "\nYour matrix is NOT triangular!\n\n";
               }
               break;                             
               case 10:
               {
                    cout << "Create the size of the matrix: ";
                    cin >> x >> y;
                    Matrix mat1(x, y);
                    mat1.fill();
                    cout << "This is your matrix\n\n";
                    mat1.show_content();
                    cout << "Your determinant is: " << determinant(mat1) << "\n\n";               
               }
               break;
               case 11:
               {
                    cout << "Create the size of the matrix: ";
                    cin >> x >> y;
                    Matrix mat1(x, y);
                    mat1.fill();
                    cout << "This is your matrix\n\n";
                    mat1.show_content();
                    Matrix result(x, y);
                    cout << "The inverse of this matrix is:\n\n";
                    result = inverse(mat1);
                    result.show_content();          
               }
               break;
               case 12:
               {
                    x = 32 - rand()%30;
                    cout << "This is your matrix\n\n";
                    Matrix mat1(x, x);
                    mat1.fill_random();
                    mat1.show_content();
                    cout << "\n\nThis is the determinant: " << determinant(mat1);
                    cout << "\nThis is the transpose\n\n";
                    Matrix temp = mat1;
                    temp.transpose();
                    temp.show_content();
                    cout << "\n\nThis is the inverse\n\n";
                    temp = inverse(mat1);
                    temp.show_content();
               }
               break;
        }
    } while(option);
    
    system("PAUSE");
    return 0;
}
