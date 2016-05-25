//final_project_matrixInverseCalculator.cpp
//J. Long
//
//given an n x n matrix by the user, calculates the inverse of the matrix

//calculates the inverse of an n x n matrix using formula: inverse(A) = 1/det(A) * (adjoint(A))
//adjoint(A) = transpose(cofactor(A))

#include <iostream>
//allows use of vectors
#include <vector>
//allows use of power function
#include <cmath>
//allows use of setw() function
#include <iomanip>

using namespace std;

//function takes a matrix represented by a vector of vectors of doubles, as well as an int representing the size of the matrix
//and calculates the determinant, returning a double
double calcDeterminant(vector <vector <double> > matrix, int n);

//function takes a matrix represented by a vector of vectors of doubles, calculates the coFactor matrix,
//returning the coFactor matrix represented by a vector of vectors of doubles
vector <vector <double> > calcCoFactorMatrix(vector <vector <double> > matrix);

//function takes a matrix represented by a vector of vectors of doubles, calculates the transpose of the matrix,
//returning the transpose of the matrix represented by a vector of vectors of doubles
vector <vector <double> > calcTransposeMatrix(vector <vector <double> > matrix);

//function takes a matrix represented by a vector of vectors of doubles, and prints it out, returning nothing
void printMatrix(vector <vector <double> > matrix);

int main()
{
    //takes input for the number of rows/columns in the square matrix, only accepting a positive integer less than or equal to 6
    int n = 0;
    while ((n < 1)||(n > 6))
    {
        cout << "Enter the number of rows/columns in the n x n matrix." << endl;
        cout << "Ensure n is a positive integer less than or equal to 6." << endl;
        cin >> n;
    }
    //creates a vector of vectors of doubles to represent the input matrix
    vector <vector <double> > matrix;
    //the outer vector is set to the size that was input by the user
    matrix.resize(n);
    //looping through each element in the outer vector
    for (int i = 0; i < n; i++)
    {
        //each inner vector is set to the size that was input by the user, in effect creating a square matrix
        matrix[i].resize(n);
    }
    //loops through every row in the matrix
    for (int i = 0; i < n; i++)
    {
        //loops through every column in the matrix
        for (int j = 0; j < n; j++)
        {
            //asks for input at the position i+1, j+1
            cout << "Enter the value at row " << i+1 << ", column " <<  j+1 << endl;
            //stores this value in the matrix in the appropriate position
            cin >> matrix [i] [j];
        }
    }
    //prints out the inputted matrix by calling the printMatrix() function
    cout << "The matrix entered was:" << endl;
    cout << "********************************************************************************";
    printMatrix(matrix);
    cout << "********************************************************************************" << endl;
    //calculates the determinant by calling the calcDeterminant(function)
    //sends the inputted matrix and the size of the inputted matrix
    double det = calcDeterminant(matrix,n);
    //if the determinant is 0, the inverse does not exist, so prints an appropriate message and ends the program
    if (det == 0)
    {
        cout << "The inverse of this matrix does not exist because the determinant is 0." << endl;
        cout << "The matrix is a singular matrix." << endl;
        return 0;
    }
    //creates another matrix container using a vector of vectors of doubles, in order to store the inverse matrix
    vector <vector <double> > invMatrix;
    //the inverse matrix is sized to be the same as the inputted matrix
    invMatrix.resize(n);
    for (int i = 0; i < n; i++)
    {
        invMatrix[i].resize(n);
    }
    //loops through every row in the inverse matrix
    for (int i = 0; i < n; i++)
    {
        //loops through every column in the inverse matrix
        for (int j = 0; j < n; j++)
        {
            //sets the inverse matrix at row i, column j, equal to the reciprocal of the determinant, calculated previously,
            //multiplied by the value in the transpose of the cofactor matrix of the inputted matrix at row i, column j
            //calls the calcCoFactorMatrix() function inside of the calcTransposeMatrix() function
            invMatrix[i][j] = 1/det * calcTransposeMatrix(calcCoFactorMatrix(matrix))[i][j];
        }
    }
    //prints out the inverse matrix using the printMatrix() function
    cout << "The inverse matrix is:" << endl;
    cout << "********************************************************************************";
    printMatrix(invMatrix);
    cout << "********************************************************************************" << endl;
    return 0;
}

//function takes a matrix represented by a vector of vectors of doubles, as well as an int representing the size of the matrix
//and calculates the determinant, returning a double
double calcDeterminant(vector <vector <double> > matrix, int n)
{
    //creates a variable to store the determinant
    double det = 0;
    //if it is a 1x1 matrix, the determinant is the only value in the matrix
    if (n == 1)
    {
        det = matrix[0][0];
        //returns the determinant right away
        return det;
    }
    //if it is a larger matrix, calculate the determinant using recursive formula
    else if (n > 1)
    {
        //initial sum is 0
        double sum = 0;
        //loops through every column in the matrix
        for (int mColumn = 0; mColumn < n; mColumn++)
        {
            //for each column, make a new square matrix, sized one smaller than the current matrix
            vector <vector <double> > newmatrix;
            newmatrix.resize(n-1);
            for (int newMRow = 0; newMRow < n-1; newMRow++)
            {
                newmatrix[newMRow].resize(n-1);
            }
            //loops through each row in the new matrix
            for (int newMRow = 0; newMRow < n-1; newMRow++)   //loops through rows in new matrix
            {
                //second variable to refer to the column of the old matrix, acts as a loop starting at 0
                int mColumn2 = 0;
                //loops through each column in the new matrix
                for (int newMColumn = 0; newMColumn < n-1; newMColumn++)
                {
                    //if the column in the old matrix referred to in the outer loop is equal to the current column of the new matrix
                    if (mColumn == newMColumn)
                    {
                        //increment the variable to refer to the column number of the old matrix by 1, effectively skipping that column
                        mColumn2++;
                    }
                    //sets the new matrix at row newMRow and column newMColumn equal to the old matrix
                    //at row newMRow + 1 so that no value in the first row of the old matrix is reassigned, and column mColumn2
                    //effectively creates a new matrix with every element absent from the row and column of the element referred to
                    //in the first row of the outer matrix
                    newmatrix[newMRow][newMColumn] = matrix[newMRow+1][mColumn2];
                    //increases mColumn2 by 1, acting as a loop
                    mColumn2++;
                }
            }
            //adds the value of the matrix at position mColumn multiplied by (-1) to the power of the column number
            //multiplied by the determinant of the new matrix of size n-1
            sum += matrix[0][mColumn]*pow(-1,mColumn)*calcDeterminant(newmatrix, n-1);
            //sets the determinant equal to the sum
            det = sum;
        }
    }
    //returns the value of the determinant after the term has been added for each element in the 1st row of the old matrix
    return det;
}

//function takes a matrix represented by a vector of vectors of doubles, calculates the coFactor matrix,
//returning the coFactor matrix represented by a vector of vectors of doubles
vector <vector <double> > calcCoFactorMatrix(vector <vector <double> > matrix)
{
    //calculates the size of the external vector from the vector of vectors of doubles that was entered
    //this is the size of the matrix as only square matrices are allowed
    int n = matrix.size();
    //creates a new matrix of size n, to store each value in the cofactor matrix
    vector <vector <double> > coFactorMatrix;
    coFactorMatrix.resize(n);
    for (int cofMRow = 0; cofMRow < n; cofMRow++)
    {
        coFactorMatrix[cofMRow].resize(n);
    }
    //though not part of the definition of a cofactor, for purposes of the inverse calculator, the cofactor of any 1 x 1 matrix is 1
    //this allows the inverse of a 1 x 1 matrix to be calculated
    if (n == 1)
    {
        coFactorMatrix[0][0] = 1;
    }
    //if the size of the matrix is not 1 (ie bigger than 1, because negative matrix sizes are not allowed)
    else
    {
        //for each element in the matrix
        for (int mRow = 0; mRow < n; mRow++)
        {
            for (int mColumn = 0; mColumn < n; mColumn++)
            {
                //create a matrix of size one less than the inputted matrix
                vector <vector <double> > newmatrix;
                newmatrix.resize(n-1);
                for (int newMRow = 0; newMRow < n-1; newMRow++)
                {
                    newmatrix[newMRow].resize(n-1);
                }
                //using similar code as in the determinant calculator, creates a matrix with every element not in the row and column
                //referred to in the outer loop
                //creates a variable to reference the row of the old matrix
                int mRow2 = 0;
                //loops through each row of the new matrix
                for (int newMRow = 0; newMRow < n-1; newMRow++)
                {
                    //skips assigning the row of the old matrix if it is equal to the new row
                    if (mRow == newMRow)
                    {
                        mRow2++;
                    }
                    //creates a variable to reference the column of the old matrix
                    int mColumn2 = 0;
                    //loops through each column of the new matrix
                    for (int newMColumn = 0; newMColumn < n-1; newMColumn++)
                    {
                        //skips assigning the column of the old matrix if it is equal to the new column
                        if (mColumn == newMColumn)
                        {
                            mColumn2++;
                        }
                        //assigns the value of the old matrix to the new matrix, skipping columns and rows referred to in the outer loop
                        newmatrix[newMRow][newMColumn] = matrix[mRow2][mColumn2];
                        //increment mColumn2 by 1 to act as a looping mechanism
                        mColumn2++;
                    }
                    //increment mRow2 by 1 to act as a looping mechanism
                    mRow2++;
                }
                //calculates the determinant of the new matrix by calling the calcDeterminant() function
                double det = calcDeterminant(newmatrix, n-1);
                //sets the coFactorMatrix at the row and column referred to in the outer loop equal to the determinant of the new matrix
                //multiplied by (-1) to the power of the row number + the column number
                coFactorMatrix[mRow][mColumn] = det*pow(-1,mRow+mColumn);
            }
        }
    }
    //returns the entire cofactor matrix
    return coFactorMatrix;
}

//function takes a matrix represented by a vector of vectors of doubles, calculates the transpose of the matrix,
//returning the transpose of the matrix represented by a vector of vectors of doubles
vector <vector <double> > calcTransposeMatrix(vector <vector <double> > matrix)
{
    //calculates the size of the matrix, and creates a matrix of the same size to store the transpose of the matrix
    int n = matrix.size();
    vector <vector <double> > transposeMatrix;
    transposeMatrix.resize(n);
    for (int transposeMRow = 0; transposeMRow < n; transposeMRow++)
    {
        transposeMatrix[transposeMRow].resize(n);
    }
    //loop through the row of the new matrix
    for (int transposeMRow = 0; transposeMRow < n; transposeMRow++)
    {
        //loop through the column of the new matrix
        for (int transposeMColumn = 0; transposeMColumn < n; transposeMColumn++)
        {
            //assigns the new matrix the value of the old matrix at the reversed column and row index, essentially flipping the matrix
            transposeMatrix[transposeMRow][transposeMColumn] = matrix[transposeMColumn][transposeMRow];
        }
    }
    //returns the entire transpose matrix
    return transposeMatrix;
}

//function takes a matrix represented by a vector of vectors of doubles, and prints it out, returning nothing
void printMatrix(vector <vector <double> > matrix)
{
    //calculates the size of the outer vector, representing the size of the entire matrix as only square matrices are allowed
    int n = matrix.size();
    //loop through each row of the matrix
    for (int row = 0; row < n; row++)
    {
        //loop through each column of the matrix
        for (int column = 0; column < n; column++)
        {
            //prints the value of the matrix at the correct row and column, using setw() to format output into a table
            //setw() is set to 13 as that is the largest number of digits that the computer will take from the user
            cout << setw(13) << matrix[row][column];
        }
        //ends the line after printing out all the columns in the row
        cout << endl;
    }
}
