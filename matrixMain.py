import numpy as np
class myMatrix:
    """
    This matrix class builds a matrix object for easier string visualization.
    """

    def __init__(self,matrix):
        """
        Initialize the matrix attributes
        :param matrix:
        """
        self.rows = len(matrix[0])
        self.columns = len(matrix)
        self.matrix = matrix
    def __str__(self):
        """
        Overrides python code to represent matrix in string form.
        :return: Returns a string object of the matrix
        """
        mString = ""
        for x in range(0,self.columns):
            mString += "\n"
            for y in range(0,self.rows):
                mString += str(self.matrix[x][y]) + " "
        print("\n")
        return mString


def mCopy(matrix):
    """
    Creates a copy of a given matrix.
    :param matrix: matrix object in the form of 2d list array [][].
    :return:
    """
    columns = len(matrix)
    rows = len(matrix[0])
    a = [[0 for x in range(rows)] for y in range(columns)]

    for i in range(columns):
        for j in range(rows):
            a[i][j] = matrix[i][j]


    return a


def getMDeterminantHelper(m, i, j):
    """
    Helper for the getMatrixDeterminant in 2x2 matrix form.
    :param m: matrix
    :param i: index i
    :param j: index j
    :return: list
    """
    return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]

def getMatrixDeternminant(m):
    """
    Recursive method to get the determinant of the matrix using 2x2 matrix.
    :param m:
    :return:
    """
    #Base Case
    if len(m) == 2:
        return m[0][0]*m[1][1]-m[0][1]*m[1][0]

    determinant = 0
    for c in range(len(m)):
        #Recusive method here using the helper
        determinant += ((-1)**c)*m[0][c]*getMatrixDeternminant(getMDeterminantHelper(m, 0, c))
    return determinant

def gaussJordan (matrix):
    """
    This is gaussJordan method to solve a matrix and return the identity matrix. This will solve the system of linear
    equations. Matrix must be in the standard format of a 2d array.
    :param matrix: List of list matrix [[]]
    :return: identity matirx or solution to system of linear equations
    """


    (columns, rows) = (len(matrix), len(matrix[0]))
    for y in range(0, columns):
        maxrow = y
        #Finds the max pivot
        for y2 in range(y + 1, columns):
            if abs(matrix[y2][y]) > abs(matrix[maxrow][y]):
                maxrow = y2
        (matrix[y], matrix[maxrow]) = (matrix[maxrow], matrix[y])
        # Destroys column
        for y2 in range(y + 1, columns):
            c = matrix[y2][y] / matrix[y][y]
            for x in range(y, rows):
                matrix[y2][x] -= matrix[y][x] * c
    # Backsubstitute for the gaussjordan algorithm
    for y in range(columns - 1, 0 - 1, -1):
        c = matrix[y][y]
        for y2 in range(0, y):
            for x in range(rows - 1, y - 1, -1):
                matrix[y2][x] -= matrix[y][x] * matrix[y2][y] / c
        matrix[y][y] /= c
        for x in range(columns, rows):  # Normalize row y
            matrix[y][x] /= c
    return matrix

def gauss (matrix):
    """
    This gauss solver will return a matrix in the upper-diagnol matrix form. This will not solve the system of
    linear equations. Must use the gaussJordan for that.
    """


    (column, rows) = (len(matrix), len(matrix[0]))
    for y in range(0, column):
        maxrow = y
        # Finds the max pivot
        for y2 in range(y + 1, column):
            if abs(matrix[y2][y]) > abs(matrix[maxrow][y]):
                maxrow = y2
        (matrix[y], matrix[maxrow]) = (matrix[maxrow], matrix[y])
        # Destroys column
        for y2 in range(y + 1, column):
            c = matrix[y2][y] / matrix[y][y]
            for x in range(y, rows):
                matrix[y2][x] -= matrix[y][x] * c
    return matrix


def transposeMatrix(matrix):
    """
    This will find the transpose for the given matrix. This method is used to find the inverse.
    :param matrix: List of list matrix [[]]
    :return:
    """
    rows = len(matrix[0])
    columns = len(matrix)
    tranposed = [[0 for x in range(rows)] for y in range(columns)]

    for i in range(columns):
        for j in range(rows):
            # Swaps the i and j columns
            tranposed[j][i] = matrix[i][j]
    return tranposed


def getMatrixInverse(m):
    determinant = getMatrixDeternminant(m)
    #Easy method for 2x2 matrix determinant
    if len(m) == 2:
        return [[m[1][1]/determinant, -1*m[0][1]/determinant],
                [-1*m[1][0]/determinant, m[0][0]/determinant]]

    #Create cofactors of the matrix to find the determanent bfore getting the transpose
    cofactors = []
    for r in range(len(m)):
        cofactorRow = []
        for c in range(len(m)):
            minor = getMDeterminantHelper(m, r, c)
            cofactorRow.append(((-1)**(r+c)) * getMatrixDeternminant(minor))
        cofactors.append(cofactorRow)

    cofactors = transposeMatrix(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c] = cofactors[r][c]/determinant
    return cofactors

# def augment(matrix,vector):
#     # Gets the sizes of the matrix and vector to ensure matrix properties can be applied.
#     mColumnSize = len(matrix)
#     mRowSize = len(matrix[0])
#     vColumnSize = len(vector)
#     vRowsSize = len(vector[0])
#
#     # Creates an empty matrix for the expected matrix n x m
#
#     newMatrix = [[0 for i in range(vRowsSize)] for j in range(mColumnSize+vColumnSize)]
#     # print("Rows:\t{},Columns:\t{}".format(mRowSize,mColumnSize))
#     if vRowsSize != mColumnSize:
#         print("This cannot be augmented. VectorRows:{}\t MatrixColumns{}".format(vRowsSize,mColumnSize))
#     else:
#         # print("Success. VectorRows:{}\t MatrixColumns{}".format(vRowsSize,mColumnSize))
#         # print(newMatrix)
#
#         #This will augment matrix and vector
#         for x in range(0,vRowsSize):
#             print("\n")
#             for y in range(0,mColumnSize):
#                 newMatrix[x][y] = matrix[x][y]
#                 # Sets the vector x+1 adds the column for the vector
#                 newMatrix[x+1][y] = vector[0][y]
#     # newMatrix =  myMatrix(newMatrix)
#
#
#     return newMatrix

def main():
    matrix = [[1,0,2,1],[2,-1,3,-1],[4,1,8,2]]
    matrix_12 = mCopy(matrix)
    matrix2 = [[1,-2,0],[-1,2,1],[0,-1,-2]]



    gj = gaussJordan(matrix)
    g = gauss(matrix_12)
    print("----gauss Jordan---\n")
    print(myMatrix(gj))
    print("The vector solution is: \n{}\n{}\n{}".format(gj[0][3],gj[1][3],gj[2][3]))
    print("---gauss----\n")
    print(myMatrix(g))
    print("-------\n")
    print("The Determinant of the gauss is {}".format(getMatrixDeternminant(matrix2)))
    print("-------\n")
    print("The Inverse of the gauss is {}".format(myMatrix(getMatrixInverse(matrix2))))

main()