#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    if (a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }

    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            result.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return result;
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    if (a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }

    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            result.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
    return result;
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    if (a.cols != b.rows)
    {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }

    Matrix result = create_matrix(a.rows, b.cols);
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < b.cols; j++)
        {
            result.data[i][j] = 0;
            for (int k = 0; k < a.cols; k++)
            {
                result.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
    return result;
}

Matrix scale_matrix(Matrix a, double k)
{
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            result.data[i][j] = a.data[i][j] * k;
        }
    }
    return result;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix result = create_matrix(a.cols, a.rows);
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            result.data[j][i] = a.data[i][j];
        }
    }
    return result;
}

double det_matrix(Matrix a)
{
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }

    if (a.rows == 1)
        return a.data[0][0];
    if (a.rows == 2)
        return a.data[0][0] * a.data[1][1] - a.data[0][1] * a.data[1][0];
    double det = 0;
    // 对行列式第一行进行Laplace展开
    for (int p = 0; p < a.cols; p++) 
    {
        Matrix sub = create_matrix(a.rows - 1, a.cols - 1);
        double sgn = (p % 2 == 0) ? 1 : -1;
        for (int i = 1; i < a.rows; i++)
        {
            int sub_col = 0;
            for (int j = 0; j < a.cols; j++)
            {
                if (j == p) continue;
                sub.data[i - 1][sub_col++] = a.data[i][j];
            }
        }
        det += sgn * a.data[0][p] * det_matrix(sub);
    }
    return det;
}

Matrix inv_matrix(Matrix a) // 利用伴随矩阵求逆
{
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }

    double det = det_matrix(a);
    if (fabs(det) < 1e-8)
    {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    // 伴随矩阵adj
    int n = a.rows;
    Matrix adj = create_matrix(n, n); 

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            // (j,i) 位置的代数余子式
            Matrix sub = create_matrix(n - 1, n - 1);
            for (int k = 0, sub_i = 0; k < n; k++) // 行遍历j
            {
                if (k == j) continue;
                for (int l = 0, sub_j = 0; l < n; l++) // 列遍历i
                {
                    if (l == i) continue;
                    sub.data[sub_i][sub_j++] = a.data[k][l];
                }
                sub_i++;
            }

            double sgn = ((i + j) % 2 == 0) ? 1 : -1;
            double cofactor = sgn * det_matrix(sub); // 代数余子式A(j.i)的值
            adj.data[i][j] = cofactor; 
        }
    }

    return scale_matrix(adj, 1.0 / det);
}

int rank_matrix(Matrix a)
{
    int m = a.rows;
    int n = a.cols;
    
    Matrix temp = create_matrix(m, n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            temp.data[i][j] = a.data[i][j];
        }
    }

    int rank = 0;  
    double eps = 1e-8;  
    
    for (int col = 0; col < n; col++) {
        // 在当前列从rank行开始寻找主元
        int pivot_row = rank;
        while (pivot_row < m && fabs(temp.data[pivot_row][col]) < eps) {
            pivot_row++;
        }
        
        // 如果找到了非零元素
        if (pivot_row < m) {
            // 增加秩
            rank++;
            
            // 如果主元不在当前行，交换行
            if (pivot_row != rank - 1) {
                for (int j = 0; j < n; j++) {
                    double swap = temp.data[rank - 1][j];
                    temp.data[rank - 1][j] = temp.data[pivot_row][j];
                    temp.data[pivot_row][j] = swap;
                }
            }
            
            // 将主元行归一化，即将主元行整体除以主元，让主元化为1，接下来所有代码合起来即R_j = R_j - (a_ji/a_ii) * R_i的实现
            double pivot = temp.data[rank - 1][col];
            for (int j = col; j < n; j++) {
                temp.data[rank - 1][j] /= pivot;
            }
            
            // 消去其他行在当前列的元素
            for (int i = 0; i < m; i++) {
                if (i != rank - 1) {
                    double factor = temp.data[i][col];
                    for (int j = col; j < n; j++) {
                        temp.data[i][j] -= factor * temp.data[rank - 1][j];
                    }
                }
            }
        }
    }
    
    return rank;
}

double trace_matrix(Matrix a)
{
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    double trace = 0;
    for (int i = 0; i < a.rows; i++)
    {
        trace += a.data[i][i];
    }
    return trace;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}