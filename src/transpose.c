void transpose(double *x, int nrow, int ncol, double *y)
{
    int i, j;

    for (i = 0; i < ncol; i++)
        for (j = 0; j < nrow; j++)
            y[j*ncol + i] = x[i*nrow + j];
}
