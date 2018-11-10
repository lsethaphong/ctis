function b = ssp_mxv ( m, n, nz_num, row, col, a, x )

%% SSP_MXV multiplies an SSP matrix times a vector.
%
%  Discussion:
%
%    The SSP storage format stores the row, column and value of each nonzero
%    entry of a sparse matrix.
%
%  Modified:
%
%    13 February 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, N, the number of rows and columns of the matrix.
%
%    Input, integer NZ_NUM, the number of nonzero elements in the matrix.
%
%    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column indices
%    of the nonzero elements.
%
%    Input, real A(NZ_NUM), the nonzero elements of the matrix.
%
%    Input, real X(N), the vector to be multiplied by A.
%
%    Output, real B(M), the product vector A*X.
%
  b(1:m) = 0.0E+00;

  for ( k = 1 : nz_num )

    i = row(k);
    j = col(k);
    b(i) = b(i) + a(k) * x(j);

  end
