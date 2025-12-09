import pandas as pd 
df = pd.read_csv("~/pacakges-set/expresion-matrix.csv",index_col = 0)
print(f"original shape :{df.shape}")
df_transposed = df.T
print (f"Transposed shape :{df_transposed.shape}")
df_transposed.to_csv("~/pacakges-set/expr_matrix_transposed.csv")
print(df.transposed.head())
