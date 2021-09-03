import pandas as pd
import numpy as np

df1 = pd.read_csv(r"C:\Users\angus\OneDrive\Desktop\urop\migration_bias\migration_bias_10x.csv")
df2 = pd.read_csv(r"C:\Users\angus\OneDrive\Desktop\urop\migration_bias\migration_bias_10x_left.csv")

df3 = pd.concat([df2, df1.iloc[:,3:]], axis = 1)

cols_names = list(df3.columns)
cols_names[0] = ''
df3.columns = cols_names
df3.index = df3.iloc[:,0]
df3 = df3.iloc[:,1:]

df3.to_csv(r"C:\Users\angus\OneDrive\Desktop\urop\migration_bias\migration_bias_10x_updated.csv")
