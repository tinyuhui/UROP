import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(r"C:\Users\angus\OneDrive\Desktop\urop\migration_bias.csv", index_col=0)

import seaborn as sns
cm = sns.light_palette("green", as_cmap=True)

df2 = df.style.background_gradient(cmap=cm)
