import plotly
import math
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('info_step_size_0.0035.csv')

fig = px.histogram(df, x=(df['resid']),nbins=50)
mean = round((np.mean(df['resid'])),5)
st = round((np.std(df['resid'])),5)
plt.text(0.08,50, r'$\mu=$'+str(mean)+' $\sigma=$ '+str(st))
fig.show() 
