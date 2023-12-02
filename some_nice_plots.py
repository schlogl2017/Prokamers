import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

#file = pd.read_csv('path/to/csv_file.csv')
#df = pd.DataFrame(file, columns=['bin_end', 'TTCAAG', 'dKO_Log2FC'])

#bin_end = df['bin_end']
#TTCAAG = df['TTCAAG']
#dKO_Log2FC = df['dKO_Log2FC']

#fig, ax = plt.subplots()
#ax2 = ax.twinx()

#sns.barplot(x=bin_end, y=dKO_Log2FC, ax=ax, color="blue", data=df)

#sns.scatterplot(x=bin_end, y=TTCAAG, ax=ax2, color="red", data=df)

#plt.title('Histone Position in TS559 vs dKO')
#plt.xlabel('Genomic Position (Bin = 1000nt)', fontsize=10)
#plt.xticks([])
#plt.ylabel('Log2 Fold Change', fontsize=10)

#plt.show()


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from PIL import Image

npoints=200
xRange=np.arange(0,npoints,1)
randomdata0=np.abs(np.random.normal(0,1,npoints))
randomdata1=np.random.normal(10,1,npoints)
axtick=[7,10,14]
ax2tick=[0,1.5,3]

fig0=plt.figure(0)
ax=fig0.gca()
ax2=ax.twinx()
sns.scatterplot(x=xRange,y=randomdata1,ax=ax)
ax.set_yticks(axtick)
ax.set_ylim([6,15])
ax2.set_yticks(ax2tick)
ax2.set_ylim([0,3.5])
plt.xticks([])

canvas0 = FigureCanvas(fig0)
s, (width, height) = canvas0.print_to_buffer()
X0 = Image.frombytes("RGBA", (width, height), s) #Contains the data of the first plot


fig1=plt.figure(1)
ax=fig1.gca()
ax2=ax.twinx()
sns.barplot(x=xRange,y=randomdata0,ax=ax2)
ax.set_yticks(axtick)
ax.set_ylim([6,15])
ax2.set_yticks(ax2tick)
ax2.set_ylim([0,3.5])
plt.xticks([])

canvas1 = FigureCanvas(fig1)
s, (width, height) = canvas1.print_to_buffer()
X1 = Image.frombytes("RGBA", (width, height), s) #Contains the data of the second plot

plt.figure(13,figsize=(10,10))
plt.imshow(Image.blend(X0,X1,0.5),interpolation='gaussian')
Axes=plt.gca()
Axes.spines['top'].set_visible(False)
Axes.spines['right'].set_visible(False)
Axes.spines['bottom'].set_visible(False)
Axes.spines['left'].set_visible(False)
Axes.set_xticks([])
Axes.set_yticks([])

plt.show()




