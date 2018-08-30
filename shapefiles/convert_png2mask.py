import numpy as np 
import matplotlib.pyplot as plt
from PIL import Image


im = Image.open('test.png')
imarray = np.array(im)
plt.imshow(imarray)
plt.show()