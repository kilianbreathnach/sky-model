import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import numpy as np


def compare(data, model, res, chi, filename):

    head = np.percentile(data, 99.9)
    toe = np.percentile(data, 0.1)
    chirange = np.percentile(data, 5.0)
    resrange = np.percentile(data, 50.0)

    width, height = np.shape(data)
    figw, figh = (width/60.)*2, (height/60.)*2

    fig = plt.figure(figsize=(figw, figh))
#    fig = plt.figure(dpi=60.)

    axr = fig.add_subplot(221)
    axr.set_axis_off()
    axr.set_title(r"Data", fontsize=36)

    caxr = axr.imshow(data, cmap=plt.cm.gray, vmin=toe, vmax=head)
    divr = make_axes_locatable(plt.gca())
    cbarr = divr.append_axes("right", "5%", pad="3%")
    fig.colorbar(caxr, cax=cbarr)

    axf = fig.add_subplot(222)
    axf.set_axis_off()
    axf.set_title(r"Model Image", fontsize=36)

    caxf = axf.imshow(model, cmap=plt.cm.gray, vmin=toe, vmax=head)
    divf = make_axes_locatable(plt.gca())
    cbarf = divf.append_axes("right", "5%", pad="3%")
    fig.colorbar(caxf, cax=cbarf)

    axres = fig.add_subplot(223)
    axres.set_axis_off()
    axres.set_title(r"Residual Image",fontsize=36)

    caxres = axres.imshow(res, cmap=plt.cm.gray, vmin=-resrange, vmax=resrange)
    divres = make_axes_locatable(plt.gca())
    cbarres = divres.append_axes("right", "5%", pad="3%")
    fig.colorbar(caxres, cax=cbarres)

    axchi = fig.add_subplot(224)
    axchi.set_axis_off()
    axchi.set_title(r"Chi Image", fontsize=36)

    caxchi = axchi.imshow(chi, cmap=plt.cm.gray, vmin=-chirange, vmax=chirange)
    divchi = make_axes_locatable(plt.gca())
    cbarchi = divchi.append_axes("right", "5%", pad="3%")
    fig.colorbar(caxchi, cax=cbarchi)

    plt.tight_layout()

    fig.savefig(filename)

