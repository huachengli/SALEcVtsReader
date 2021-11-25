import numpy as np
import matplotlib.pyplot as plt
import profile_plot

if __name__ == "__main__":

    PlotStep = 83
    PlotPath = "/home/huacheng/Documents/Github/data/s_dg_a70_bs2000I"
    data = profile_plot.read_data(PlotPath, PlotStep, 2, rx_=[-1.0e5, 0.4e5], ry_=[-0.5e5, 0.5e5], rz_=[-2.0e4, 2.0e4])
    # data = profile_plot.read_data(PlotPath, PlotStep, 2)
    data.load()
    fig = plt.figure(1, figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.set_aspect(1.0)
    cx = data.plot(ax)

    plt.show()