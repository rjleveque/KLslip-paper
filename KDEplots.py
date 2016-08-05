
from pylab import *

def joint_plot(x,y,rho,xlimits=None,ylimits=None,xname=None,yname=None,fig=None):

    """
    Plot joint probability density rho along with marginal of each variable.
    x,y,rho should be 2d arrays of the same shape
    """

    # marginals:
    rho_x = rho.sum(0)
    rho_y = rho.sum(1)

    if fig is None:
        fig = figure(figsize=(8,8))

    axes([.1,.1,.55,.55])
    clines = rho.max()*linspace(.1,.9,9)
    contour(x,y,rho,clines,colors='k')
    if xlimits is not None:
        xlim(xlimits)
    if ylimits is not None:
        ylim(ylimits)

    if xname is not None:
        xlabel(xname,fontsize=14)
    if yname is not None:
        ylabel(yname,fontsize=14)

    axes([.1,.7,.55,.13])
    plot(x[0,:],rho_x,'k')
    if xlimits is not None:
        xlim(xlimits)
    ylim(0,1.1*rho_x.max())
    yticks([])

    axes([.7,.1,.13,.55])
    xticks([]); 
    plot(rho_y,y[:,0],'k')
    xlim(0,1.1*rho_y.max())
    if ylimits is not None:
        ylim(ylimits)

    return fig

