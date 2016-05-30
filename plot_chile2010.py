
from pylab import *
import os
from clawpack.geoclaw import topotools,dtopotools
from clawpack.clawutil.data import get_remote_file

subdir = 'figures'
os.system('mkdir -p %s' % subdir)

def savefigp(fname):
    fname = os.path.join(subdir,fname)
    savefig(fname, bbox_inches='tight')
    print "Created ",fname

filename = 'pacific_shorelines_east_4min.npy'
url = 'http://www.geoclaw.org/topo/' + filename
get_remote_file(url=url, output_dir='.', force=True, verbose=True)
shore = load(filename)

fault_geometry_file = 'chile2010_usgs.txt'
column_map = {"latitude":0, "longitude":1, "depth":2, "slip":3, "rake":4, "strike":5, "dip":6}
defaults = {'length': 30, 'width':20}
coordinate_specification = 'top center'
input_units = {'slip': 'cm', 'depth': 'km', 'length': 'km', 'width': 'km'}
rupture_type = 'static'
skiprows = 11
delimiter = None # white space is the default

fault = dtopotools.Fault()
fault.read(fault_geometry_file, column_map, coordinate_specification,
           rupture_type,skiprows, delimiter, input_units, defaults)
print "There are %s subfaults" % len(fault.subfaults)

x,y = fault.create_dtopo_xy()
dtopo = fault.create_dtopography(x,y,verbose=True)

figure(figsize=(7,5))
ax1 = subplot(111)
fault.plot_subfaults(axes=ax1,slip_color=True)
axis([-75,-70,-39,-32]);
savefigp('chile2010_slip.png')

figure(figsize=(7,5))
ax1 = subplot(111)
dtopo.plot_dZ_colors(2., axes=ax1)
plot(shore[:,0]-360., shore[:,1], 'g', linewidth=2)
axis([-75,-70,-39,-32]);
savefigp('chile2010_dz.png')

