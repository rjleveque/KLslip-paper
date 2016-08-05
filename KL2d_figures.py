
from pylab import *
from clawpack.geoclaw import topotools, dtopotools
from clawpack.visclaw import colormaps
import os,sys

from scipy import stats
import KDEplots # local module for plotting kde

from numpy import random
import pandas as pd


import fetch_etopo  # will fetch topo file if not already here

subdir = 'figs_CSZ'
os.system('mkdir -p %s' % subdir)

def savefigp(fname):
    fname = fname.replace('png','jpg') # to output as jpg
    fname = os.path.join(subdir,fname)
    savefig(fname,dpi=250,bbox_inches='tight')
    print "Created ",fname

def compute_subfault_distances(fault):
    """
    Estimate the distance between subfaults i and j for every pair in the list
    fault.subfaults.

    :Inputs:
      -  *fault* of class dtopotools.Fault or some subclass,
    
    :Outputs:
      - *D* array of Euclidean distances based on longitudes, latitudes, and depths
      - *Dstrike* array of estimated distances along strike direction
      - *Ddip* array of estimated distances along dip direction
    with D**2 = Dstrike**2 + Ddip**2 to within roundoff.

    For each array, the [i,j] entry is distance from subfault i to j when
    ordered in the order the subfaults appear in the list fault.subfaults.

    Distance in dip direction based on differences in depth.  

    """

    import numpy
    from numpy import pi,sqrt,sin,cos,tan
    rad = pi/180.       # conversion factor from degrees to radians
    rr = 6.378e6        # radius of earth
    lat2meter = rr*rad  # conversion factor from degrees latitude to meters

    nsubfaults = len(fault.subfaults)
    D = numpy.zeros((nsubfaults,nsubfaults))
    Dstrike = numpy.zeros((nsubfaults,nsubfaults))
    Ddip = numpy.zeros((nsubfaults,nsubfaults))
    for i,si in enumerate(fault.subfaults):
        xi = si.longitude
        yi = si.latitude
        zi = si.depth
        for j,sj in enumerate(fault.subfaults):
            xj = sj.longitude
            yj = sj.latitude
            zj = sj.depth
            dx = abs(xi-xj)*cos(0.5*(yi+yj)*pi/180.) * lat2meter
            dy = abs(yi-yj) * lat2meter
            dz = abs(zi-zj)

            # Euclidean distance:
            D[i,j] = sqrt(dx**2 + dy**2 + dz**2)
            
            # estimate distance down-dip based on depths:
            dip = 0.5*(si.dip + sj.dip)
            ddip1 = dz / sin(dip*pi/180.)
            Ddip[i,j] = ddip1 
            if Ddip[i,j] > D[i,j]:
                # should not happen...
                if 0:
                    print "i,j,dx,dy,dz: ",i,j,dx,dy,dz
                    print "*** Ddip = %s, D = %s" % (Ddip[i,j], D[i,j])

            # compute distance in strike direction to sum up properly:
            dstrike2 = max(D[i,j]**2 - Ddip[i,j]**2, 0.)
            Dstrike[i,j] = sqrt(dstrike2)
                
    return D,Dstrike,Ddip
        

# read topo file for CSZ region:
#topo = topotools.Topography()
#topo.read('etopo1_-130_-120_35_55_4min.tt3',3)

topo = topotools.Topography()
topo.read('etopo1_-128_-122_40_50_1min.tt3',3)

# Crescent City location:
xcc = -124.1838
ycc = 41.7456

column_map = {"longitude":1, "latitude":2, "depth":3, "strike":4, 
              "length":5, "width":6, "dip":7}
defaults = {'rake': 90, 'slip':1.0}
coordinate_specification = 'top center'
input_units = {'slip': 'm', 'depth': 'km', 'length': 'km', 'width': 'km'}
rupture_type = 'static'
skiprows = 1
delimiter = ','

fault = dtopotools.Fault()
fault.read('CSZe01.csv', column_map, coordinate_specification,
           rupture_type,skiprows, delimiter, input_units, defaults)
print "There are %s subfaults" % len(fault.subfaults)

for s in fault.subfaults:
    s.longitude = s.longitude - 360.  # adjust to W coordinates
    
figure(figsize=(12,6))
ax = subplot(121);
fault.plot_subfaults(ax)
xticks(range(-128,-123));

# Now subdivide each subfault further

new_subfaults = []  # to accumulate all new subfaults

phi_plate = 60.  # angle oceanic plate moves clockwise from north, to set rake

for subfault in fault.subfaults:
    subfault.rake = subfault.strike - phi_plate - 180.
    # subdivide into nstrike x ndip subfaults, based on the dimensions of the
    # fault:
    nstrike = int(subfault.length/12000)
    ndip = int(subfault.width/10000)
    f = dtopotools.SubdividedPlaneFault(subfault, nstrike, ndip)
    new_subfaults = new_subfaults + f.subfaults

# reset fault.subfaults to the new list of all subfaults after subdividing:
new_fault = dtopotools.Fault(subfaults = new_subfaults)
n = len(new_fault.subfaults)
print "Subdivided fault has %s subfaults" % n

ax = subplot(122);
new_fault.plot_subfaults(ax)
xticks(range(-128,-123));


figure(figsize=(6,10))
ax = subplot(111)
contourf(topo.X,topo.Y,topo.Z,[0,20000],colors=[[.3,1,.3]])
fault.plot_subfaults(ax)
axis((-128,-122,40,50))

D, Dstrike, Ddip = compute_subfault_distances(new_fault)

# make correlation matrix:
# Gaussian with correlation lengths Lstrike and Ldip:
Lstrike = 400e3
Ldip = 40e3

print "Correlation lengths: Lstrike = %g, Ldip = %g" % (Lstrike,Ldip)
r = sqrt((Dstrike/Lstrike)**2 + (Ddip/Ldip)**2)
C = exp(-r)


lengths = array([s.length for s in fault.subfaults])
widths = array([s.width for s in fault.subfaults])
areas = lengths * widths
total_area = sum(areas)

Mw_desired = 9.0
Mo_desired = 10.**(1.5*Mw_desired + 9.05)
mean_slip = Mo_desired / (fault.subfaults[0].mu * total_area)
print "mean_slip %g meters required for Mw %s" % (mean_slip, Mw_desired)

# Turn this into a constant vector:
mean_slip = mean_slip * ones(n)

alpha = 0.5
sigma_slip = alpha * mean_slip


## Lognormal:
Cov_g = log((sigma_slip/mean_slip) * (C*(sigma_slip/mean_slip)).T + 1.)
mean_slip_g = log(mean_slip) - diag(Cov_g)/2.

## This should be the same:
Cov_g = log(alpha**2 * C + 1.)

# Find eigenvalues, and eigenvector matrix.
# Columns V[:,k] are eigenvectors.

print "Finding eigenmodes from %s by %s matrix C" % (n,n)
lam, V = eig(Cov_g)
    
eigenvals = real(lam)  # imaginary parts should be at rounding level
V = real(V)

# Sort eigenvalues:
i = list(argsort(lam))
i.reverse()
lam = lam[i]
V = V[:,i]


figure(figsize=(12,6))

ni = 1; nj = 4;
ax = axes((.1,.1,.15,.8))
contourf(topo.X,topo.Y,topo.Z,[0,20000],colors=[[.3,1,.3]])
contour(topo.X,topo.Y,topo.Z,[0],colors='g')
fault.plot_subfaults(ax)
axis((-128,-122,40,50))
xticks(fontsize=12)
yticks(fontsize=12)

cmap_slip = colormaps.make_colormap({0:'g',0.5:'w',1.:'m'})

for ii in range(ni):
    for jj in range(nj):
        pij = ii*nj + jj
        if jj<3:
            ax = axes((.3 + jj*0.15,.12,.12,.76))
            shrink = 0
        else:
            ax = axes((.3 + jj*0.15,.12,.17,.76))
            shrink = 0.
        V_amp = sqrt(sum(V[:,pij]**2))    # abs(V[:,pij]).max()
        #weight = sqrt(eigenvals[pij]) * V_amp / mean_amp
        for j,s in enumerate(new_fault.subfaults):
            s.slip = -V[j,pij] * 18.

        new_fault.plot_subfaults(ax,slip_color=True,cmin_slip=-1,cmax_slip=1,
                plot_box=0., cmap_slip=cmap_slip, colorbar_shrink=shrink)
        title('Mode %s' % pij, fontsize=18)
        axis('off')

savefigp('CSZ.png')







#==========================================
# Select only southern subfaults:

# read topo file for CA:
topo = topotools.Topography()
topo.read('etopo1_-126_-123_39_45_1min.tt3',3)

fault.subfaults = fault.subfaults[:8]

if 0:
    figure(figsize=(12,6))
    ax = subplot(121);
    fault.plot_subfaults(ax)
    xticks(range(-126,-123));
    contourf(topo.X,topo.Y,topo.Z,[0,20000],colors=[[.3,1,.3]])

# Now subdivide each subfault further

new_subfaults = []  # to accumulate all new subfaults

phi_plate = 60.  # angle oceanic plate moves clockwise from north, to set rake

for subfault in fault.subfaults:
    subfault.rake = subfault.strike - phi_plate - 180.
    # subdivide into nstrike x ndip subfaults, based on the dimensions of the
    # fault:
    nstrike = int(subfault.length/8000)
    ndip = int(subfault.width/8000)
    f = dtopotools.SubdividedPlaneFault(subfault, nstrike, ndip)
    new_subfaults = new_subfaults + f.subfaults

# reset fault.subfaults to the new list of all subfaults after subdividing:
new_fault = dtopotools.Fault(subfaults = new_subfaults)
n = len(new_fault.subfaults)
print "Subdivided fault has %s subfaults" % n

if 0:
    ax = subplot(122);
    new_fault.plot_subfaults(ax)
    xticks(range(-126,-123));




D, Dstrike, Ddip = compute_subfault_distances(new_fault)

# make correlation matrix:
# Gaussian with correlation lengths Lstrike and Ldip:
Lstrike = 130e3
Ldip = 40e3

print "Correlation lengths: Lstrike = %g, Ldip = %g" % (Lstrike,Ldip)
r = sqrt((Dstrike/Lstrike)**2 + (Ddip/Ldip)**2)
C = exp(-r)


lengths = array([s.length for s in fault.subfaults])
widths = array([s.width for s in fault.subfaults])
areas = lengths * widths
total_area = sum(areas)

Mw_desired = 8.8
Mo_desired = 10.**(1.5*Mw_desired + 9.05)
mean_slip = Mo_desired / (fault.subfaults[0].mu * total_area)
print "mean_slip %g meters required for Mw %s" % (mean_slip, Mw_desired)

# Turn this into a constant vector:
mean_slip = mean_slip * ones(n)

alpha = 0.5
sigma_slip = alpha * mean_slip


## Lognormal:
Cov_g = log((sigma_slip/mean_slip) * (C*(sigma_slip/mean_slip)).T + 1.)
mean_slip_g = log(mean_slip) - diag(Cov_g)/2.

## This should be the same:
Cov_g = log(alpha**2 * C + 1.)

# Find eigenvalues, and eigenvector matrix.
# Columns V[:,k] are eigenvectors.

print "Finding eigenmodes from %s by %s matrix C" % (n,n)
lam, V = eig(Cov_g)
    
eigenvals = real(lam)  # imaginary parts should be at rounding level
V = real(V)

# Sort eigenvalues:
i = list(argsort(lam))
i.reverse()
lam = lam[i]
V = V[:,i]


if 0:
    figure(figsize=(14,6))

    cmap_slip = colormaps.make_colormap({0:'g',0.5:'w',1.:'m'})

    ni = 1; nj = 4;
    ax = axes((.1,.1,.2,.8))
    contourf(topo.X,topo.Y,topo.Z,[0,20000],colors=[[.3,1,.3]])
    fault.plot_subfaults(ax)
    axis((-126,-123,40,44))

    for ii in range(ni):
        for jj in range(nj):
            pij = ii*nj + jj
            if jj<3:
                ax = axes((.35 + jj*0.14,.12,.12,.76))
                shrink = 0
            else:
                ax = axes((.35 + jj*0.14,.12,.17,.76))
                shrink = 0.7
            V_amp = sqrt(sum(V[:,pij]**2))    # abs(V[:,pij]).max()
            #weight = sqrt(eigenvals[pij]) * V_amp / mean_amp
            for j,s in enumerate(new_fault.subfaults):
                s.slip = V[j,pij] * 15.

            new_fault.plot_subfaults(ax,slip_color=True,cmin_slip=-1,cmax_slip=1,
                    plot_box=0., cmap_slip=cmap_slip, colorbar_shrink=shrink)
            title('Mode %s' % pij, fontsize=18)
            axis('off')

    savefigp('CSZmodes.png')



# More modes:

figure(figsize=(14,12))

cmap_slip = colormaps.make_colormap({0:'g',0.5:'w',1.:'m'})

ni = 2; nj = 4;
ax = axes((.1,.3,.2,.4))
contourf(topo.X,topo.Y,topo.Z,[0,20000],colors=[[.3,1,.3]])
contour(topo.X,topo.Y,topo.Z,[0],colors='g')
plot([xcc],[ycc],'wo')
plot([xcc],[ycc],'kx')
fault.plot_subfaults(ax)
axis((-126,-123,40,44))

for ii in range(ni):
    for jj in range(nj):
        pij = ii*nj + jj
        if ii==0:
            ax = axes((.35 + jj*0.14,.5,.12,.38))
        else:
            ax = axes((.35 + jj*0.14,.1,.12,.38))
            
        V_amp = sqrt(sum(V[:,pij]**2))    # abs(V[:,pij]).max()
        #weight = sqrt(eigenvals[pij]) * V_amp / mean_amp
        for j,s in enumerate(new_fault.subfaults):
            s.slip = V[j,pij] * 15.

        new_fault.plot_subfaults(ax,slip_color=True,cmin_slip=-1,cmax_slip=1,
                plot_box=0., cmap_slip=cmap_slip, colorbar_shrink=0)
        title('Mode %s' % pij,fontsize=18)
        axis('off')

savefigp('CSZmodes.png')



# Taper:

max_depth = 20000.
tau = lambda d: 1. - exp((d - max_depth)*20/max_depth)


## Realizations


def KL(z):
    KL_slip = 0.*mean_slip_g.copy()  # drop the mean slip and rescale later
    # add in the terms in the K-L expansion:  (dropping V[:,0])
    for k in range(1,len(z)):
        KL_slip += z[k] * sqrt(lam[k]) * V[:,k]
    
    ## Lognormal:
    KL_slip = exp(KL_slip)
    
    # Set the fault slip for the resulting realization:
    for j,s in enumerate(new_fault.subfaults):
        s.slip = KL_slip[j] * tau(s.depth)
        
    # Rescale to have desired magnitude:
    Mo = new_fault.Mo()
    KL_slip *= Mo_desired/Mo
    for j,s in enumerate(new_fault.subfaults):
        s.slip = KL_slip[j] * tau(s.depth)
    
    return KL_slip

# grid on which to compute deformation:
nx_dtopo = 181
ny_dtopo = 361
x_dtopo = linspace(-126,-123,nx_dtopo)
y_dtopo = linspace(39,45,ny_dtopo)

n_subfaults = len(new_fault.subfaults)
dZ = zeros((ny_dtopo, nx_dtopo, n_subfaults)) # to store sea floor deformation corresponding to each mode V[:,j]

for j in range(n_subfaults):
    sfault = dtopotools.Fault(subfaults = [new_fault.subfaults[j]])
    sfault.subfaults[0].slip = 1.
    dtopo = sfault.create_dtopography(x_dtopo,y_dtopo,times=[1.], verbose=False)
    sys.stdout.write('%i...' % j)
    sys.stdout.flush()
    dZ[:,:,j] = dtopo.dZ[0,:,:]



def PotentialEnergy(dZr):
    dy = 1./60. * 111.e3  # m
    dx = dy * cos(topo.Y * pi/180.)  # m
    grav = 9.81  # m/s^2
    rho_water = 1000  # kg/m^3
    eta = ma.masked_where(topo.Z>0, dZr)
    Energy = sum(eta**2 * dx * dy) * grav * rho_water * 1e-15  # PetaJoules
    return Energy



i1cc = find(dtopo.x<xcc).max()
j1cc = find(dtopo.y<ycc).max()
a1cc = (xcc-dtopo.x[i1cc])/(dtopo.x[i1cc+1]-dtopo.x[i1cc])
a2cc = (ycc-dtopo.y[j1cc])/(dtopo.y[j1cc+1]-dtopo.y[j1cc])
if (a1cc<0.) or (a1cc>1.) or (a2cc<0.) or (a2cc>1.):
    print '*** Interpolation to CC not correct!'

def dZ_CrescentCity(dZr):
    dzy1 = (1.-a1cc)*dZr[j1cc,i1cc] + a1cc*dZr[j1cc,i1cc+1]
    dzy2 = (1.-a1cc)*dZr[j1cc+1,i1cc] + a1cc*dZr[j1cc+1,i1cc+1]
    dzcc = (1.-a2cc)*dzy2 + a2cc*dzy1
    return dzcc

seed = 13579   # so random number generator gives repeatable results
random.seed(seed)

nterms = 60
nterms2 = 7
figure(figsize=(10,12))
for i in range(1,6):
    z = randn(nterms)
    KL_slip = KL(z)
    ax = subplot(4,5,i)
    new_fault.plot_subfaults(ax, slip_color=True, cmax_slip=20., plot_box=False, colorbar_shrink=0)
    axis('off')
    title('Realization %i\n %i terms' % (i,nterms), fontsize=14)
    #dtopo = new_fault.create_dtopography(x_dtopo,y_dtopo,times=[1.], verbose=False)
    dZr = dot(dZ,KL_slip)  # linear combination of dZ from unit sources
    ax = subplot(4,5,5+i)
    dtopotools.plot_dZ_colors(dtopo.X,dtopo.Y,dZr, axes=ax, cmax_dZ = 8., \
                              dZ_interval = 1., add_colorbar=False)
    ylim(39.5,44.5)
    plot([xcc],[ycc],'wo')
    plot([xcc],[ycc],'kx')
    title('E=%4.2f,\n dB=%5.2f' \
        % (PotentialEnergy(dZr),dZ_CrescentCity(dZr)), fontsize=14)
    axis('off')
    
    z = z[:nterms2]
    KL_slip = KL(z)
    ax = subplot(4,5,10+i)
    new_fault.plot_subfaults(ax, slip_color=True, cmax_slip=20., plot_box=False, colorbar_shrink=0)
    axis('off')
    title('%i terms' % nterms2, fontsize=14)
    #dtopo = new_fault.create_dtopography(x_dtopo,y_dtopo,times=[1.], verbose=False)
    dZr = dot(dZ,KL_slip)  # linear combination of dZ from unit sources
    ax = subplot(4,5,15+i)
    dtopotools.plot_dZ_colors(dtopo.X,dtopo.Y,dZr, axes=ax, cmax_dZ = 8., \
                              dZ_interval = 1., add_colorbar=False)
    ylim(39.5,44.5)
    plot([xcc],[ycc],'wo')
    plot([xcc],[ycc],'kx')
    title('E=%4.2f,\n dB=%5.2f' \
        % (PotentialEnergy(dZr),dZ_CrescentCity(dZr)), fontsize=14)
    axis('off')
    
savefigp('CSZrealizations.png')




### Densities computed from large number of realizations:


def test(ntrials = 10000, nterms=60):
    Energy = zeros(ntrials)
    Amplitude = zeros(ntrials)
    z_shore = zeros(ntrials)
    EtaMax = zeros(ntrials)
    
    zvals = zeros((ntrials,nterms+1))
    for j in range(ntrials):
        z = randn(nterms+1)  # choose random z for this realization
        zvals[j,:] = z
        KL_slip = KL(z)
        dZr = dot(dZ,KL_slip)  # linear combination of dZ from unit sources
        Energy[j] = PotentialEnergy(dZr)
        z_offshore = where(topo.Z < 0, dZr, 0.)
        Amplitude[j] = z_offshore.max() - z_offshore.min()
        z_shore[j] = dZ_CrescentCity(dZr)
        EtaMax[j] = z_offshore.max()
    return Energy, Amplitude, z_shore, EtaMax, zvals

random.seed(12345)
ntrials = 20000
print "Generating %s samples..." % ntrials
Energy, Amplitude, z_shore, EtaMax, zvals = test(ntrials = ntrials, nterms=60)

DepthProxy = EtaMax - z_shore
realizations = pd.DataFrame()
realizations['Energy'] = Energy
realizations['Amplitude'] = Amplitude
realizations['subsidence / uplift'] = z_shore
realizations['EtaMax'] = EtaMax
realizations['depth proxy'] = DepthProxy

#realizations.to_pickle('data/realizations_2d_' + str(nterms) + '.pkl')


random.seed(12345)
ntrials = 20000
nterms2 = 7
print "Generating %s samples..." % ntrials
Energy, Amplitude, z_shore, EtaMax, zvals = test(ntrials = ntrials, nterms=nterms2)

DepthProxy = EtaMax - z_shore
realizations2 = pd.DataFrame()
realizations2['Energy'] = Energy
realizations2['Amplitude'] = Amplitude
realizations2['subsidence / uplift'] = z_shore
realizations2['EtaMax'] = EtaMax
realizations2['depth proxy'] = DepthProxy

#realizations2.to_pickle('data/realizations_2d_' + str(nterms2) + '.pkl')

# Compare KDE plots:

Q1 = 'EtaMax'
Q2 = 'subsidence / uplift'

Nx = 200
Ny = 201

x0 = 3.
x1 = 12.

y0 = -2.
y1 = 4.

x = np.linspace(x0,x1,Nx)
y = np.linspace(y0,y1,Ny)

X,Y = np.meshgrid(x,y)

xy = np.vstack((X.flatten(),Y.flatten()))


## 60 terms:

rvals = vstack((realizations[Q1].T, realizations[Q2].T))
kde = stats.gaussian_kde(rvals)
p = kde.pdf(xy) * (x1-x0)*(y1-y0)/float(Nx*Ny)
rho = reshape(p,(Ny,Nx))

KDEplots.joint_plot(X,Y,rho,xname='EtaMax (m)',yname='subsidence/uplift (m)')
savefigp('joint_Eta_DBshore_60b.jpg')

## 7 terms:

rvals = vstack((realizations2[Q1].T, realizations2[Q2].T))
kde2 = stats.gaussian_kde(rvals)
p2 = kde2.pdf(xy) * (x1-x0)*(y1-y0)/float(Nx*Ny)
rho2 = reshape(p2,(Ny,Nx))

KDEplots.joint_plot(X,Y,rho2,xname='EtaMax (m)',yname='subsidence/uplift (m)')
savefigp('joint_Eta_DBshore_7b.jpg')


# Compare KDE plots:

Q1 = 'EtaMax'
Q2 = 'Energy'

Nx = 200
Ny = 201

x0 = 3.
x1 = 12.

y0 = 0.8
y1 = 1.6

x = np.linspace(x0,x1,Nx)
y = np.linspace(y0,y1,Ny)

X,Y = np.meshgrid(x,y)

xy = np.vstack((X.flatten(),Y.flatten()))


## 60 terms:

rvals = vstack((realizations[Q1].T, realizations[Q2].T))
kde = stats.gaussian_kde(rvals)
p = kde.pdf(xy) * (x1-x0)*(y1-y0)/float(Nx*Ny)
rho = reshape(p,(Ny,Nx))

KDEplots.joint_plot(X,Y,rho,xname='EtaMax (m)',yname='Energy (PetaJoules)')
savefigp('joint_Eta_Energy_60b.jpg')

## 7 terms:

rvals = vstack((realizations2[Q1].T, realizations2[Q2].T))
kde2 = stats.gaussian_kde(rvals)
p2 = kde2.pdf(xy) * (x1-x0)*(y1-y0)/float(Nx*Ny)
rho2 = reshape(p2,(Ny,Nx))

KDEplots.joint_plot(X,Y,rho2,xname='EtaMax (m)',yname='Energy (PetaJoules)')
savefigp('joint_Eta_Energy_7b.jpg')


