from pylab import *
import sys, os

from clawpack.geoclaw import dtopotools 

import pandas as pd
from numpy import random

# this way of importing doesn't change matplotlib plot style:
import seaborn.apionly as sns 

cmap_black = sns.dark_palette('black', as_cmap=True)

subdir = 'figures'
os.system('mkdir -p %s' % subdir)

def savefigp(fname):
    fname = os.path.join(subdir,fname)
    savefig(fname, bbox_inches='tight')
    print "Created ",fname
    
# Parameters to adjust...
dip = 13.     # dip angle in degrees
depth0 = 10e3  # depth of top of fault in meters

W = 100e3  # width of fault plane (down-dip) in meters
x = linspace(0,W,201)  # meters
x = x[:-1]  # drop last point so x values are top of subfaults
xkm = x / 1e3  # km for plotting

xkm_shore = 75

ndip = len(x)
m = ndip # dimension of matrices


cdip_fraction = 0.4  # correlation length
cdip = cdip_fraction*W  # correlation length in dip

mean_slip = 10. # desired average slip (meters)
L = 1000e3  # length of fault to use in computing deformation (very long)

rigidity = 3.55e10
Mo = L*W*mean_slip * rigidity
Mw = 2./3. * (log10(Mo) - 9.05)
print "mean magnitude Mw = %5.2f" % Mw

# Rigidity was calculated via:
#Mw_desired = 9.0
#Mo_desired = 10.**(1.5*Mw_desired + 9.05)
#rigidity = Mo_desired / (L*W*mean_slip)

def corr(x1,x2):
    c = exp(-abs(x1-x2)/cdip)  # exponential
    #c = exp(-((x1-x2)/cdip)**2) # Gaussian
    return c

# Correlation matrix with no sigma's:
Chat = eye(m)
for i in range(m):
    for j in range(m):
        Chat[i,j] = corr(x[i],x[j])

taper = 'exp_depth'

if taper=='none':
    tau = lambda x: ones(x.shape)

elif taper=='cospower':
    power = 10     # used in taper tau
    tau = lambda x: 1. - (1. - 0.5*(1. + cos(2*pi*(x+0.5*W)/W)))**power

elif taper=='cospower_downdip':
    power = 10     # used in taper tau
    tau = lambda x: where(x<W/2, 1., 1. - (1. - 0.5*(1. + cos(2*pi*(x+0.5*W)/W)))**power)
        
    
elif taper=='WangHe':
    # Wang and He (2008) taper:
    broadness = 0.25  #0.25
    qskew = 0.65
    def delta(xp):
        d1 = (12./qskew**3) * xp**2 * (qskew/2. - xp/3.)
        #d2 = (12./(1.-qskew)**3) * xp**2 * ((1.-qskew)/2. - (1.-xp)/3.)
        #d = where(xp <= qskew, d1, d2)
        dq = 2. 
        d2 = dq + (12./(1.-qskew)**3) * \
              ((xp**3 / 3. - xp**2 *(1+qskew)/2. + qskew*xp) \
             - (qskew**3 / 3. - qskew**2*(1+qskew)/2. + qskew**2))
        d = where(xp <= qskew, d1, d2)
        #import pdb; pdb.set_trace()
        return d
    def tau(x):
        xp = (x-x[0])/(x[-1]-x[0])
        #import pdb; pdb.set_trace()
        return delta(xp) * (1. + sin(broadness * pi * delta(xp)))

elif taper=='exp_depth':
    max_depth = 32500.
    tau_depth = lambda d: 1. - exp((d - max_depth)*20/max_depth)
    tau = lambda x: tau_depth(depth0 + x*sin(dip*pi/180.))


if 0:
    plot(xkm, tau(x))
    xlabel('km down-dip')
    title('taper')
    ylim(-0.1,1.1);

total_slip = mean_slip * len(x)
mu = total_slip * tau(x) / sum(tau(x)) # mean slip with taper

alpha = 0.75
sigma = alpha * mu
S = diag(sigma)
C = dot(dot(S,Chat),S)

lam, V = eig(C)
# C is symmetric so eigenvalues should be real
print "Max imaginary part of lam: ",abs(imag(lam)).max()
lam = real(lam)
V = real(V)

# Normalize eigenvectors to have 2-norm = 1:
# Should be normalized already!
for k in range(m):
    V[:,k] = V[:,k] / norm(V[:,k],2)

# Negate modes if necessary to facilitate comparisons:
# Choose sign so they all are initially increasing...
for k in range(m):
    if V[10,k] < 0:
        V[:,k] = -V[:,k]

# Sort eigenvalues:
i = list(argsort(lam))
i.reverse()
lam = lam[i]
V = V[:,i]

kplot=20
figure(figsize=(5,4))
loglog(range(1,kplot+1),sqrt(lam[:kplot]), 'o',label='sqrt(lambda)')
k = linspace(1,kplot,kplot)
numer = int(round(0.5*sqrt(lam[0])))

plot(k,numer/k,'r',label='%i/k' % numer)
legend(fontsize=15)
xticks(fontsize=15); yticks(fontsize=15);
xlim(1,20);
title("First %i eigenvalues of Covariance matrix" % kplot, fontsize=15);
savefigp('eigenvalues.jpg')

figure(figsize=(10,5))
plot(xkm,tau(x)/norm(tau(x),2), 'k--', label='taper')
plot(xkm,V[:,0],'b',label='k = 0')
plot(xkm,V[:,1],'r',linewidth=2,label='k = 1')
plot(xkm,V[:,2],'g',linewidth=2,label='k = 2')
plot(xkm,V[:,3],'m',linewidth=1,label='k = 3')
plot(xkm,V[:,4],'c',linewidth=1,label='k = 4')
legend(loc='lower right',fontsize=13)
xlim(0,140)
xlabel('km down-dip')
title('Eigenmodes')
savefigp('eigenmodes.jpg')

# Replace first eigenvector by mean:
V[:,0] = mu

#input_units = {'slip': 'm', 'depth': 'km', 'length': 'km', 'width': 'km'}
coordinate_specification = "top center"

# subfault parameters for full fault width:
subfault = dtopotools.SubFault()
subfault.coordinate_specification = coordinate_specification
subfault.longitude = 0.
subfault.latitude = 0.
subfault.depth = depth0
subfault.strike = 0.
subfault.length = 1000.e3  # very long since we want 1d cross section
subfault.width = W
subfault.dip = dip
subfault.rake = 90.
subfault.slip = 1.  # will be reset for each realization

# split up into fault with ndip subfaults
#base_fault = dtopotools.Fault([subfault], input_units, coordinate_specification) # needed to set units
fault = dtopotools.SubdividedPlaneFault(subfault, nstrike=1, ndip=ndip)

# grid on which to compute deformation:
n_dtopo = 1001  
x_dtopo = linspace(-2,2,n_dtopo)
xkm_dtopo = x_dtopo * 111.  # convert from degrees to km
y_dtopo = array([-1,0.,1.])  # for 1d slice

nterms = 20
dZ = zeros((n_dtopo,nterms+1)) # to store sea floor deformation corresponding to each mode V[:,j]

for j in range(nterms+1):
    for k,s in enumerate(fault.subfaults):
        s.slip = V[k,j]

    dtopo = fault.create_dtopography(x_dtopo,y_dtopo,times=[1.], verbose=False)
    print "Applying Okada to V[:,%i] to compute dZ[:,%i]" % (j,j)
    dZ[:,j] = dtopo.dZ[0,1,:]

def KL(z):
    #total_slip = 10. * len(x)
    slip = mu  # sqrt(lam[0])*V[:,0].copy()
    for k in range(1,len(z)):
        slip = slip + z[k]*sqrt(lam[k])*V[:,k]
    #slip = where(slip>0, slip, 0.)
    multiplier = 1.  # total_slip / norm(slip,1)
    slip = slip * multiplier
    z_dtopo = dZ[:,0].copy()
    for k in range(1,len(z)):
        z_dtopo += z[k]*sqrt(lam[k])*dZ[:,k]
    z_dtopo = z_dtopo * multiplier
    return slip, z_dtopo
        

figure(figsize=(14,3))
subplot(121)
slip = mu
plot(xkm, slip, 'k', lw=2,label="mean slip")
legend(fontsize=12)
#plot(xkm,mu,'g')
ylim(-5,25)

Mo = L*W*sum(abs(slip))/len(slip) * rigidity
Mw = 2./3. * (log10(Mo) - 9.05)
text(130,5,"Mw = %5.2f" % Mw, fontsize=12, color='k')
xlim(-50,200)
xlabel('km down-dip')
title('slip on fault plane')

subplot(122)
z_dtopo = V[:,0]
plot(xkm_dtopo, dZ[:,0], 'k', lw=2,label="mean deformation")
legend(fontsize=12)
xlim(-50,200)
ylim(-5,10)
xlabel('km from trench towards shore')
title('sea floor deformation')
savefigp('sample_realization_00.jpg')
random.seed(12333)  #12345 

nterms = 20
for rno in range(1,10):
    z = randn(nterms+1)
    slip, z_dtopo = KL(z)

    nterms3 = 3
    slip3, z_dtopo3 = KL(z[:nterms3])

    figure(figsize=(14,3))
    subplot(121)
    plot(xkm,slip, 'b', lw=2,label="%i terms" % nterms)
    plot(xkm,slip3,'r', lw=2,label="%i terms" % nterms3)
    legend(fontsize=12)
    #plot(xkm,mu,'g')
    ylim(-5,25)

    Mo = L*W*sum(abs(slip))/len(slip) * rigidity
    Mw = 2./3. * (log10(Mo) - 9.05)
    text(130,5,"Mw = %5.2f" % Mw, fontsize=12, color='b')
    Mo = L*W*sum(abs(slip[:nterms3]))/len(slip[:nterms3]) * rigidity
    Mw = 2./3. * (log10(Mo) - 9.05)
    text(130,2,"Mw = %5.2f" % Mw, fontsize=12, color='r')

    xlim(-50,200)
    xlabel('km down-dip')
    title('slip on fault plane')

    subplot(122)
    plot(xkm_dtopo, z_dtopo, 'b', lw=2, label="%i terms" % nterms)
    plot(xkm_dtopo, z_dtopo3, 'r', lw=2,label="%i terms" % nterms3)
    #plot(xkm_dtopo, z_dtopo - z_dtopo3, 'g')
    legend(fontsize=12)
    #plot(xkm_dtopo, dZ[:,0], 'g')
    xlim(-50,200)
    ylim(-5,10)
    xlabel('km from trench towards shore')
    title('sea floor deformation')
    savefigp('sample_realization_%s.jpg' % str(rno).zfill(2))


if 0:
    # check that results agree with running Okada after KL sum:
    for k,s in enumerate(fault.subfaults):
        s.slip = slip[k]
    dtopo = fault.create_dtopography(x_dtopo,y_dtopo,times=[1.], verbose=False)
    plot(xkm_dtopo, dtopo.dZ[0,1,:],'g')



j_shore = find(xkm_dtopo < xkm_shore).max()
print "with xkm_shore = %g, the shore is at j = %i in the dtopo.dZ arrays of length %i" \
        % (xkm_shore, j_shore, dZ.shape[0])
dZ_shore_mean = dZ[j_shore,0]
print "dZ at the shore for the mean slip is %g" % dZ_shore_mean
b = dZ[j_shore,:] * sqrt(lam[:dZ.shape[1]])
#print "dZ at shore for any z is then %g + b^Tz where z[0]=0 and " % dZ[j_shore,0]
#print "   b = ",b

def normal_density(x, mean=0, variance=1):
    sigma = sqrt(variance)
    print "sigma = ",sigma
    return 1./(sigma*sqrt(2*pi)) * exp(-(x-mean)**2 / (2*variance))
dzs = linspace(-5,5,1001)

dzs_rho1 = normal_density(dzs, dZ_shore_mean, sum(b[1:2]**2))
dzs_rho2 = normal_density(dzs, dZ_shore_mean, sum(b[1:3]**2))
dzs_rho3 = normal_density(dzs, dZ_shore_mean, sum(b[1:4]**2))
dzs_rho20 = normal_density(dzs, dZ_shore_mean, sum(b[1:21]**2))

figure()
plot(dzs, dzs_rho1, 'm', label='1 terms')
plot(dzs, dzs_rho2, 'g', label='2 terms')
plot(dzs, dzs_rho3, 'r', label='3 terms')
plot(dzs, dzs_rho20, 'b', label='20 terms')

xlim(-4,4)
ylim(0,1)
legend()
xlabel('Meters')
title("True density for dZ at shore")
savefigp("DBshore-true.jpg")


def test(ntrials = 10000, nterms=10, xkm_shore=300):
    grav = 9.81 # m/s**2
    rho_water = 1000.   # kg/m**3
    dx = x[1] - x[0]  # for integrating potential energy
    Energy = zeros(ntrials)
    Amplitude = zeros(ntrials)
    z_shore = zeros(ntrials)
    EtaMax = zeros(ntrials)
    
    scorr = zeros(V[:,0].shape)
    zvals = zeros((ntrials,nterms+1))
    for j in range(ntrials):
        z = randn(nterms+1)  # choose random z for this realization
        zvals[j,:] = z
        slip, z_dtopo = KL(z)
        scorr += (slip[50] - mu[50])*(slip - mu) / (sigma[50]*sigma)
        z_offshore = where(xkm_dtopo < xkm_shore, z_dtopo, 0.)
        Energy[j] = sum(z_offshore**2) * grav * rho_water * \
            dx*L * 1e-15  # PetaJoules
        Amplitude[j] = z_offshore.max() - z_offshore.min()
        z_shore[j] = z_dtopo[j_shore]  # uplift or subsidence at shoreline
        EtaMax[j] = z_offshore.max()
    return Energy, Amplitude, z_shore, EtaMax, zvals


xkm_shore = 75
random.seed(12345)
ntrials = 20000
print "Generating %s samples..." % ntrials
Energy, Amplitude, z_shore, EtaMax, zvals = test(ntrials = ntrials, nterms=20, xkm_shore=xkm_shore)

DepthProxy = EtaMax - z_shore
realizations = pd.DataFrame()
realizations['Energy'] = Energy
realizations['Amplitude'] = Amplitude
realizations['subsidence / uplift'] = z_shore
realizations['EtaMax'] = EtaMax
realizations['depth proxy'] = DepthProxy


random.seed(12345)
ntrials = 20000
print "Generating %s samples..." % ntrials
Energy1, Amplitude1, z_shore1, EtaMax1, zvals1 = test(ntrials = ntrials, nterms=1, xkm_shore=xkm_shore)
realizations1 = pd.DataFrame()
realizations1['Energy'] = Energy1
realizations1['Amplitude'] = Amplitude1
realizations1['subsidence / uplift'] = z_shore1
realizations1['EtaMax'] = EtaMax1
DepthProxy1 = EtaMax1 - z_shore1
realizations1['depth proxy'] = DepthProxy1


random.seed(12345)
ntrials = 20000
print "Generating %s samples..." % ntrials
Energy2, Amplitude2, z_shore2, EtaMax2, zvals2 = test(ntrials = ntrials, nterms=2, xkm_shore=xkm_shore)
realizations2 = pd.DataFrame()
realizations2['Energy'] = Energy2
realizations2['Amplitude'] = Amplitude2
realizations2['subsidence / uplift'] = z_shore2
realizations2['EtaMax'] = EtaMax2
DepthProxy2 = EtaMax2 - z_shore2
realizations2['depth proxy'] = DepthProxy2


random.seed(12345)
ntrials = 20000
print "Generating %s samples..." % ntrials
Energy3, Amplitude3, z_shore3, EtaMax3, zvals3 = test(ntrials = ntrials, nterms=3, xkm_shore=xkm_shore)
realizations3 = pd.DataFrame()
realizations3['Energy'] = Energy3
realizations3['Amplitude'] = Amplitude3
realizations3['subsidence / uplift'] = z_shore3
realizations3['EtaMax'] = EtaMax3
DepthProxy3 = EtaMax3 - z_shore3
realizations3['depth proxy'] = DepthProxy3


close('all')

def hazard_curve(realizations, zetai):
    counts = zeros(zetai.shape)
    for j in range(len(zetai)):
        i = realizations['depth proxy'] > zetai[j]
        counts[j] = sum(i)
    prob = counts / len(realizations)

    return prob

zetai = linspace(0,12,121)
prob = hazard_curve(realizations, zetai)
prob3 = hazard_curve(realizations3, zetai)
prob2 = hazard_curve(realizations2, zetai)
prob1 = hazard_curve(realizations1, zetai)
figure(figsize=(10,5))
semilogy(zetai, prob, 'b', label='20 terms')
semilogy(zetai, prob3, 'r', label='3 terms')
semilogy(zetai, prob2, 'g', label='2 terms')
semilogy(zetai, prob1, 'm', label='1 terms')
legend(loc='lower left')
title('Hazard curve for depth proxy')
xlabel('Exceedance value (meters)')
ylabel('probability')
savefigp('hazcurves1.jpg')


def plot_gaussian_contours():
    z1 = linspace(-3, 3, 100)
    Z1,Z2 = meshgrid(z1,z1)
    G = 1./(2*pi) * exp(-(Z1**2 + Z2**2)/2.)
    contour(Z1,Z2,G,5, colors='k')
    axis('scaled')

figure()
j = DepthProxy3 > 8
plot_gaussian_contours()
plot(zvals3[j,1],zvals3[j,2],'.')
savefigp('z12_bigdepth.jpg')


if 0:
    # try scatter plot where color indicates z[3]
    # not very illuminating...
    z3 = zvals3[j,3]
    c = where(z3>1,'r','b')
    c = where(z3<-1,'g',c)
    scatter(zvals3[j,1],zvals3[j,2],marker='.',color=c)


figure()
j = Energy3 > 9.5
plot_gaussian_contours()
plot(zvals3[j,1],zvals3[j,2],'.')
savefigp('z12_bigenergy.jpg')

## seaborn plots -- importing change matplotlib style
#import seaborn as sns

figure()
sns.kdeplot(realizations['subsidence / uplift'],label='20 terms')
sns.kdeplot(realizations3['subsidence / uplift'],label='3 terms')
plot(dzs, dzs_rho20, 'r', label='true density')
legend()
xlim(-4,4)
ylim(0,1.0)
title("Kernel density estimates for dZ at shore")
xlabel('Meters')
savefigp('DBshore-samples.jpg')

cmap_black = sns.dark_palette('black', as_cmap=True)


figure(figsize=(5,5))
#sns.jointplot('EtaMax','Energy',data=realizations,kind='kde',
#        xlim=(2,10),ylim=(0.4e7,1.1e7), stat_func=None)

g = sns.JointGrid('EtaMax','Energy',data=realizations,\
    xlim=(2,10),ylim=(4,11))
g = g.plot_joint(sns.kdeplot, cmap=cmap_black)
g = g.plot_marginals(sns.kdeplot, color='k', shade=True)
savefigp('joint_Eta_Energy_20.jpg')

figure(figsize=(5,5))
#sns.jointplot('EtaMax','Energy',data=realizations3,kind='kde',
#       xlim=(2,10),ylim=(0.4e7,1.1e7), stat_func=None)

g = sns.JointGrid('EtaMax','Energy',data=realizations3,\
    xlim=(2,10),ylim=(4,11))
g = g.plot_joint(sns.kdeplot, cmap=cmap_black)
g = g.plot_marginals(sns.kdeplot, color='k', shade=True)

savefigp('joint_Eta_Energy_3.jpg')


#figure(figsize=(5,5))
#sns.jointplot('EtaMax','subsidence / uplift',data=realizations,kind='kde',
#        xlim=(2,10),ylim=(-3,3), stat_func=None)
g = sns.JointGrid('EtaMax','subsidence / uplift',data=realizations,\
    xlim=(2,10),ylim=(-3,3))
g = g.plot_joint(sns.kdeplot, cmap=cmap_black)
g = g.plot_marginals(sns.kdeplot, color='k', shade=True)
savefigp('joint_Eta_DBshore_20.jpg')

#figure(figsize=(5,5))
#sns.jointplot('EtaMax','subsidence / uplift',data=realizations3,kind='kde',
#        xlim=(2,10),ylim=(-3,3), stat_func=None)
g = sns.JointGrid('EtaMax','subsidence / uplift',data=realizations3,\
    xlim=(2,10),ylim=(-3,3))
g = g.plot_joint(sns.kdeplot, cmap=cmap_black)
g = g.plot_marginals(sns.kdeplot, color='k', shade=True)
savefigp('joint_Eta_DBshore_3.jpg')


figure(figsize=(5,5))
sns.kdeplot(realizations['Energy'], label='20 terms')
sns.kdeplot(realizations3['Energy'], label='3 terms')
sns.kdeplot(realizations2['Energy'], label='2 terms')
sns.kdeplot(realizations1['Energy'], label='1 terms')
xlim(0,20)
xlabel('PetaJoules')
title('Potential energy density')
savefigp('kde_Energy.jpg')


figure(figsize=(5,5))
sns.kdeplot(realizations['EtaMax'], label='20 terms')
sns.kdeplot(realizations3['EtaMax'], label='3 terms')
sns.kdeplot(realizations2['EtaMax'], label='2 terms')
sns.kdeplot(realizations1['EtaMax'], label='1 terms')
xlabel('Meters')
title('EtaMax density')
savefigp('kde_EtaMax.jpg')


