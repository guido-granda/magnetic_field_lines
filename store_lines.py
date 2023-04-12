## code to generate hdf5 files containing the mag lines ##
from data_helper import *
from mag_lines_mod import *
import numpy as np
import h5py
from yt.units import *

# location of data cube in hdf5 format #
loc='./'
# location of output #
out='./'

#input and output names
basename='reg1_cg_'
out_name='reg1_cg_lines_'

n_s=50
n_e=50
d_n=1

nlines=12 # number of lines per plane dimension

for i in range(n_s,n_e+1,d_n):
    print('i: {}'.format(i))
    f=h5py.File(loc+basename+str(i).zfill(3)+'.h5','r')
    x=np.array(f['/grid/x']) # cm
    y=np.array(f['/grid/y']) # cm
    z=np.array(f['/grid/z']) # cm
    dens=np.array(f['/grid/dens']) # g/cm³
    dx=np.array(f['/grid/dx'][0][0][0]) # cm
    bx=np.array(f['/grid/magx'])/1e-6 # uG
    by=np.array(f['/grid/magy'])/1e-6 # uG
    bz=np.array(f['/grid/magz'])/1e-6 # uG
    
    f.close()

    x_min=np.amin(x)#*cm).in_units('pc')
    x_max=np.amax(x)#*cm).in_units('pc')
    y_min=np.amin(y)#*cm).in_units('pc')
    y_max=np.amax(y)#*cm).in_units('pc')
    z_min=np.amin(z)#*cm).in_units('pc')
    z_max=np.amax(z)#*cm).in_units('pc')
    print("x_min, x_max: {} ,{}".format(x_min,x_max))
    print("y_min, y_max: {} ,{}".format(y_min,y_max))
    print("z_min, z_max: {} ,{}".format(z_min,z_max))
    print("dx: {}".format(dx))
    nx=int(round((x_max-x_min+dx)/dx))
    ny=int(round((y_max-y_min+dx)/dx))
    nz=int(round((z_max-z_min+dx)/dx))
    print('nx: {}'.format(nx))
    print('ny: {}'.format(ny))
    print('nz: {}'.format(nz))

    x_lin=np.linspace(x_min,x_max,nx,endpoint=True)
    y_lin=np.linspace(y_min,y_max,ny,endpoint=True)
    z_lin=np.linspace(z_min,z_max,nz,endpoint=True)
    
    print("x_lin shape: {} ".format(x_lin.shape))
    print("y lin shape: {}".format(y_lin.shape))
    print("z lin shape: {}".format(z_lin.shape))

    tidy_dens=tidy_3d_field(x,y,z,dens,dx,dx,dx)
    tidy_bx  =tidy_3d_field(x,y,z,bx,dx,dx,dx) 
    tidy_by  =tidy_3d_field(x,y,z,by,dx,dx,dx) 
    tidy_bz  =tidy_3d_field(x,y,z,bz,dx,dx,dx)

    # points in cm in the left
    ypoints_l=np.linspace(y_min+10*dx,y_max-10*dx,num=nlines)
    zpoints_l=np.linspace(z_min+10*dx,z_max-10*dx,num=nlines)
    xx_l=np.ones(nlines**2)*(x_min+3*dx)
     
    yy_l,zz_l=np.meshgrid(ypoints_l,zpoints_l,indexing='ij')
    # points in cm in the right 
    ypoints_r=np.linspace(y_min+10*dx,y_max-10*dx,num=nlines)
    zpoints_r=np.linspace(z_min+10*dx,z_max-10*dx,num=nlines)
    xx_r=np.ones(nlines**2)*(x_max-3*dx)
    yy_r,zz_r=np.meshgrid(ypoints_r,zpoints_r,indexing='ij')

    xx=np.append(xx_l,xx_r)
    yy=np.append(yy_l.reshape(nlines**2),yy_r.reshape(nlines**2))
    zz=np.append(zz_l.reshape(nlines**2),zz_r.reshape(nlines**2))

    points=np.stack((xx,yy,zz),axis=-1)

    lineas_flujo(xg=x_lin,yg=y_lin,zg=z_lin,vx_in=tidy_bx,vy_in=tidy_by,vz_in=tidy_bz,start_points=points,outfile=out+out_name+str(i).zfill(3)+'.h5')

