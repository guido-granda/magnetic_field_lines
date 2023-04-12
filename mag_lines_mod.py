import numpy as np
import h5py
from scipy.interpolate import RegularGridInterpolator

def lineas_flujo( xg, yg, zg, vx_in, vy_in, vz_in,start_points=None,outfile='lines_out.hdf5'):
    """Funtion to obtain stream lines. Input xg, yg, zg coordinates arrays with "shapes m1 , m2 , m3, vx_in,vy_in,vz_in 
    vector field arrays with dimensions (m1,m2,m3), start_points logical or numpy array with dimensins (n,3) where n is the number of starting points, if it is None, it indicates if the start
    point is the center of the computational domain provided otherwise provide an array with the coordinates of the points with dimensions
    (npoints,3) where npoints is the number of points, lint  interpolation lenght. """
                    
    #https://docs.scipy.org/doc/scipy6/reference/generated/scipy.interpolate.RegularGridInterpolator.html
    # Define parallel interpolation functions 
    interpol_vx_p = RegularGridInterpolator((xg, yg, zg),vx_in,bounds_error=False, fill_value=0.)
    interpol_vy_p = RegularGridInterpolator((xg, yg, zg),vy_in,bounds_error=False, fill_value=0.)
    interpol_vz_p = RegularGridInterpolator((xg, yg, zg),vz_in,bounds_error=False, fill_value=0.)
    # Define antiparallel interpolation functions
    interpol_vx_ap = RegularGridInterpolator((xg, yg, zg),-vx_in,bounds_error=False, fill_value=0.)
    interpol_vy_ap = RegularGridInterpolator((xg, yg, zg),-vy_in,bounds_error=False, fill_value=0.)
    interpol_vz_ap = RegularGridInterpolator((xg, yg, zg),-vz_in,bounds_error=False, fill_value=0.)

    nx = len(xg)
    ny = len(yg)
    nz = len(zg)

    ## getting domain of the cube
    xmin=np.amin(xg)
    xmax=np.amax(xg)

    ymin=np.amin(yg)
    ymax=np.amax(yg)

    zmin=np.amin(zg)
    zmax=np.amax(zg)

    ds = xg[1]-xg[0]

    # Starting points for the lines 
    if start_points is not None:
        xs = start_points[:,0]
        ys = start_points[:,1]
        zs = start_points[:,2]

    else:
        xs = np.array([xg[np.int(nx/2.)]])
        ys = np.array([yg[np.int(ny/2.)]])
        zs = np.array([zg[np.int(nz/2.)]])
    #number of initial points for the lines
    npoints=len(xs)
    # output file for the lines
    f=h5py.File(outfile,'w')
    f.create_dataset('nlines',data=npoints,dtype='i4')
    
    for i in range(npoints):
        # starting coordinates for each stream line
        xp=np.array([xs[i]])
        yp=np.array([ys[i]])
        zp=np.array([zs[i]])
        # starting fields for each stream line
        vx_ar=interpol_vx_p( np.array([xp[0],yp[0],zp[0]]) )
        vy_ar=interpol_vy_p( np.array([xp[0],yp[0],zp[0]]) )
        vz_ar=interpol_vz_p( np.array([xp[0],yp[0],zp[0]]) )

        # parallel lines
        while ( (((xmin+ds)<xp[-1]) and (xp[-1]<(xmax-ds))) and (((ymin+ds)<yp[-1]) and (yp[-1]<(ymax-ds))) and (((zmin+ds)<zp[-1]) and (zp[-1]<(zmax-ds))) ):
            # parallel interpolation 
            vx_p = interpol_vx_p( np.array([xp[-1],yp[-1],zp[-1]]) )
            vy_p = interpol_vy_p( np.array([xp[-1],yp[-1],zp[-1]]) )
            vz_p = interpol_vz_p( np.array([xp[-1],yp[-1],zp[-1]]) )
        
            vel_p = np.sqrt( vx_p**2 + vy_p**2 + vz_p**2 ) + 1e-80
            dt_p  = 0.5 * ds/vel_p
            # adding forward points
            xp = np.insert( xp, len(xp), xp[-1]+vx_p*dt_p )
            yp = np.insert( yp, len(yp), yp[-1]+vy_p*dt_p )
            zp = np.insert( zp, len(zp), zp[-1]+vz_p*dt_p )
            # storing interpolated fields
            #forward field components
            vx_ar=np.insert(vx_ar,len(vx_ar),vx_p)
            vy_ar=np.insert(vy_ar,len(vy_ar),vy_p)
            vz_ar=np.insert(vz_ar,len(vz_ar),vz_p)
      
        # antiparallel lines    
        while ( (((xmin+ds)<xp[0]) and (xp[0]<(xmax-ds))) and (((ymin+ds)<yp[0]) and (yp[0]<(ymax-ds))) and (((zmin+ds)<zp[0]) and (zp[0]<(zmax-ds))) ):

            # antiparallel interpolation 
            vx_ap = interpol_vx_ap(np.array([xp[0],yp[0],zp[0]]))
            vy_ap = interpol_vy_ap(np.array([xp[0],yp[0],zp[0]]))
            vz_ap = interpol_vz_ap(np.array([xp[0],yp[0],zp[0]]))

            vel_ap = np.sqrt( vx_ap**2 + vy_ap**2 + vz_ap**2 ) + 1e-80
            dt_ap  = 0.5 * ds/vel_ap
            # adding backward points
            xp = np.insert( xp,0,xp[0]+vx_ap*dt_ap )
            yp = np.insert( yp,0,yp[0]+vy_ap*dt_ap )
            zp = np.insert( zp,0,zp[0]+vz_ap*dt_ap )
            #backward field components
            vx_ar=np.insert(vx_ar,0,-vx_ap)
            vy_ar=np.insert(vy_ar,0,-vy_ap)
            vz_ar=np.insert(vz_ar,0,-vz_ap)

        grp=f.create_group('line_'+str(i))
        
        d1=grp.create_dataset('x',data=xp)
        d1.attrs['units']='cm'
        
        d2=grp.create_dataset('y',data=yp)
        d2.attrs['units']='cm'
        
        d3=grp.create_dataset('z',data=zp)
        d3.attrs['units']='cm'
        
        d4=grp.create_dataset('x point',data=xs)
        d4.attrs['units']='cm'

        d5=grp.create_dataset('y point',data=ys)
        d5.attrs['units']='cm'
        
        d6=grp.create_dataset('z point',data=zs)
        d6.attrs['units']='cm'

        d7=grp.create_dataset('Bx',data=vx_ar)
        d7.attrs['units']='uG'
        
        d8=grp.create_dataset('By',data=vy_ar)
        d8.attrs['units']='uG'

        d9=grp.create_dataset('Bz',data=vz_ar)
        d9.attrs['units']='uG'
      
        d10=grp.create_dataset('B',data=np.sqrt(vx_ar**2+vy_ar**2+vz_ar**2))    
        d10.attrs['units']='uG'

    f.close()

def lineas_flujo_2d( xg, yg, vx_in, vy_in,start_points=None,outfile='lines_out.hdf5',max_interpol=2000):
    """Funtion to obtain stream lines. Input xg, yg coordinates arrays with "shapes m1 , m2, vx_in,vy_in 
    vector field arrays with dimensions (m1,m2), start_points logical or numpy array, if it is None, it indicates if the start
    point is the center of the computational domain provided otherwise provide an array with the coordinates of the points with dimensions
    (npoints,2) where npoints is the number of points, lint  interpolation lenght. """
                    
    #https://docs.scipy.org/doc/scipy6/reference/generated/scipy.interpolate.RegularGridInterpolator.html
    # Define las funciones de interpolaciónparalelas
    interpol_vx_p = RegularGridInterpolator((xg, yg),vx_in,bounds_error=False, fill_value=0.)
    interpol_vy_p = RegularGridInterpolator((xg, yg),vy_in,bounds_error=False, fill_value=0.)
    # Define las funciones de interpolació antiparalelas
    interpol_vx_ap = RegularGridInterpolator((xg, yg),-vx_in,bounds_error=False, fill_value=0.)
    interpol_vy_ap = RegularGridInterpolator((xg, yg),-vy_in,bounds_error=False, fill_value=0.)

    nx = len(xg)
    ny = len(yg)

    ## getting domain of the cube
    xmin=np.amin(xg)
    xmax=np.amax(xg)

    ymin=np.amin(yg)
    ymax=np.amax(yg)

    ds = xg[1]-xg[0]

    # Puntos iniciales
    if start_points is not None:
        xs = start_points[:,0]
        ys = start_points[:,1]

    else:
        xs = np.array([xg[np.int(nx/2.)]])
        ys = np.array([yg[np.int(ny/2.)]])
    #number of initial points for the lines
    npoints=len(xs)
    #  output file for the lines  
    f=h5py.File(outfile,'w')
    f.create_dataset('nlines',data=npoints,dtype='i4')
    for i in range(npoints):
        print('line {}'.format(i))
        # starting coordinates for each stream line
        xp=np.array([xs[i]])
        yp=np.array([ys[i]])
        # starting fields for each stream line
        vx_ar=interpol_vx_p( np.array([xp[0],yp[0]]) )
        vy_ar=interpol_vy_p( np.array([xp[0],yp[0]]) )
        # parallel lines
        j=0
        while ( (((xmin+ds)<xp[-1]) and (xp[-1]<(xmax-ds))) and (((ymin+ds)<yp[-1]) and (yp[-1]<(ymax-ds))) and (j<=max_interpol)):
            j+= 1
            # parallel interpolation 
            vx_p = interpol_vx_p( np.array([xp[-1],yp[-1]]) )
            vy_p = interpol_vy_p( np.array([xp[-1],yp[-1]]) )
        
            vel_p = np.sqrt( vx_p**2 + vy_p**2) + 1e-80
            dt_p  = 0.5 * ds/vel_p
        
            # adding forward points
            xp = np.insert( xp, len(xp), xp[-1]+vx_p*dt_p )
            yp = np.insert( yp, len(yp), yp[-1]+vy_p*dt_p )
            
            # storing interpolated fields
            vx_ar=np.insert(vx_ar,len(vx_ar),vx_p)
            vy_ar=np.insert(vy_ar,len(vy_ar),vy_p)

        # antiparallel lines
        while ( (((xmin+ds)<xp[0]) and (xp[0]<(xmax-ds))) and (((ymin+ds)<yp[0]) and (yp[0]<(ymax-ds))) and (j<=max_interpol)):
            j+=1
            # antiparallel interpolation 
            vx_ap = interpol_vx_ap(np.array([xp[0],yp[0]]))
            vy_ap = interpol_vy_ap(np.array([xp[0],yp[0]]))

            vel_ap = np.sqrt(vx_ap**2 + vy_ap**2) + 1e-80
            dt_ap  = 0.5 * ds/vel_ap
            # adding backward points
            xp = np.insert( xp,0,xp[0]+vx_ap*dt_ap )
            yp = np.insert( yp,0,yp[0]+vy_ap*dt_ap )
            #backward field components
            vx_ar=np.insert(vx_ar,0,-vx_ap)
            vy_ar=np.insert(vy_ar,0,-vy_ap)

        grp=f.create_group('line_'+str(i))
        
        d1=grp.create_dataset('x',data=xp)
        d1.attrs['units']='cm'
        
        d2=grp.create_dataset('y',data=yp)
        d2.attrs['units']='cm'
        
        d4=grp.create_dataset('x point',data=xs)
        d4.attrs['units']='cm'

        d5=grp.create_dataset('y point',data=ys)
        d5.attrs['units']='cm'
        
        d7=grp.create_dataset('Bx',data=vx_ar)
        d7.attrs['units']='uG'
        
        d8=grp.create_dataset('By',data=vy_ar)
        d8.attrs['units']='uG'

        d10=grp.create_dataset('B',data=np.sqrt(vx_ar**2+vy_ar**2))    
        d10.attrs['units']='uG'

    f.close()

