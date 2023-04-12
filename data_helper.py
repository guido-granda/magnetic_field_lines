import sys
import numpy as np

def find_nearest(r_ar,target,res):
    """Function to find the nearest value inside the r_ar array to the target value."""
    idx=((res>=(r_ar-target)) & ((r_ar-target)>=-res))
    new_ar=r_ar[idx]
    return new_ar,idx
def tidy_image_simplified(ar1,dar1,n1,n2,ar1_o,prop_arr):
    """Code to order scalar and vector fields belonging to an slice. These fields are 
    coming out from yt.save_as_dataset in a strange order. Imputs: ar1 array containing the coordinate that goes
    on the Y axis of the slice, dar1 resolution on the Y axis, n1 the dimension of the ar1 array, n2 the dimension
    of the array that goes on the X axis of the slice, ar1_o increasing ordered version of ar1, prop_arr field
    to reorder. It returns the ordered version of prop_arr as out_arr"""
    out_arr=np.zeros((n1,n2))
    for i in range(n1):
        id1=(((ar1_o[i]-dar1) <= ar1) & (ar1 <= (ar1_o[i]+dar1)))
        out_arr[i,:]=prop_arr[id1]
    return out_arr
def tidy_3d_field(xarr,yarr,zarr,prop_arr,dx,dy,dz):
    """Code to tidy up an array coming from yt.save as dataset. Inputs: xarr, yarr, and zarr coordinate arrays, prop_arr
    is the untidy property array, dx, dy ,and dz are the resolution along each axis. It calls find nearest and tidy_image routines"""
    # x array
    xmin=np.amin(xarr)
    xmax=np.amax(xarr)
    diff_x=(xmax-xmin)/dx+1
    nx=int(round(diff_x))    
    # y array
    ymin=np.amin(yarr)
    ymax=np.amax(yarr)
    diff_y=(ymax-ymin)/dy+1
    ny=int(round(diff_y))
    # z array
    zmin=np.amin(zarr)
    zmax=np.amax(zarr)
    diff_z=(zmax-zmin)/dz+1
    nz=int(round(diff_z))
    # ordered version of the y array
    y_ar_o=np.linspace(ymin,ymax,ny,endpoint=True)
    prop_arr_tidy=np.zeros((nx,ny,nz))
    for i in range(nx):
        xs_i,id_i=find_nearest(xarr,xmin+i*dx,dx/2)
        ys_i=yarr[id_i] # y coordinates that correspond to xs_i
        zs_i=zarr[id_i] # z coordinates that correspond to xs_i
        prop_arr_i=prop_arr[id_i] # property array that correspond to these coordinates
        prop_arr_tidy[i,:,:]=tidy_image_simplified(ys_i,dy/2,ny,nz,y_ar_o,prop_arr_i)
    return prop_arr_tidy
def tidy_2d_field(xarr,yarr,prop_arr,dx,dy):
    """Code to tidy up an array coming from yt.save as dataset. Inputs: xarr, yarr, and zarr coordinate arrays, prop_arr
    is the untidy property array, dx, dy ,and dz are the resolution along each axis. It calls find nearest and tidy_image routines"""
    # x array
    xmin=np.amin(xarr)
    xmax=np.amax(xarr)
    diff_x=(xmax-xmin)/dx+1
    nx=int(round(diff_x))    
    # y array
    ymin=np.amin(yarr)
    ymax=np.amax(yarr)
    diff_y=(ymax-ymin)/dy+1
    ny=int(round(diff_y))
    # ordered version of the y array
    y_ar_o=np.linspace(ymin,ymax,ny,endpoint=True)

    prop_arr_tidy=np.zeros((nx,ny))
    for i in range(nx):
        xs_i,id_i=find_nearest(xarr,xmin+i*dx,dx/2)
        ys_i=yarr[id_i] # y coordinates that correspond to xs_i
        #print('ys_i shape: {}'.format(ys_i.shape))
        prop_arr_i=prop_arr[id_i] # property array that correspond to these coordinates
        #print('prop_arr shape: {}'.format(prop_arr_i.shape)) 
        prop_arr_tidy[i,:]=tidy_image_simplified_2d(ys_i,dy/2,ny,y_ar_o,prop_arr_i)
    return prop_arr_tidy

