import nibabel
import numpy
from sklearn.utils.extmath import cartesian
from scipy import stats

# Simulate some data
def data(dims=(8,8,8,100),nsubj=20):
    nx, ny, nz, nt = dims
    x = numpy.linspace(-(nx-1)/2.,(nx-1)/2., nx)
    y = numpy.linspace(-(ny-1)/2.,(ny-1)/2., ny)
    z = numpy.linspace(-(nz-1)/2.,(nz-1)/2., nz)
    ijk = cartesian([x,y,z])
    ijk_p = stats.norm.sf(abs(ijk))
    p = numpy.mean(ijk_p, axis=1).reshape(nx,ny,nz,1)
    affine = numpy.diag([3,3,3,1])

    X = []
    for i in range(nsubj):
        r = numpy.random.random([nx,ny,nz,nt])
        array_noise = numpy.random.normal(loc=0., scale=1., size=[nx,ny,nz,nt])
        array_signal = r < p
        array_data = numpy.zeros_like(array_noise)
        array_data[...,0:nt/2] = array_noise[...,0:nt/2] + array_signal[...,nt/2:nt]
        array_data[...,nt/2:nt] = array_noise[...,nt/2:nt]
        array_img = nibabel.Nifti1Image(array_data, affine)
        X.append(array_img)
        #nibabel.save(array_img, 'sim{i:02d}.nii'.format(i=i))

    return X

# Generate the CV blocks
def cv(ncv=10, nt=100):
    r = numpy.ceil(nt/float(ncv))
    CV = numpy.tile(numpy.identity(ncv),(r,1))[0:nt,] == 1
    return CV

def targets(nt=100):
    y = numpy.zeros(nt, dtype=numpy.bool)
    y[0:nt/2] = True
    return y
