import nibabel
import numpy
import generate_data as generate

GenerateNewData = False

if GenerateNewData:
    X = generate.data(dims=(8,8,8,100),nsubj=20)
    for i,array_img in enumerate(X):
        nibabel.save(array_img, 'sim{i:02d}.nii'.format(i=i))

else:
    X = []
    nsubj = 20
    for i in range(nsubj):
        X.append(nibabel.load('sim{i:02d}.nii'.format(i=i)))

CV = generate.cv(ncv=10, nt=100)
y = generate.targets(nt=100)


