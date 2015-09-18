#### import the simple module from the paraview
from paraview.simple import *

# read in fluid data
import os.path

fluid = PVDReader( FileName='fluid.vti.pvd' )
print "read fluid.vti"

# get animation scene
animationScene1 = GetAnimationScene()
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'Slice'
slice1 = Slice()
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.SliceType.Origin = [195.5, 50.5, 1.0]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]
slice1.Crinkleslice = 1

# get vorticity
calc1 = PythonCalculator()
calc1.Expression = 'curl(Velocity)[:,2]*100'
calc1.ArrayName = 'Omega'

# get fluid display
fluidDisplay = Show()
ColorBy(fluidDisplay, ('POINTS', 'Omega'))
omegaLUT = GetColorTransferFunction('Omega')
omegaLUT.RGBPoints= [-10, 1, 0, 0,
                      -1, 1, 0.5, 0,
                    -0.3, 1, 1, 0.5,
                       0, 1, 1, 1,
                     0.3, 0.5, 1, 1,
                       1, 0, 0.5, 1,
                      10, 0, 0, 1]
omegaLUT.ColorSpace = 'RGB'
fluidDisplay.Ambient = 0.2

# get view
view = GetRenderView()
view.ViewSize = [1280,640]
view.CameraPosition = [50, 25, 150]
view.CameraFocalPoint = [50, 25, 0]
#view.ResetCamera()
#camera = GetActiveCamera()
#camera.Dolly(2)

# create a new 'PVD Reader'
bodyFvtipvd = PVDReader(FileName='bodyF.vti.pvd')
print "read body.vti"

# create a new 'Contour'
contour1 = Contour(Input=bodyFvtipvd)
contour1.ContourBy = ['POINTS', 'Pressure']
contour1.Isosurfaces = [0.0]
contour1.PointMergeMethod = 'Uniform Binning'

# show data in view
contour1Display = Show()
contour1Display.ColorArrayName = [None, '']
contour1Display.DiffuseColor = [1.0, 1.0, 0.0]
contour1Display.Ambient = 0.2

# get image
animationScene1.AnimationTime=220
SaveScreenshot('out1.png')
WriteAnimation('omega.avi', FrameRate=15.0, Compression=True)
