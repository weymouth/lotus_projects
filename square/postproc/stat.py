#### import the simple module from the paraview
from paraview.simple import *

# read in fluid data
import os.path

fluid = PVDReader( FileName='fluid.vtr.pvd' )
print "read fluid.vtr"

# get animation scene
animationScene1 = GetAnimationScene()
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get vorticity
calc1 = PythonCalculator()
calc1.Expression = 'curl(Velocity)[:,2]*50'
calc1.ArrayName = 'Omega'

# get fluid display
fluidDisplay = Show()
ColorBy(fluidDisplay, ('POINTS', 'Omega'))
omegaLUT = GetColorTransferFunction('Omega')
omegaLUT.RGBPoints= [-4, 1, 0, 0,
                      4, 0, 0, 1]
#omegaLUT.ColorSpace = 'RGB'
fluidDisplay.Ambient = 0.2

# get view
view = GetRenderView()
view.ViewSize = [1280,640]
view.CameraPosition = [100, 0, 250]
view.CameraFocalPoint = [100, 0, 0]
#view.ResetCamera()
#camera = GetActiveCamera()
#camera.Dolly(2)

# create a new 'PVD Reader'
bodyFvtipvd = PVDReader(FileName='bodyF.vtr.pvd')
print "read body.vtr"

# show data in view
bodyDisplay = Show()
ColorBy(bodyDisplay, ('POINTS', 'Pressure'))
pressureLUT = GetColorTransferFunction('Pressure')
pressureLUT.RGBPoints = [-1, 0, 0, 0,
                          1, 0, 0, 0]

pressureLUT.EnableOpacityMapping = 1
pressurePWF = GetOpacityTransferFunction('Pressure')
pressurePWF.Points = [-1, 1, 0.5, 0,
                       1, 0, 0.5, 0]
# get image
animationScene1.AnimationTime=220
SaveScreenshot('out1.png')
#WriteAnimation('omega.avi', FrameRate=15.0, Compression=True)
