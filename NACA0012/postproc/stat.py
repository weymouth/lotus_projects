#### import the simple module from the paraview
from paraview.simple import *

# read in fluid data
import os.path
if os.path.exists('flu2d.vtr.pvd'):
  fluid = PVDReader( FileName='flu2d.vtr.pvd' )
  print "read flu2d.vtr"
else:
  fluid = PVDReader( FileName='fluid.vtr.pvd' )
  print "read fluid.vtr"
	
# get animation scene
animationScene1 = GetAnimationScene()
animationScene1.UpdateAnimationUsingDataTimeSteps()
animationScene1.GoToLast()

# get vorticity
calc1 = PythonCalculator()
calc1.Expression = 'curl(Velocity)[:,2]*164'
calc1.ArrayName = 'Omega'

# get fluid display
fluidDisplay = Show()
ColorBy(fluidDisplay, ('POINTS', 'Omega'))
omegaLUT = GetColorTransferFunction('Omega')
omegaLUT.RGBPoints= [-25, 0, 0, 1, 0, 1, 1, 1, 25, 1, 0, 0]
fluidDisplay.Ambient = 0.3

# read in body data
#body = PVDReader( FileName='bodyF.vtr.pvd' )
#bodyDisplay = Show()
#ColorBy(bodyDisplay, ('POINTS', 'Pressure'))

# Modify transfer functions
pressureLUT = GetColorTransferFunction('Pressure')
pressureLUT.RGBPoints= [-2, 0, 0, 0, 2, 0, 0, 0]
pressureLUT.EnableOpacityMapping = 1

pressurePWF = GetOpacityTransferFunction('Pressure')
pressurePWF.Points = [-2, 1, 0.5, 0, 2, 0, 0.5, 0]

# get view
view = GetRenderView()
view.ViewSize = [1280,640]
view.ResetCamera()
camera = GetActiveCamera()
camera.Dolly(4)

SaveScreenshot('out.png')
#WriteAnimation('temp.avi', FrameRate=15.0, Compression=True)

