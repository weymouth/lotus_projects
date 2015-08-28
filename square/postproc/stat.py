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

# get vorticity
calc1 = PythonCalculator()
calc1.Expression = 'curl(Velocity)[:,2]*100'
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
view.CameraPosition = [400, 0, 1000]
view.CameraFocalPoint = [400, 0, 0]
#view.ResetCamera()
#camera = GetActiveCamera()
#camera.Dolly(2)

# get image
animationScene1.GoToLast()
#animationScene1.AnimationTime=220
SaveScreenshot('out1.png')
#WriteAnimation('omega.avi', FrameRate=15.0, Compression=True)

# read in mixing data
if os.path.exists('fkepr.vtr.pvd'):
  mix = PVDReader( FileName='fkepr.vtr.pvd' )
  print "read fkepr.vtr"

  # get mix kinetic energy display
  mixDisplay = Show()
  ColorBy(mixDisplay, ('POINTS', 'Pressure'))
  pressureLUT = GetColorTransferFunction('Pressure')
  pressureLUT.RGBPoints= [0, 0, 0, 0,
                          0.33, 1, 0, 0,
                          0.66, 1, 1, 0,
                          1, 1, 1, 1]
  pressureLUT.ColorSpace = 'RGB'
  mixDisplay.RescaleTransferFunctionToDataRange()
#  mixDisplay.SetScalarBarVisibility ( view , True)
#  pressureLUT.RescaleTransferFunction(0, .25)
#  pressureLUT.MapControlPointsToLogSpace()
#  pressureLUT.UseLogScale = 1

  # get image
  SaveScreenshot('out2.png')
#  WriteAnimation('fke.avi', FrameRate=15.0, Compression=True)
