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
calc1.Expression = 'curl(Velocity)[:,2]*160'
calc1.ArrayName = 'Omega'

# get fluid display
fluidDisplay = Show()
ColorBy(fluidDisplay, ('POINTS', 'Omega'))
omegaLUT = GetColorTransferFunction('Omega')
omegaLUT.RGBPoints= [-15, 1, 0, 0, 0, 1, 1, 1, 15, 0, 0, 1]
fluidDisplay.Ambient = 0.3

# get view
view = GetRenderView()
view.ViewSize = [1280,640]
view.CameraPosition = [200, 0, 750]
view.CameraFocalPoint = [200, 0, 0]

# get image
SaveScreenshot('out1.png')
#WriteAnimation('omega.avi', FrameRate=15.0, Compression=True)

# read in mixing data
if os.path.exists('fkepr.vtr.pvd'):
  Hide()
  mix = PVDReader( FileName='fkepr.vtr.pvd' )
  print "read fkepr.vtr"

  # get mix kinetic energy display
  mixDisplay = Show()
  ColorBy(mixDisplay, ('POINTS', 'Pressure'))
  pressureLUT = GetColorTransferFunction('Pressure')
  pressureLUT.RGBPoints= [0, 0, 0, 0, 0.33, 1, 0, 0, 0.66, 1, 1, 0, 1, 1, 1, 1]
#  mixDisplay.RescaleTransferFunctionToDataRange()
#  mixDisplay.SetScalarBarVisibility ( view , True)
  pressureLUT.RescaleTransferFunction(0, 0.138)
#  pressureLUT.MapControlPointsToLogSpace()
#  pressureLUT.UseLogScale = 1

  # get image
  SaveScreenshot('out2.png')
#  WriteAnimation('fke.avi', FrameRate=15.0, Compression=True)

# read in second invarient
if os.path.exists('invt2.vtr.pvd'):
  Hide()
  inv = PVDReader( FileName='invt2.vtr.pvd' )
  print "read fkepr.vtr"

  # get mix kinetic energy display
  invDisplay = Show()
  ColorBy(invDisplay, ('POINTS', 'Pressure'))
  pressureLUT.RGBPoints= [-1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0]
  invDisplay.Ambient = 0.3

  # get image
  SaveScreenshot('out3.png')
