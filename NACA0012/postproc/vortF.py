#### import the simple module from the paraview
from paraview.simple import *

# read in fluid data
fluid = PVDReader( FileName='vortF.vti.pvd' )

# get animation scene
animationScene1 = GetAnimationScene()
animationScene1.UpdateAnimationUsingDataTimeSteps()
animationScene1.GoToLast()

# get sqrt(-lam2)
calc1 = Calculator()
calc1.ResultArrayName = 'Lam2'
calc1.Function = 'sqrt(max(0,-Pressure))'

# change representation type
fluidDisplay = Show()
fluidDisplay.SetRepresentationType('Volume')
ColorBy(fluidDisplay, ('POINTS', 'Lam2'))

lam2LUT = GetColorTransferFunction('Lam2')
lam2LUT.RGBPoints = [0, 0, 0, 0, 
                     0.33, 1, 0, 0, 
                     0.66, 1, 1, 0, 
                     1, 1, 1, 1]

lam2PWF = GetOpacityTransferFunction('Lam2')
lam2PWF.Points = [0, 0, 0.5, 0, 
                  0.1, 0, 0.5, 0, 
                  1, 1, 0.5, 0]

lam2LUT.RescaleTransferFunction(0.0, 0.25)
lam2PWF.RescaleTransferFunction(0.0, 0.25)

# create a new 'PVD Reader'
body = PVDReader(FileName='bodyF.vti.pvd')
clip1 = Clip()
clip1.ClipType = 'Scalar'
clip1.Scalars = ['POINTS', 'Pressure']
clip1.Value = 0.0
clip1.InsideOut = 1
clip1Display = Show()
ColorBy(clip1Display, None)
clip1Display.DiffuseColor = [1.0, 1.0, 0.0]

# get view
view = GetRenderView()
view.ViewSize = [1280,640]
view.ResetCamera()
camera = GetActiveCamera()
camera.Dolly(2)
camera.Elevation(20)
camera.Pitch(-5)

# write out
SaveScreenshot('vortF.png')
WriteAnimation('vortF.avi', FrameRate=15.0, Compression=True)
