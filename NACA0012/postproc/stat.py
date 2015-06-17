#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# read in fluid data
fluid = PVDReader( FileName='fluid.vti.pvd' )

# get animation scene
animationScene1 = GetAnimationScene()
animationScene1.GoToLast()

# get vorticity
calc1 = PythonCalculator()
calc1.Expression = 'curl(Velocity)[:,2]*325'
calc1.ArrayName = 'Omega'

# get fluid display
fluidDisplay = Show()
ColorBy(fluidDisplay, ('POINTS', 'Omega'))
omegaLUT = GetColorTransferFunction('Omega')
omegaLUT.RescaleTransferFunction(-50.0, 50.0)

# read in body data
body = PVDReader( FileName='bodyF.vti.pvd' )

# get contours
Contour1 = Contour()
Contour1.ContourBy = ['POINTS', 'Pressure']
Contour1.Isosurfaces = [-2.0, -1.0, 0.0, 1.0, 2.0]
contourRep = GetDisplayProperties(Contour1)
contourRep.LineWidth = 3
Show()

# get view
view = GetRenderView()
view.ViewSize = [1280,640]
#view.ResetCamera()
view.CameraFocalPoint = [175, 20, 0]
view.CameraPosition = [175, 20, 175]

ExportView('out.pdf')

