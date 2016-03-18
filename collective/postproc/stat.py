#### import the simple module from the paraview
from paraview.simple import *

# read in body data
body = PVDReader( FileName='bodyF.vti.pvd' )
print "read body.vti"

# get animation scene
animationScene1 = GetAnimationScene()
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'Reflect'
reflectp5 = Reflect(Input=body)
reflectp5.Plane = 'Y Min'
reflect1 = Reflect(Input=reflectp5)
reflect1.Plane = 'Y Min'

# create a new 'Contour'
contour1 = Contour(Input=reflect1)
contour1.ContourBy = ['POINTS', 'Pressure']
contour1.Isosurfaces = [0.0]
contour1.PointMergeMethod = 'Uniform Binning'

# change solid color
contour1Display = Show()
contour1Display.DiffuseColor = [0.7, 0.87, 0.55]

# read in vorticity data
vort = PVDReader( FileName='vortF.vtr.pvd' )
print "read vort.vtr"

# create a new 'Reflect'
reflect1p5 = Reflect(Input=vort)
reflect1p5.Plane = 'Y Min'
reflect2 = Reflect(Input=reflect1p5)
reflect2.Plane = 'Y Min'

# create a new 'Contour'
contour2 = Contour(Input=reflect2)
contour2.ContourBy = ['POINTS', 'Pressure']
contour2.Isosurfaces = [-0.0005]
contour2.PointMergeMethod = 'Uniform Binning'

contour2Display = Show()
contour2Display.DiffuseColor = [0.3, 0.6, 0.9]
contour2Display.Ambient = 0.2

# get view
view = GetRenderView()
view.Background = [1.0, 1.0, 1.0]
view.ViewSize = [1280,640]
view.CameraPosition = [170, -64, -700]
view.CameraFocalPoint = [170,-64, -10]
view.CameraParallelScale = 400
#view.ResetCamera()
#camera = GetActiveCamera()
#camera.Dolly(2)

# get single image for report
#animationScene1.GoToLast()
#SaveScreenshot('out1.png')

# get animation images
WriteAnimation('lambda2/lambda2.png',FrameRate=15.0, Compression=True)
