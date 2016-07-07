#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
import glob
import re

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

outfile = '/home/alex/Videos/keyfrtest.ogv'
func = 'e'
# anim_duration = 40.0
anim_duration = 5.0
pause_after_duration = 5.0
total_duration = anim_duration + pause_after_duration
fileNames = sorted(glob.glob('/home/alex/code/ngsapps/precip_gauss_40/*.vtk'), key=numericalSort)
# print(fileNames)
# create a new 'Legacy VTK Reader'
reader = LegacyVTKReader(FileNames=fileNames)
framerate = len(reader.TimestepValues) / anim_duration
print 'Resulting framerate: ' + framerate

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [1920, 1080]
# change interaction mode for render view
renderView1.InteractionMode = '3D'

renderView1.OrientationAxesVisibility = 0

if func == 'e':
    # create a new 'Warp By Scalar'
    warpByScalar1 = WarpByScalar(Input=reader)
    warpByScalar1.Scalars = ['POINTS', func]
    warpByScalar1.ScaleFactor = 10.0

    # show data in view
    display = Show(warpByScalar1, renderView1)
else:
    # show data in view
    display = Show(reader, renderView1)

# set scalar coloring
ColorBy(display, ('POINTS', func))

# # show color bar/color legend
# display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for func
LUT = GetColorTransferFunction(func)
# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
LUT.ApplyPreset('jet', True)

# get opacity transfer function/opacity map for func
PWF = GetOpacityTransferFunction(func)

# Rescale transfer function
if func == 'e':
    PWF.RescaleTransferFunction(0.0, 1.0)
    LUT.RescaleTransferFunction(0.0, 1.0)
else:
    PWF.RescaleTransferFunction(-0.4, 0.1)
    LUT.RescaleTransferFunction(-0.4, 0.1)

# # get color legend/bar for LUT in view renderView1
# LUTColorBar = GetScalarBar(LUT, renderView1)
# # Properties modified on LUTColorBar
# LUTColorBar.DrawSubTickMarks = 0
# LUTColorBar.RangeLabelFormat = '%.1f'
# LUTColorBar.AutomaticLabelFormat = 0
# LUTColorBar.LabelFormat = '%-#6.2g'

renderView1.AxesGrid.Visibility = 1
# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.FacesToRender = 36
renderView1.AxesGrid.XTitle = ''
renderView1.AxesGrid.YTitle = ''
renderView1.AxesGrid.ZTitle = ''

# current camera placement for renderView1
if func == 'e':
    renderView1.CameraPosition = [200, -50, 450]
    renderView1.CameraFocalPoint = [200, 100, -20]
else:
    renderView1.CameraPosition = [100, 100, 400]
    renderView1.CameraFocalPoint = [100, 100, 0]

# cam = GetActiveCamera()
# cam.Elevation(-10)
# print cam.GetPosition()
# print cam.GetViewUp()
# cam.Azimuth(-10)

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

keyf0 = CompositeKeyFrame()
keyf0.Interpolation = 'Ramp'
keyf0.KeyTime = 0
keyf0.KeyValues = [0]

keyf1 = CompositeKeyFrame()
keyf1.Interpolation = 'Ramp'
keyf1.KeyTime = anim_duration / total_duration
keyf1.KeyValues = [2000]

keyf2 = CompositeKeyFrame()
keyf2.Interpolation = 'Ramp'
keyf2.KeyTime = 1
keyf2.KeyValues = [2000]

tt = GetTimeTrack()
tt.KeyFrames = [keyf0, keyf1, keyf2]
tt.UseAnimationTime = 0

# # Properties modified on animationScene1
# animationScene1.PlayMode = 'Sequence'

# # Properties modified on animationScene1
# animationScene1.NumberOfFrames = int(total_duration * frame_rate)

# # save animation images/movie
# WriteAnimation(outfile, Magnification=1, FrameRate=frame_rate, Compression=False)

# # Properties modified on animationScene1
# animationScene1.PlayMode = 'Real Time'

# # Properties modified on animationScene1
# animationScene1.Duration = int(total_duration)

# animationScene1.Play()
