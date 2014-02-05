import glob, os
from pylab import*
from os.path import expanduser, join
from mayavi import mlab
mlab.close(all=True)


# Data Types ------------------------------------------------------------------
geometricHeader = dtype([("nCores", float),("origo",float,(3,)),
                     ("xPoints", float),("xLimits",float,(2,)),
                     ("yPoints", float),("yLimits",float,(2,)),
                     ("zPoints", float),("zLimits",float,(2,)),
                     ])
                     
systemHeader = dtype([("id", float),("charge", float), 
                     ("position", float, (3,))])

dataType = dtype([("density", float)])

#Read files--------------------------------------------------------------------
def readFiles():
    rawDataPath = "/home/milad/kurs/qmd/density"
    stateFiles = glob.glob1(rawDataPath,'*.bin')
    nStateFiles = len(stateFiles)
    print "# state files: ", nStateFiles

    rawData = [None] * nStateFiles
    atomList= [None] * nStateFiles
    geometricData = [None] * nStateFiles
    
    for i in range(0, nStateFiles):
        cubeFile = join(rawDataPath, "cubeFile" + "%04d" % i  + ".bin")
        if os.path.exists(cubeFile):
            header, atoms, densityData = loadCubeFile(cubeFile)
            
            rawData[i]  = densityData
            atomList[i] = atoms
            geometricData[i] = header
    
    return geometricData, atomList,rawData

# File loader------------------------------------------------------------------
def loadCubeFile(fileName):
    fileName = expanduser(fileName)
    cubeFile = open(fileName, "rb")
    
    #Header with geometric data
    geometricData = fromfile(cubeFile, dtype=geometricHeader, count = 1)
    nCores = int(geometricData[0][0])
    nX = int(geometricData[0][2])
    nY = int(geometricData[0][4])
    nZ = int(geometricData[0][6])
        
    #Header with core data
    atoms = fromfile(cubeFile, dtype=systemHeader, count = nCores)

    #Density data
    densityData = fromfile(cubeFile, dtype=dataType)["density"]
    densityData = array(densityData).reshape(nX,nY,nZ)
    cubeFile.close()
    
    return geometricData, atoms, densityData
    
# Display cores----------------------------------------------------------------
def displayCores(atoms):
    nCores = len(atoms) 
    corePositions = zeros((nCores,3))
    coreCharge = zeros(nCores)
    
    for i in range(0,nCores):
        corePositions[i] =  array((atoms[i][2][1], atoms[i][2][0], atoms[i][2][2]))
        coreCharge[i]    =  sqrt(atoms[i][0])
        
    corePlot = mlab.points3d(corePositions[:,0],corePositions[:,1],corePositions[:,2],
                             coreCharge,
                             scale_factor=0.25,
                             resolution=20,
                             colormap = "prism",
                             opacity = 0.8)
        
    return corePlot.mlab_source
    
# Display density--------------------------------------------------------------
def displayDensity(geometricData, densityData): 
    xLim = array((geometricData[0][3]))
    yLim = array((geometricData[0][5]))
    zLim = array((geometricData[0][7]))
    nX = complex(geometricData[0][2])
    nY = complex(geometricData[0][4])
    nZ = complex(geometricData[0][6])

    # Display the electron localization function   
    x, y, z = mgrid[xLim[0]:xLim[1]:nX, yLim[0]:yLim[1]:nY, zLim[0]:zLim[1]:nZ]
    source = mlab.pipeline.scalar_field(x, y, z, densityData)
    vmax=densityData.max()
    vmin= densityData.min()

    densityPlot = mlab.pipeline.volume(source,
                                vmin=vmin,
                                vmax=vmax*0.05)
    densityPlot = setColormap(densityPlot)
    
    return densityPlot.mlab_source


# Save figure------------------------------------------------------------------
def saveFigure(i):
    rawDataPath = "/home/milad/kurs/qmd/density"
    figure = join(rawDataPath, "fig/density" + "%04d" % i  + ".png")
    mlab.savefig(figure)
    
    
# Color map--------------------------------------------------------------------
def setColormap(densityPlot):
    from tvtk.util.ctf import ColorTransferFunction
    ctf = ColorTransferFunction() 
   
   # Add points to CTF 
    ctf.add_rgb_point(0, 1.0, 1.0, 1.0) 
    ctf.add_rgb_point(1, 0.0, 0.0, 1.0) 
    
    densityPlot._volume_property.set_color(ctf) 
    densityPlot._ctf = ctf 
    densityPlot.update_ctf = True
    
    return densityPlot

#Animater---------------------------------------------------------------------
@mlab.animate(delay=10)
def anim(atomList, rawData, coreSource, densitySource,
         render = True):
     
     f = mlab.gcf()
     i = 0
     
     while 1:
         #Update core positions
         if coreSource!= None:
             nStep = len(atomList)
             atoms = atomList[i]
             nCores = len(atoms)
             corePositions = zeros((nCores,3))
             coreCharge = zeros(nCores)
 
             for j in range(0,nCores):
                 corePositions[j] =  array((atoms[j][2][1], atoms[j][2][0],atoms[j][2][1]))
                 coreCharge[j]    =  atoms[j][0] 
                 
             coreSource.set(x = corePositions[:,0],y = corePositions[:,1], z = corePositions[:,2])
         
         #Update density
         if densitySource!= None:
             nStep = len(rawData)
             density = rawData[i]
             densitySource.set(scalars=density)
         
         f.scene.camera.azimuth(1)
         f.scene.render()
         i+=1
         
         if render:
             saveFigure(i)
         if i >= nStep:
             i = 0
         yield
#------------------------------------------------------------------------------



densityFig = mlab.figure("charge density", bgcolor=(0, 0, 0), size=(350, 350))
geometricData, atomList, rawData = readFiles()

showCores = True
showDensity = True
animate = False



if showCores:
    coreSource = displayCores(atomList[0])
else:
    coreSource = None

if showDensity:
    densitySource = displayDensity(geometricData[0],rawData[0])
else:
    densitySource = None

if animate:
    anim(atomList, rawData, coreSource, densitySource) 











