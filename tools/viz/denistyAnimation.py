import glob, os
from pylab import*
from os.path import expanduser, join
from mayavi import mlab
import volumeSlicer as vs
mlab.close(all=True)


" Data Types******************************************************************" 
geometricHeader = dtype([("nCores", float),("origo",float,(3,)),
                     ("xPoints", float),("xLimits",float,(2,)),
                     ("yPoints", float),("yLimits",float,(2,)),
                     ("zPoints", float),("zLimits",float,(2,)),
                     ])
                     
systemHeader = dtype([("id", float),("charge", float), 
                     ("position", float, (3,))])

dataType = dtype([("density", float)])

"Read files*******************************************************************" 
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

"File loader *****************************************************************" 
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
            
"*****************************************************************************"   
class cubeVizualizer:
    """
    Superclass for vizualization of cube files, inlcuding:
        
        1. Core positions
        2. Charge density
    """ 
    def __init__(self, geometricData, atomList, densityData,
                 bgcolor, showCores, showDensity,
                 animate, saveFigure):

        self.geometricData = geometricData
        self.atomList = atomList
        self.densityData  = densityData    
        self.showCores = showCores
        self.showDensity = showDensity
        self.bgcolor = bgcolor
        self.animate = animate
        self.saveFigure = saveFigure
        
        self.vmax = self.densityData[0].max()
        self.vmin = self.densityData[0].min()
        
        self.xLim = array((self.geometricData[0][0][3]))
        self.yLim = array((self.geometricData[0][0][5]))
        self.zLim = array((self.geometricData[0][0][7]))
        self.nX = complex(self.geometricData[0][0][2])
        self.nY = complex(self.geometricData[0][0][4])
        self.nZ = complex(self.geometricData[0][0][6])

        
        self.x, self.y, self.z = mgrid[self.xLim[0]:self.xLim[1]:self.nX, 
                                       self.yLim[0]:self.yLim[1]:self.nY,
                                       self.zLim[0]:self.zLim[1]:self.nZ]
        
        print "# points: ",int(self.nX.real),int(self.nY.real),int(self.nZ.real)
        print "limits: ", self.xLim, self.yLim, self.zLim
                                       
    
    def vizualize(self):
        
        # Set background color
        if self.bgcolor == "w":
            bgcolor = (1,1,1)
        elif self.bgcolor == "b":
            bgcolor = (0,0,0)
        else:
            bgcolor = None
        
        mlab.figure("cubeViz", bgcolor= bgcolor, size=(350, 350))
        
        #Show cores?        
        if self.showCores:
            self.coreSource = self.displayCores(self.atomList[0])
        else:
            self.coreSource = None
            
        
        #Show density?    
        if self.showDensity:
            self.densitySource = self.displayDensity(self.densityData[0])
        else:
            self.densitySource = None
            
        
        #Animation?    
        if self.animate:
           self.anim(self)
        
    def displayCores(self,atoms):
        nCores = len(atoms) 
        corePositions = zeros((nCores,3))
        coreCharge = zeros(nCores)
        
        for i in range(0,nCores):
            corePositions[i] =  array((atoms[i][2][1], atoms[i][2][0], 
                                        atoms[i][2][2]))
            coreCharge[i]    =  atoms[i][0]
            
        corePlot = mlab.points3d(corePositions[:,0],corePositions[:,1],
                                 corePositions[:,2],
                                 sqrt(coreCharge),
                                 scale_factor=0.25,
                                 resolution=20,
                                 colormap = "prism",
                                 opacity = 0.8)
            
        return corePlot.mlab_source
        
    def displayDensity(self, densityData):
        raise NotImplementedError
        
        
    @mlab.animate(delay=10, ui = True)
    def anim(self): 
        f = mlab.gcf()
        i = 0;

        while 1:
         #Update core positions
         if self.coreSource!= None:
             nStep = len(self.atomList)
             atoms = self.atomList[i]
             nCores = len(atoms)
             corePositions = zeros((nCores,3))
             coreCharge = zeros(nCores)
 
             for j in range(0,nCores):
                 corePositions[j] =  array((atoms[j][2][1],
                                     atoms[j][2][0],
                                     atoms[j][2][1]))
                 coreCharge[j]    =  atoms[j][0] 
                 
             self.coreSource.set(x = corePositions[:,0],
                            y = corePositions[:,1],
                            z = corePositions[:,2])
         
         #Update density
         if self.densitySource!= None:
             nStep = len(self.densityData)
             density = self.densityData[i]
             self.densitySource.set(scalars=density)
         
#         f.scene.camera.azimuth(10)
#             f.scene.camera.elevation(10)
         f.scene.render()

         if self.saveFigure: 
             rawDataPath = "/home/milad/kurs/qmd/density"
             figure = join(rawDataPath, "fig/cubeFile" + "%04d" % i+".png")
             mlab.savefig(figure)
         
         i+=1             
         if i >= nStep:
             self.saveFigure = False
             i = 0
         yield

"*****************************************************************************"   
class contourRepresentation(cubeVizualizer):
    """
    Electron density represented by iso-surfaces 
    """ 
    
    def displayDensity(self, densityValues):
        densityPlot = mlab.contour3d(self.x, self.y, self.z, 
                                     densityValues,
                                     vmin = self.vmin,
                                     vmax = self.vmax*0.4,
                                     colormap = "Blues",
                                     opacity = 0.2,
                                     contours=50, transparent=True)
        mlab.outline()
        mlab.axes()
        return densityPlot.mlab_source
"*****************************************************************************" 
class volumeRepresentation(cubeVizualizer):
    """
    Volume rendering to display the electron density.
    """ 
    
    def displayDensity(self, densityValues):
        source = mlab.pipeline.scalar_field(self.x, self.y, self.z, 
                                            densityValues)

        densityPlot = mlab.pipeline.volume(source,
                                    vmin = self.vmin,
                                    vmax = self.vmax*0.3)
        #Set colormap                            
        densityPlot = self.setColormap(densityPlot)
    
        return densityPlot.mlab_source
        
    def setColormap(self, densityPlot):
        from tvtk.util.ctf import ColorTransferFunction
        ctf = ColorTransferFunction() 
       
       # Add points to CTF 
        ctf.add_rgb_point(0, 1.0, 1.0, 1.0) 
        ctf.add_rgb_point(0.8, 0.0, 0.0, 1.0) 
        ctf.add_rgb_point(1, 0.0, 0.0, 1.0) 
        
        densityPlot._volume_property.set_color(ctf) 
#        densityPlot._ctf = ctf 
        densityPlot.update_ctf = True
        
        return densityPlot
"*****************************************************************************" 
class slicer1Representation(cubeVizualizer):
    """
    Electron density represented by slices
    """ 
    
    def displayDensity(self, densityValues):
            source = mlab.pipeline.scalar_field(self.x, self.y, self.z, 
                                                densityValues)    
            densityPlot = mlab.pipeline.image_plane_widget(source,
                            colormap="hot",opacity = 1, transparent=True,                   
                            plane_orientation='x_axes',vmin = self.vmin,
                                    vmax = self.vmax*0.5,
                            slice_index=20,)
    
            densityPlot = mlab.pipeline.image_plane_widget(source,
                        colormap="hot",  opacity = 1,transparent=True,vmin = self.vmin,
                                    vmax = self.vmax*0.5,
                        plane_orientation='y_axes',
                        slice_index=20,)         
                        
            densityPlot = mlab.pipeline.image_plane_widget(source,
                    colormap="hot",   opacity = 1, transparent=True, vmin = self.vmin,
                                    vmax = self.vmax*0.5,                            
                    plane_orientation='z_axes',
                    slice_index=20,)    
                    
            mlab.outline()
#            mlab.axes()
            
            return densityPlot.mlab_source 

"*****************************************************************************" 
class slicer2Representation(cubeVizualizer):
    """
    Electron density represented by slices
    """ 
    
    def displayDensity(self, densityValues):
            
        m = vs.VolumeSlicer(data=densityValues)
        m.configure_traits()
        return m
        
"*****************************************************************************" 
class surf2DRepresentation(cubeVizualizer):
    """
    Electron density represented in 2D
    """ 
    
    def displayDensity(self, densityValues):
        extent = [self.xLim[0], self.xLim[1], 
                self.yLim[0], self.yLim[1],
                log10(self.vmax), log10(self.vmin)]
        densityPlot = mlab.surf(densityValues[:,:,0],opacity=0.8,
                           extent=extent)
                           
        mlab.outline()
        mlab.axes()
        return densityPlot.mlab_source       

"*****************************************************************************"
def define_command_line_options(parser=None):
    if parser is None:
        import argparse
        parser = argparse.ArgumentParser()

    parser.add_argument('--showCores', action='store_true', default=1,
                        help='show cores')  
                        
    parser.add_argument('--showDensity', action='store_true', default=True,
                        help='show charge density')  
                        
    parser.add_argument(
        '--bgcolor', type=str, default= "b" , help='background coclor')
                        
    parser.add_argument('--animate', action='store_true', default=1,
                        help='make animation')
    
    parser.add_argument('--saveFigure', action='store_true', default=0,
                        help='save figures for movie')
                        
    parser.add_argument(
        '--vizType', type=int, default= 1 , 
        help='visual representation of electron density')
                             
    return parser
"*****************************************************************************"      
def main():
    # Read input from the command line
    parser = define_command_line_options()
    args = parser.parse_args()   
    
    #Read files
    geometricData, atomList, densityData = readFiles()
    
    # visual representation:
    if args.vizType ==1:
        vizType = volumeRepresentation
    elif args.vizType ==2:
        vizType = contourRepresentation
    elif args.vizType ==3:
        vizType = slicer1Representation
    elif args.vizType ==4:
        vizType = slicer2Representation
    elif args.vizType ==5:
        vizType = surf2DRepresentation
    else:
        sys.exit("unknwon visual type")
    
     
    cubeViz = vizType(geometricData, atomList, densityData,
                             args.bgcolor, args.showCores,args.showDensity,
                             args.animate, args.saveFigure)   
                                             
    cubeViz.vizualize()

"*****************************************************************************"
if __name__ == '__main__':
    main()
    







