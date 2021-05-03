import os
import sys
import glob
import vtk

def compressVTI(fileIn, fileOut):
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(fileIn)
    reader.Update()
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(fileOut)
    writer.SetInputData(reader.GetOutput())
    writer.Write()




if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("Usage: " + sys.argv[0] + " outputDirToProcess")
        sys.exit(0)

    outputDir = sys.argv[1]

    # Python 2.7 version, non-recursive, toplevel dir only
    # for file in os.listdir(outputDir):
    #     if file.endswith(".vti"):
    #         base = os.path.splitext(file)[0]
    #         if base[-2:] == "_c":
    #             continue

    #         inpFile = os.path.join(outputDir, file)
    #         outpFile = os.path.join(outputDir, base+"_c.vti")

    #         if(os.path.isfile(outpFile)):
    #             print("File " + inpFile + " already compressed, skipping....")
    #             continue

    #         print(inpFile + " --> " + outpFile)
    #         compressVTI(inpFile, outpFile)
            
    # Python 3.5+ version, fully recursive
    for filename in glob.iglob(outputDir + '**/*.vti', recursive=True):
        base = os.path.splitext(filename)[0]
        if base[-2:] == "_c":
            continue

        inpFile = filename
        outpFile = base+"_c.vti"

        if(os.path.isfile(outpFile)):
            print("File " + inpFile + " already compressed, skipping....")
            continue

        print(inpFile + " --> " + outpFile)
        compressVTI(inpFile, outpFile)
