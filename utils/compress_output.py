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
        print("Usage: " + sys.argv[0] + " <directory to walk recursively>")
        sys.exit(0)

    outputDir = sys.argv[1]

    # Python 3.5+ is required
    for filename in glob.iglob(outputDir + '**/*.vti', recursive=True):
        base = os.path.splitext(filename)[0]
        if base[-2:] == "_c":
            continue

        inpFile = filename#os.path.join(outputDir, filename)
        outpFile = base+"_c.vti"#os.path.join(outputDir, base+"_c.vti")

        if(os.path.isfile(outpFile)):
            print("File " + inpFile + " already compressed, skipping....")
            continue

        print(inpFile + " --> " + outpFile)
        compressVTI(inpFile, outpFile)
