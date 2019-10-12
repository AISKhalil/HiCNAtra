"""
Saves mirnylab.h5dict.h5dict file to a HDF5 file with different datsets.

Usage
-----

python h5dictToHDF5.py in_hdf5_file out_hdf5_file


Each key of the h5dict dataset is converted to a key in HDF5 file.
"""


from mirnylib.h5dict import h5dict
from mirnylib.numutils import generalizedDtype
import os
import sys
import numpy
import h5py
from scipy.io import savemat






def convertFile(filename,outFilename):
    
	if not os.path.exists(filename):
	    raise IOError("File not found: %s" % filename)

	outDict = h5py.File(outFilename, mode = 'w')
	mydict  = h5dict(filename, 'r')

	selectedKeys = ['chrms1', 'chrms2', 'cuts1', 'cuts2', 'rfragIdxs1', 'rfragIdxs2', 'strands1', 'strands2','rsites1','rsites2']
	#for i in list(mydict.keys()):
	for i in list(selectedKeys):    
		keyData = mydict[i]
		if issubclass(type(keyData), numpy.ndarray):
			print(("saving numpy array", i, "to", outFilename))
			if issubclass(type(keyData[0]), numpy.bool_):
				binaryStrands = numpy.zeros(len(keyData),dtype = numpy.int8)
				indices = numpy.where(keyData == True)
				binaryStrands[indices[0]] = 1 
				outDict.create_dataset(i, data = binaryStrands) 
			else:
				outDict.create_dataset(i, data = keyData)
			continue

	txtSavefile = i+'.txt'
	if type(keyData) == str:
	    datarepr = keyData
	else:
	    datarepr = repr(keyData)
	print(("saving data", i, "to", txtSavefile))
	with open(txtSavefile, 'w') as f:
	    f.write(datarepr)
	    
	outDict.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage : python h5dictToReadableHDF5.py in_h5dict_file out_hdf5_file")
        print("Converts h5dict file to a hdf5 file")
        print("Keys ('chrms1', 'chrms2', 'cuts1', 'cuts2', 'rfragIdxs1', 'rfragIdxs2', 'strands1', 'strands2','rsites1','rsites2') are kept ib the out_hdf5_file for HiCNAtra usage")
        exit()
        
    convertFile(sys.argv[1], sys.argv[2])




