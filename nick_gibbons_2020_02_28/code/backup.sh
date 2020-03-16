#!/bin/bash

echo "Syncing WmuteSSD/..."
rsync -tarvPp --delete /home/qungibbo/ /media/qungibbo/Elements/WintermuteSSD --exclude "*.ssh" --exclude "programs/*" --exclude "sourcecode/us3d*" --exclude "sourcecode/paraview*" --exclude "sourcecode/cfcfd3/*" --exclude "sourcecode/hdf5-1.8.12/*" --exclude "sourcecode/openmpi-1.8.12/*" --exclude "*.vtk" --exclude "*.o" --exclude "*.mod" --exclude ".cache*" --exclude ".wine/*" --exclude "*.swp" --exclude "*.pvtu" --exclude "*.vtu" --exclude "*.npy"


echo "Syncing WmuteDATA/..."
rsync -tarvPp --delete /media/qungibbo/data/ /media/qungibbo/Elements/WintermuteDATA --exclude "binaries/*" --exclude "*.vtk" --exclude "*.o" --exclude "*.mod" --exclude "*.swp" --exclude "*.pvtu" --exclude "*.vtu" --exclude "*.npy"
