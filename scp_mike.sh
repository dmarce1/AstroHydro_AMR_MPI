cd Release
rm *.silo.gz
rm *.dat
make clean
cd ..
cd Debug
rm *.silo.gz
rm *.dat
make clean
cd ..
cd ..
rm amr_mpi.tar.gz
tar -cvf "amr_mpi.tar" amr_mpi --exclude *.silo --exclude *.bin --exclude *.dat
gzip -f amr_mpi.tar
mv amr_mpi.tar.gz ./amr_mpi
cd amr_mpi
scp amr_mpi.tar.gz mike.hpc.lsu.edu:
ssh mike.hpc.lsu.edu "gunzip -f amr_mpi.tar.gz"
ssh mike.hpc.lsu.edu "rm -r ./amr_mpi/*"
ssh mike.hpc.lsu.edu "tar -xvf amr_mpi.tar" 
rm amr_mpi.tar.gz

