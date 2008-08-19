if [ "$(uname -s)" == 'Darwin' ];
then
g95 -O1 -o main nrutil.f90 gaussj.f90 jacobi.f90 nrtype.f90 indexx.f90 nr.f90  main.f90;
fi;
