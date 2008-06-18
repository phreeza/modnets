if [ "$(uname -s)" == 'Darwin' ];
then
g95 -o main nrutil.f90 nrtype.f90 indexx.f90 nr.f90  main.f90;
fi;
