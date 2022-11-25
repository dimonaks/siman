@echo off
setlocal 
::set path=%path%;D:\Elements\Root\Software\mingw64v5.2\mingw64\bin
g++ -c -std=c++11 -O3 findpores.cpp -o findpores.o
g++ -shared -fPIC -Wl,-soname,libfindpores.so -o libfindpores.so findpores.o
endlocal

