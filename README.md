##### [UPDATING]: add SALEcVtsSlice
direct export slice of 3d output to vts files
usage:

    $./SALEcVtsReader ./plane.plot
example of input
```shell
[SALEc]
    input = /home/huacheng/Documents/Github/data/pdata/SALEc_20.inp
    step  = [1,1]

[Plane]
    data = /home/huacheng/Documents/Github/data/pdata/Al1100Test
    step = [range,1,2,1]
    output = Slice
   
    number = 2
    name = [xoz,xoy]
    nx = [0.0, 0.0]
    ny = [1.0, 0.0]
    nz = [0.0, 1.0]
    d  = [0.0, 0.0]
```

##### [UPDATING]: add VtpTracer.c (30/07/2024)
analysis 2d crater tracers from vtp data
usage:
    $VtpTracer -[option] [option args]
    -d [data directory]
        set the datafile directoty, DEFAULT is current path
    -f [input file name]
        set the *.inp file name, DEFAULT is sale2d.inp
    -v 
        convert the vtp tracers info to vts data
    -n 
        the max step

##### [UPDATING]: add Vts2dProfile (26/04/2024)
extract 2d crater profile from vts/vtm data
usage:

    -d [data directory] 
        set the datafile directoty, DEFAULT is current path
    -f [input file name]
        set the *.inp file name, DEFAULT is sale2d.inp
    -b 
        write the profile to binary format
    -o [output prefix]
        set the name of output file
    -e [Integer]
        set the step that loaded for profile


#### SALEcVtsReader
This is a simple tool to read vtm/vts data generated by SALEc and output a slice 
plane data. This tool is much faster than pvbatch in ParView.
The vtm/vts in SALEc is similar with citcoms and this tool can also applied to citcoms's
binary output.

##### compile
```shell
cd SALEcVtsReader-dir
mkdir build 
cd build
cmake .. && make
```
The default compile flag is -O0 and it can be edited if necessary.

#### run
SALEcVtsRead need a input file to specify the vts/vtm data dir. 
The executable file is SALEcVtsReader in build dir after compiled.
SALEcVtsReader requires a single argument. 
```shell
./SALEcVtsReader ./plane.plot
```
./plane.plot can be replaced by other input file. 
The following is an example of input.
```shell
# example input files for SALEcVtsReader by huachengli
[SALEc]
    input = /home/huacheng/Documents/Github/data/pdata/SALEc_20.inp
    step  = [1,1]

[Plane]
    data = /home/huacheng/Documents/Github/data/pdata/Al1100Test
    step = [range,1,2,1]
    output = Plane
   
    number = 2
    name = [Density,Density]
    nx = [0.0, 0.0]
    ny = [1.0, 0.0]
    nz = [0.0, 1.0]
    d  = [0.0, 0.0]
```

