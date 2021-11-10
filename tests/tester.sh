#!/bin/bash
#executing this tests the current build of ARTEMIS by running the examples and comparing to expected results

test_dir=$(pwd)
home_dir=$(pwd|sed 's_[/]tests[/]*$__')
ARTEMIS=$home_dir/bin/artemis
echo $home_dir
echo $ARTEMIS
gfortran --version
## NED, FOR HOME COMPUTER, COMPILE ARTEMIS USING:
##make FC=gfortran-mp-10

test_dirs=$(ls -d */)
echo $test_dirs
for dir in $test_dirs; do
    cd $test_dir/$dir
    echo $(pwd)
    $ARTEMIS -f param.in
done

