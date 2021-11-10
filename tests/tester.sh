#!/bin/bash
#executing this tests the current build of ARTEMIS by running the examples and comparing to expected results

test_dir=$(pwd)
home_dir=$(pwd|sed 's_[/]tests[/]*$__')
ARTEMIS=$home_dir/bin/artemis
echo $home_dir
echo $ARTEMIS

test_dirs=$(ls -d */)
echo $test_dirs
for dir in $test_dirs; do
    cd $test_dir/$dir
    echo $(pwd)
    #if [[ "$dir" != "identify_terminations/" ]]; then
    #	continue
    #fi
    rm diff.txt
    $ARTEMIS -f param.in >> "test.out"
    POSCARS=$(find D* -name "POSCAR*" -not -path "DCHECK/*" -not -path "POSCAR*")
    for POSCAR in $POSCARS; do
	path=$(echo ${POSCAR} | sed 's_^[A-Za-z]*/_DCHECK/_')
	diff $POSCAR $path  >> "diff.txt"
    done
    if [ ! -s "diff.txt" ]; then
	echo "This test is consistent"
    else
	echo "This test has differences with what is expected, please check diff.txt and resolve the issues"
    fi
    echo
    
done

