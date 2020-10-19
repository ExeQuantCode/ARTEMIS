#!/bin/bash
#executing this tars up the code directory.

scriptdir="$(readlink -f "$(dirname "${BASH_SOURCE[0]}")")/"
codedir=$(sed 's|tools/||' <<<"$scriptdir")
codeend=$(sed 's|/.*/\(.*\)/|\1|' <<<"$codedir")
tarname=$(sed 's|/$||;s|.*/D\?||' <<<"$codedir")".tar.gz"
tardir=$(sed 's|/$||;s|\(.*/\).*|\1|' <<<"$codedir")
tarpath=$tardir$tarname

#cd $codedir

#echo "tar -czvf $tarpath -C $tardir $codeend"
#tar --exclude-backups -czvf $tarpath -C $tardir $codeend #$scriptdir
tar --exclude="*#" --exclude="*.#*" --exclude="*~" --exclude="TODO" --exclude="OLD_*" --exclude="*.aux" --exclude="*.log" --exclude="*.out" --exclude="*.toc" --exclude="*.txt" --exclude="${codeend}/bin" --exclude="*.tar.gz" --exclude="${codeend}/src/archive" --exclude="${codeend}/obj" -czvf $tarpath -C $tardir $codeend #$scriptdir
mv $tarpath ../.
echo "made tar file $tarname"
